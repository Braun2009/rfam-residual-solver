#!/usr/bin/env python3
"""Build and screen a reaction-retained C3Mech candidate with Cantera.

The training partition chooses the reaction-score threshold. The holdout
partition is not opened until that threshold and the retained mechanism have
been frozen in memory. This is detailed-mechanism surrogate evidence only; it
does not make the candidate an RFAM runtime mechanism.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence


class CandidateError(RuntimeError):
    pass


@dataclass(frozen=True)
class Case:
    case_id: str
    temperature: float
    pressure: float
    equivalence_ratio: float
    fuel: str
    oxidizer: str


@dataclass(frozen=True)
class Tolerances:
    idt_relative: float
    inverse_temperature_time_relative: float
    heat_release_nrmse: float
    species_mass_fraction_absolute: float


@dataclass
class Trajectory:
    complete: bool
    failure: str
    ignition_delay: float
    samples: list[dict[str, Any]]
    reaction_scores: list[float]
    accepted_steps: int
    maximum_temperature: float


OBSERVABLE_SPECIES = ("C3H8", "O2", "CO", "CO2", "H2O", "N2")
TEMPERATURE_RISE_NODES = (0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0)


def reject_duplicate_keys(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise CandidateError("duplicate-json-key:" + key)
        result[key] = value
    return result


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=lambda token: (_ for _ in ()).throw(
                CandidateError("nonfinite-json:" + token)
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CandidateError("json-read:" + str(path) + ":" + str(error)) from error
    if not isinstance(value, dict):
        raise CandidateError("json-root-not-object:" + str(path))
    return value


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("ascii")


def finite_positive(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value)) and float(value) > 0.0


def load_cases(path: Path, expected_partition: str) -> list[Case]:
    document = load_json(path)
    if document.get("schema") != "rfam.r2-c3mech-case-partition/v1":
        raise CandidateError("case-schema:" + str(path))
    if document.get("partition") != expected_partition:
        raise CandidateError("case-partition:" + str(path))
    rows = document.get("cases")
    if not isinstance(rows, list) or not rows:
        raise CandidateError("case-list:" + str(path))
    result: list[Case] = []
    seen: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            raise CandidateError("case-row-type:" + str(path))
        case_id = row.get("id")
        if not isinstance(case_id, str) or not case_id or case_id in seen:
            raise CandidateError("case-id:" + str(path))
        seen.add(case_id)
        values = (
            row.get("temperature_k"),
            row.get("pressure_pa"),
            row.get("equivalence_ratio"),
        )
        if not all(finite_positive(value) for value in values):
            raise CandidateError("case-state:" + case_id)
        fuel = row.get("fuel")
        oxidizer = row.get("oxidizer")
        if not isinstance(fuel, str) or not fuel or not isinstance(oxidizer, str) or not oxidizer:
            raise CandidateError("case-composition:" + case_id)
        result.append(
            Case(case_id, float(values[0]), float(values[1]), float(values[2]), fuel, oxidizer)
        )
    return result


def reaction_species(reaction: Any) -> set[str]:
    names = set(reaction.reactants) | set(reaction.products)
    names.update(getattr(reaction, "orders", {}) or {})
    third_body = getattr(reaction, "third_body", None)
    if third_body is not None:
        names.update(getattr(third_body, "efficiencies", {}) or {})
    return names


def close_selected_reactions(
    reactions: Sequence[Any], selected_indices: Iterable[int]
) -> tuple[list[int], set[str]]:
    selected = {int(index) for index in selected_indices}
    if any(index < 0 or index >= len(reactions) for index in selected):
        raise CandidateError("reaction-index-out-of-range")

    duplicate_equations = {
        reactions[index].equation
        for index in selected
        if bool(getattr(reactions[index], "duplicate", False))
    }
    if duplicate_equations:
        for index, reaction in enumerate(reactions):
            if bool(getattr(reaction, "duplicate", False)) and reaction.equation in duplicate_equations:
                selected.add(index)

    species: set[str] = set(OBSERVABLE_SPECIES)
    for index in selected:
        species.update(reaction_species(reactions[index]))
    return sorted(selected), species


NONEXECUTING_METADATA_KEYS = frozenset(("note", "description"))


def normalized_input_data(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): normalized_input_data(item)
            for key, item in sorted(value.items())
            if str(key) not in NONEXECUTING_METADATA_KEYS
        }
    if isinstance(value, (list, tuple)):
        return [normalized_input_data(item) for item in value]
    if hasattr(value, "tolist"):
        return normalized_input_data(value.tolist())
    if isinstance(value, float):
        if not math.isfinite(value):
            raise CandidateError("nonfinite-cantera-input-data")
        return format(value, ".14g")
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return str(value)


def semantic_digest(items: Iterable[Any]) -> str:
    payload = [normalized_input_data(item.input_data) for item in items]
    return sha256_bytes(canonical_bytes(payload))


def canonicalize_cantera_yaml(path: Path) -> tuple[bytes, str]:
    try:
        import yaml
    except ImportError as error:
        raise CandidateError("pyyaml-unavailable") from error
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        raise CandidateError("generated-yaml-read:" + str(error)) from error
    if not isinstance(document, dict):
        raise CandidateError("generated-yaml-root-not-object")
    document["date"] = "1970-01-01T00:00:00Z"
    try:
        raw = yaml.safe_dump(
            document,
            allow_unicode=False,
            default_flow_style=False,
            sort_keys=False,
            width=120,
        ).encode("ascii")
        path.write_bytes(raw)
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        raise CandidateError("generated-yaml-canonicalize:" + str(error)) from error
    return raw, str(yaml.__version__)


def semantic_numeric_difference(left: Any, right: Any, path: str = "root") -> float:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        left_keys = {str(key) for key in left if str(key) not in NONEXECUTING_METADATA_KEYS}
        right_keys = {str(key) for key in right if str(key) not in NONEXECUTING_METADATA_KEYS}
        if left_keys != right_keys:
            raise CandidateError("semantic-key-mismatch:" + path)
        return max(
            (
                semantic_numeric_difference(left[key], right[key], path + "." + str(key))
                for key in left
                if str(key) not in NONEXECUTING_METADATA_KEYS
            ),
            default=0.0,
        )
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            raise CandidateError("semantic-list-size-mismatch:" + path)
        return max(
            (
                semantic_numeric_difference(first, second, path + "[" + str(index) + "]")
                for index, (first, second) in enumerate(zip(left, right))
            ),
            default=0.0,
        )
    if isinstance(left, (int, float)) and not isinstance(left, bool) and isinstance(right, (int, float)) and not isinstance(right, bool):
        first = float(left)
        second = float(right)
        if not math.isfinite(first) or not math.isfinite(second):
            raise CandidateError("semantic-nonfinite:" + path)
        return abs(first - second) / max(abs(first), abs(second), 1.0e-300)
    if left != right:
        raise CandidateError("semantic-value-mismatch:" + path)
    return 0.0


def compare_semantic_items(left: Sequence[Any], right: Sequence[Any], label: str) -> float:
    if len(left) != len(right):
        raise CandidateError(label + "-count-mismatch")
    maximum = 0.0
    for index, (first, second) in enumerate(zip(left, right)):
        maximum = max(
            maximum,
            semantic_numeric_difference(
                first.input_data, second.input_data, label + "[" + str(index) + "]"
            ),
        )
    if maximum > 5.0e-14:
        raise CandidateError(label + "-numeric-drift:" + format(maximum, ".17g"))
    return maximum


def interpolate_sample(left: dict[str, Any], right: dict[str, Any], target_temperature: float) -> dict[str, Any]:
    denominator = right["temperature"] - left["temperature"]
    if not finite_positive(denominator):
        fraction = 1.0
    else:
        fraction = min(1.0, max(0.0, (target_temperature - left["temperature"]) / denominator))

    def blend(first: float, second: float) -> float:
        return first + fraction * (second - first)

    return {
        "temperature_rise": target_temperature - left["initial_temperature"],
        "time": blend(left["time"], right["time"]),
        "temperature": target_temperature,
        "heat_release_rate": blend(left["heat_release_rate"], right["heat_release_rate"]),
        "mass_fractions": {
            name: blend(left["mass_fractions"][name], right["mass_fractions"][name])
            for name in OBSERVABLE_SPECIES
        },
    }


def snapshot(gas: Any, time: float, initial_temperature: float) -> dict[str, Any]:
    fractions = {
        name: float(gas[name].Y[0]) if name in gas.species_names else 0.0
        for name in OBSERVABLE_SPECIES
    }
    return {
        "time": float(time),
        "initial_temperature": initial_temperature,
        "temperature": float(gas.T),
        "heat_release_rate": float(gas.heat_release_rate),
        "mass_fractions": fractions,
    }


def run_trajectory(
    gas_factory: Callable[[], Any],
    case: Case,
    capture_reaction_scores: bool,
    maximum_time: float,
) -> Trajectory:
    import cantera as ct

    gas = gas_factory()
    missing = [name for name in ("C3H8", "O2") if name not in gas.species_names]
    if missing:
        return Trajectory(False, "missing-initial-species:" + ",".join(missing), math.nan, [], [], 0, math.nan)
    try:
        gas.TP = case.temperature, case.pressure
        gas.set_equivalence_ratio(case.equivalence_ratio, case.fuel, case.oxidizer)
        gas.TP = case.temperature, case.pressure
        reactor = ct.IdealGasConstPressureReactor(gas, energy="on")
        network = ct.ReactorNet([reactor])
        network.rtol = 1.0e-10
        network.atol = 1.0e-18
        network.max_time_step = min(1.0e-4, maximum_time / 20.0)
    except ct.CanteraError as error:
        return Trajectory(False, "reactor-setup:" + str(error), math.nan, [], [], 0, math.nan)

    initial_temperature = case.temperature
    target_temperature = initial_temperature + TEMPERATURE_RISE_NODES[-1]
    if target_temperature >= 1995.0:
        return Trajectory(False, "target-temperature-outside-frozen-thermo-margin", math.nan, [], [], 0, initial_temperature)

    reaction_scores = [0.0] * gas.n_reactions if capture_reaction_scores else []
    previous = snapshot(gas, 0.0, initial_temperature)
    samples = [interpolate_sample(previous, previous, initial_temperature)]
    next_node = 1
    accepted_steps = 0
    ignition_delay = math.nan
    try:
        while network.time < maximum_time and accepted_steps < 250000:
            network.step()
            accepted_steps += 1
            current = snapshot(reactor.thermo, network.time, initial_temperature)
            if capture_reaction_scores:
                rates = [abs(float(value)) for value in reactor.thermo.net_rates_of_progress]
                scale = max(rates, default=0.0)
                if scale > 0.0 and math.isfinite(scale):
                    for index, value in enumerate(rates):
                        if math.isfinite(value):
                            reaction_scores[index] = max(reaction_scores[index], value / scale)
            while next_node < len(TEMPERATURE_RISE_NODES):
                node_temperature = initial_temperature + TEMPERATURE_RISE_NODES[next_node]
                if current["temperature"] < node_temperature:
                    break
                samples.append(interpolate_sample(previous, current, node_temperature))
                if next_node == len(TEMPERATURE_RISE_NODES) - 1:
                    ignition_delay = samples[-1]["time"]
                next_node += 1
            if next_node == len(TEMPERATURE_RISE_NODES):
                break
            if not all(
                math.isfinite(value)
                for value in (network.time, current["temperature"], current["heat_release_rate"])
            ):
                return Trajectory(False, "nonfinite-trajectory", math.nan, samples, reaction_scores, accepted_steps, current["temperature"])
            previous = current
    except ct.CanteraError as error:
        return Trajectory(False, "cantera-integration:" + str(error), math.nan, samples, reaction_scores, accepted_steps, previous["temperature"])

    complete = len(samples) == len(TEMPERATURE_RISE_NODES) and finite_positive(ignition_delay)
    failure = "" if complete else "ignition-rise-not-reached"
    return Trajectory(complete, failure, ignition_delay, samples, reaction_scores, accepted_steps, previous["temperature"])


def case_metrics(reference: Trajectory, candidate: Trajectory, tolerances: Tolerances) -> dict[str, Any]:
    if not reference.complete:
        raise CandidateError("reference-incomplete:" + reference.failure)
    if not candidate.complete:
        return {"pass": False, "failure": candidate.failure}
    if len(reference.samples) != len(candidate.samples):
        return {"pass": False, "failure": "sample-count-mismatch"}
    idt_relative = abs(candidate.ignition_delay / reference.ignition_delay - 1.0)
    inverse_time = max(
        abs(left["time"] - right["time"]) / reference.ignition_delay
        for left, right in zip(reference.samples[1:], candidate.samples[1:])
    )
    qdot_scale = max(abs(sample["heat_release_rate"]) for sample in reference.samples)
    qdot_scale = max(qdot_scale, 1.0)
    qdot_nrmse = math.sqrt(
        math.fsum(
            ((left["heat_release_rate"] - right["heat_release_rate"]) / qdot_scale) ** 2
            for left, right in zip(reference.samples, candidate.samples)
        )
        / len(reference.samples)
    )
    species_errors = {
        name: max(
            abs(left["mass_fractions"][name] - right["mass_fractions"][name])
            for left, right in zip(reference.samples, candidate.samples)
        )
        for name in OBSERVABLE_SPECIES
    }
    species_maximum = max(species_errors.values())
    passed = (
        idt_relative <= tolerances.idt_relative
        and inverse_time <= tolerances.inverse_temperature_time_relative
        and qdot_nrmse <= tolerances.heat_release_nrmse
        and species_maximum <= tolerances.species_mass_fraction_absolute
    )
    return {
        "pass": passed,
        "idt_relative": idt_relative,
        "inverse_temperature_time_relative": inverse_time,
        "heat_release_nrmse": qdot_nrmse,
        "species_mass_fraction_absolute": species_errors,
        "species_mass_fraction_absolute_maximum": species_maximum,
        "reference_ignition_delay_s": reference.ignition_delay,
        "candidate_ignition_delay_s": candidate.ignition_delay,
        "reference_steps": reference.accepted_steps,
        "candidate_steps": candidate.accepted_steps,
    }


def validate_thresholds(values: Sequence[float]) -> tuple[float, ...]:
    if not values or any(not math.isfinite(value) or value < 0.0 or value > 1.0 for value in values):
        raise CandidateError("threshold-range")
    result = tuple(values)
    if any(result[index] <= result[index + 1] for index in range(len(result) - 1)):
        raise CandidateError("threshold-order-must-be-strictly-descending")
    if result[-1] != 0.0:
        raise CandidateError("thresholds-must-end-in-zero-fallback")
    return result


def build_candidate(source: Any, scores: Sequence[float], threshold: float) -> tuple[Any, list[int], list[str]]:
    import cantera as ct

    reactions = source.reactions()
    if len(scores) != len(reactions):
        raise CandidateError("reaction-score-count")
    initial = [index for index, value in enumerate(scores) if value >= threshold]
    selected_indices, required_species = close_selected_reactions(reactions, initial)
    required_species.update(("C3H8", "O2"))
    unknown = sorted(required_species - set(source.species_names))
    if unknown:
        raise CandidateError("candidate-unknown-species:" + ",".join(unknown))
    selected_names = [name for name in source.species_names if name in required_species]
    selected_species = [source.species(name) for name in selected_names]
    selected_reactions = [reactions[index] for index in selected_indices]
    try:
        candidate = ct.Solution(
            thermo=source.thermo_model,
            kinetics=source.kinetics_model,
            species=selected_species,
            reactions=selected_reactions,
            name="gas",
        )
    except ct.CanteraError as error:
        raise CandidateError("candidate-build:" + str(error)) from error
    return candidate, selected_indices, selected_names


def rate_type_counts(reactions: Sequence[Any]) -> dict[str, int]:
    result: dict[str, int] = {}
    for reaction in reactions:
        name = str(reaction.reaction_type)
        result[name] = result.get(name, 0) + 1
    return dict(sorted(result.items()))


def screen_cases(
    gas_factory: Callable[[], Any],
    cases: Sequence[Case],
    references: Mapping[str, Trajectory],
    tolerances: Tolerances,
    maximum_time: float,
) -> tuple[bool, dict[str, Any]]:
    rows: dict[str, Any] = {}
    passed = True
    for case in cases:
        candidate = run_trajectory(gas_factory, case, False, maximum_time)
        metrics = case_metrics(references[case.case_id], candidate, tolerances)
        rows[case.case_id] = metrics
        passed = passed and bool(metrics["pass"])
    return passed, rows


def parse_thresholds(text: str) -> tuple[float, ...]:
    try:
        return validate_thresholds(tuple(float(token) for token in text.split(",")))
    except ValueError as error:
        raise CandidateError("threshold-number") from error


def fresh_output(path: Path, label: str) -> Path:
    absolute = Path(os.path.abspath(os.fspath(path)))
    if not absolute.parent.is_dir():
        raise CandidateError(label + "-parent-not-directory")
    if os.path.lexists(absolute):
        raise CandidateError(label + "-already-exists")
    return absolute


def execute(arguments: argparse.Namespace) -> tuple[bytes, bytes]:
    import cantera as ct

    mechanism = arguments.mechanism.resolve(strict=True)
    training_path = arguments.training.resolve(strict=True)
    holdout_path = arguments.holdout.resolve(strict=True)
    training = load_cases(training_path, "protocol-training")
    thresholds = parse_thresholds(arguments.thresholds)
    tolerances = Tolerances(
        arguments.idt_relative_tolerance,
        arguments.inverse_temperature_time_tolerance,
        arguments.heat_release_nrmse_tolerance,
        arguments.species_absolute_tolerance,
    )
    if not all(finite_positive(value) for value in tolerances.__dict__.values()):
        raise CandidateError("tolerance-range")
    if not finite_positive(arguments.maximum_time):
        raise CandidateError("maximum-time")

    source = ct.Solution(str(mechanism))
    for name in ("C3H8", "O2", *OBSERVABLE_SPECIES):
        if name not in source.species_names:
            raise CandidateError("source-missing-species:" + name)

    training_references: dict[str, Trajectory] = {}
    scores = [0.0] * source.n_reactions
    for case in training:
        trajectory = run_trajectory(
            lambda: ct.Solution(str(mechanism)), case, True, arguments.maximum_time
        )
        if not trajectory.complete:
            raise CandidateError("training-reference-incomplete:" + case.case_id + ":" + trajectory.failure)
        training_references[case.case_id] = trajectory
        scores = [max(left, right) for left, right in zip(scores, trajectory.reaction_scores)]

    search: list[dict[str, Any]] = []
    selected: Optional[tuple[Any, list[int], list[str], float, dict[str, Any]]] = None
    for threshold in thresholds:
        candidate, indices, names = build_candidate(source, scores, threshold)
        passed, metrics = screen_cases(
            lambda candidate=candidate: ct.Solution(
                thermo=candidate.thermo_model,
                kinetics=candidate.kinetics_model,
                species=candidate.species(),
                reactions=candidate.reactions(),
                name="gas",
            ),
            training,
            training_references,
            tolerances,
            arguments.maximum_time,
        )
        row = {
            "threshold": threshold,
            "species_count": candidate.n_species,
            "reaction_count": candidate.n_reactions,
            "pass": passed,
            "case_metrics": metrics,
        }
        search.append(row)
        if passed:
            selected = (candidate, indices, names, threshold, metrics)
            break
    if selected is None:
        raise CandidateError("no-training-candidate")

    candidate, selected_indices, selected_names, threshold, training_metrics = selected
    with tempfile.TemporaryDirectory(prefix="rfam-c3mech-retained-", dir=str(arguments.mechanism_output.parent)) as temporary:
        staged_mechanism = Path(temporary) / "candidate.yaml"
        candidate.write_yaml(staged_mechanism, precision=17, skip_user_defined=False)
        candidate_raw, yaml_version = canonicalize_cantera_yaml(staged_mechanism)
        reloaded = ct.Solution(str(staged_mechanism))
        source_reactions = [source.reaction(index) for index in selected_indices]
        source_species = [source.species(name) for name in selected_names]
        reloaded_reactions = reloaded.reactions()
        reloaded_species = reloaded.species()
        reaction_semantics = semantic_digest(source_reactions)
        species_semantics = semantic_digest(source_species)
        reaction_roundtrip_error = compare_semantic_items(
            source_reactions, reloaded_reactions, "reloaded-reaction-semantics"
        )
        species_roundtrip_error = compare_semantic_items(
            source_species, reloaded_species, "reloaded-species-thermo"
        )
        if semantic_digest(reloaded_reactions) != reaction_semantics:
            raise CandidateError("reloaded-reaction-semantic-digest-mismatch")
        if semantic_digest(reloaded_species) != species_semantics:
            raise CandidateError("reloaded-species-semantic-digest-mismatch")

        # The holdout authority is opened only after selection and a staged,
        # byte-stable candidate have been created and semantically checked.
        holdout = load_cases(holdout_path, "protocol-holdout")
        holdout_references: dict[str, Trajectory] = {}
        for case in holdout:
            reference = run_trajectory(
                lambda: ct.Solution(str(mechanism)), case, False, arguments.maximum_time
            )
            if not reference.complete:
                raise CandidateError("holdout-reference-incomplete:" + case.case_id + ":" + reference.failure)
            holdout_references[case.case_id] = reference
        holdout_pass, holdout_metrics = screen_cases(
            lambda: ct.Solution(str(staged_mechanism)),
            holdout,
            holdout_references,
            tolerances,
            arguments.maximum_time,
        )

        selected_reactions = [source.reaction(index) for index in selected_indices]
        result = {
            "schema": "rfam.r2-c3mech-retained-candidate/v1",
            "status": "internal-detailed-mechanism-surrogate-screen",
            "source": {
                "mechanism_path": str(mechanism),
                "mechanism_sha256": sha256_file(mechanism),
                "mechanism_bytes": mechanism.stat().st_size,
                "training_path": str(training_path),
                "training_sha256": sha256_file(training_path),
                "holdout_path": str(holdout_path),
                "holdout_sha256": sha256_file(holdout_path),
                "cantera_version": ct.__version__,
                "pyyaml_version": yaml_version,
                "python_version": sys.version.split()[0],
            },
            "protocol": {
                "reactor": "IdealGasConstPressureReactor",
                "energy": "on",
                "ignition_definition": "linear-interpolated-T-minus-T0-equals-300-K",
                "temperature_rise_nodes_k": list(TEMPERATURE_RISE_NODES),
                "maximum_time_s": arguments.maximum_time,
                "rtol": 1.0e-10,
                "atol": 1.0e-18,
                "tolerances": tolerances.__dict__,
                "thresholds": list(thresholds),
                "selection_partition": "training-only",
                "holdout_opened_after_selection": True,
            },
            "source_dimensions": {
                "species": source.n_species,
                "reactions": source.n_reactions,
                "reaction_types": rate_type_counts(source.reactions()),
            },
            "candidate": {
                "selected_threshold": threshold,
                "species_count": candidate.n_species,
                "reaction_count": candidate.n_reactions,
                "species_names": selected_names,
                "source_reaction_indices": selected_indices,
                "source_reaction_index_sha256": sha256_bytes(canonical_bytes(selected_indices)),
                "reaction_types": rate_type_counts(selected_reactions),
                "reversible_reactions": sum(bool(reaction.reversible) for reaction in selected_reactions),
                "duplicate_reactions": sum(bool(reaction.duplicate) for reaction in selected_reactions),
                "third_body_reactions": sum(getattr(reaction, "third_body", None) is not None for reaction in selected_reactions),
                "reaction_semantics_sha256": reaction_semantics,
                "species_thermo_semantics_sha256": species_semantics,
                "reaction_roundtrip_maximum_relative_numeric_difference": reaction_roundtrip_error,
                "species_roundtrip_maximum_relative_numeric_difference": species_roundtrip_error,
                "yaml_sha256": sha256_bytes(candidate_raw),
                "yaml_bytes": len(candidate_raw),
                "training_pass": True,
                "holdout_pass": holdout_pass,
            },
            "threshold_search": search,
            "training_metrics": training_metrics,
            "holdout_metrics": holdout_metrics,
            "claims": {
                "cantera_reaction_objects_retained_without_rate_law_simplification": True,
                "yaml_reload_semantics_verified": True,
                "mechanism_selected_for_rfam_runtime": False,
                "experimental_validation_complete": False,
                "phase_aware_chemistry_enabled": False,
                "production_ready": False,
            },
        }
        return candidate_raw, canonical_bytes(result)


def publish_pair(mechanism_output: Path, mechanism_raw: bytes, result_output: Path, result_raw: bytes) -> None:
    created: list[Path] = []
    with tempfile.TemporaryDirectory(prefix="rfam-retained-publish-", dir=str(mechanism_output.parent)) as temporary:
        staged_mechanism = Path(temporary) / "candidate.yaml"
        staged_result = Path(temporary) / "result.json"
        staged_mechanism.write_bytes(mechanism_raw)
        staged_result.write_bytes(result_raw)
        try:
            os.link(staged_mechanism, mechanism_output)
            created.append(mechanism_output)
            os.link(staged_result, result_output)
            created.append(result_output)
        except BaseException:
            for path in reversed(created):
                try:
                    path.unlink()
                except OSError:
                    pass
            raise


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--holdout", type=Path, required=True)
    parser.add_argument("--mechanism-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    parser.add_argument(
        "--thresholds",
        default="0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001,0",
    )
    parser.add_argument("--maximum-time", type=float, default=0.005)
    parser.add_argument("--idt-relative-tolerance", type=float, default=0.20)
    parser.add_argument("--inverse-temperature-time-tolerance", type=float, default=0.25)
    parser.add_argument("--heat-release-nrmse-tolerance", type=float, default=0.25)
    parser.add_argument("--species-absolute-tolerance", type=float, default=0.02)
    arguments = parser.parse_args(argv)
    try:
        mechanism_output = fresh_output(arguments.mechanism_output, "mechanism-output")
        result_output = fresh_output(arguments.result_output, "result-output")
        if mechanism_output == result_output:
            raise CandidateError("output-path-collision")
        arguments.mechanism_output = mechanism_output
        arguments.result_output = result_output
        mechanism_raw, result_raw = execute(arguments)
        publish_pair(mechanism_output, mechanism_raw, result_output, result_raw)
        result = json.loads(result_raw)
        print(
            "R2_C3MECH_RETAINED_CANDIDATE_BLOCKED"
            f" species={result['candidate']['species_count']}"
            f" reactions={result['candidate']['reaction_count']}"
            f" training_pass={int(result['candidate']['training_pass'])}"
            f" holdout_pass={int(result['candidate']['holdout_pass'])}"
            " production=0"
        )
        return 0
    except Exception as error:
        reason = str(error).replace("\r", "\\r").replace("\n", "\\n")
        print("R2_C3MECH_RETAINED_CANDIDATE_FAIL " + reason, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
