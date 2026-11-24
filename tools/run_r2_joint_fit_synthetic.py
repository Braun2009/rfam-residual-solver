#!/usr/bin/env python3
"""Execute the frozen R2 joint-fit protocol against a synthetic affine runner.

This is deliberately not a chemistry fit.  It binds the accepted v3 objective
authorities and a repository-authored synthetic execution contract, then
exercises optimizer ordering, accounting, residual assembly, and the declared
local identifiability calculations without loading Cantera or RFAM chemistry.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple


BLOCKED_MARKER = "R2_JOINT_FIT_SYNTHETIC_EXECUTION_BLOCKED"
FAIL_MARKER = "R2_JOINT_FIT_SYNTHETIC_EXECUTION_FAIL"

EXECUTION_AUTHORITY = {
    "path": "tests/reference/r2_joint_fit_synthetic_execution_contract_v1.json",
    "sha256": "4c3db0f2711a00d93d8bdff4b08a72df41107db1a930cd072cc10fa7d4b0e268",
    "byte_count": 6028,
}

FROZEN_V3_AUTHORITIES = {
    "objective_contract": {
        "path": "tests/reference/r2_joint_fit_objective_contract_v1.json",
        "sha256": "4e15ed2c7cca709b5bd54f04faa6cd2db50347e31dacc9676213b3fe375679e7",
        "byte_count": 15830,
    },
    "synthetic_fixture": {
        "path": "tests/reference/r2_joint_fit_synthetic_fixture_v1.json",
        "sha256": "6cef95b7fd3d56eee4a009214a5205d539c69681712e3088be6b7418a1b2199c",
        "byte_count": 6652,
    },
    "objective_validator": {
        "path": "tools/validate_r2_joint_fit_objective.py",
        "sha256": "7e04bc983e3f3312b0fed1afb3433ca83fcdb4d93e88fd0c73fa526d1ce5bb50",
        "byte_count": 35024,
    },
    "objective_test": {
        "path": "tests/test_r2_joint_fit_objective_contract.py",
        "sha256": "ec0005ccba7715dd00d6eddeaa9b5297ef9a2c3bea9a446c87cfff37fab5c3a8",
        "byte_count": 13622,
    },
    "objective_documentation": {
        "path": "docs/r2_joint_fit_objective_contract.md",
        "sha256": "cece29648591d19f498608dee9bf0f3510d82199f710ec506b56a4cf77197538",
        "byte_count": 6529,
    },
}

OBSERVABLE_ORDER = (
    "ignition_delay", "qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O"
)
HISTORY_OBSERVABLES = ("qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O")

EXECUTION_APPROVALS = {
    "detailed_mechanism_selected_for_fit": False,
    "experimental_training_data_selected": False,
    "experimental_holdout_data_selected": False,
    "physical_right_to_fit_resolved": False,
    "physical_right_to_redistribute_resolved": False,
    "parameter_publication_approved": False,
    "physical_validation_approved": False,
    "hem_chemistry_enabled": False,
    "fluent_execution_approved": False,
    "reacting_hot_start_approved": False,
    "production_ready": False,
}

EXPECTED_RESPONSE_ENTRIES = (
    ("R1_A_SI", "SYNTHETIC_TRAIN_001", "qdot_w_m3", 1, 100000000.0),
    ("R1_Ea_J_per_mol", "SYNTHETIC_TRAIN_001", "qdot_w_m3", 2, 100000000.0),
    ("R1_order_C3H8", "SYNTHETIC_TRAIN_001", "qdot_w_m3", 3, 100000000.0),
    ("R1_order_O2", "SYNTHETIC_TRAIN_001", "temperature_k", 1, 20.0),
    ("R2_A_SI", "SYNTHETIC_TRAIN_001", "temperature_k", 2, 20.0),
    ("R2_Ea_J_per_mol", "SYNTHETIC_TRAIN_001", "temperature_k", 3, 20.0),
)

EXECUTION_PROTOCOL = {
    "optimizer_source": "frozen-v3-objective-contract.optimizer_protocol",
    "training_partition_source": "frozen-v3-objective-contract.partitions.training_case_ids",
    "holdout_partition_source": "frozen-v3-objective-contract.partitions.holdout_case_ids",
    "holdout_runner_calls_during_selection": 0,
    "holdout_runner_calls_after_selection": 1,
    "residual_and_identifiability_source": "frozen-v3-objective-contract.identifiability_uncertainty",
    "sensitivity_runner_cache": "separate-from-optimizer-no-cross-point-cache",
    "linear_algebra": "deterministic-symmetric-gram-jacobi-v1",
    "jacobi_pair_order": "cyclic-p-ascending-then-q-ascending",
    "jacobi_relative_offdiagonal_tolerance": 1e-15,
    "jacobi_maximum_sweeps": 100,
    "eigenpair_order": "descending-eigenvalue-then-original-column-index",
    "result_float_serialization": "finite-binary64-rounded-to-15-significant-decimal-digits",
    "canonical_json": "utf8-indent-2-preserve-schema-order-trailing-newline",
}

EXECUTION_CLAIMS = {
    "synthetic_optimizer_may_execute": True,
    "physical_optimizer_executed": False,
    "physical_fit_executed": False,
    "parameters_publishable": False,
    "synthetic_identifiability_is_physical_evidence": False,
    "physical_uncertainty_claimed": False,
    "runner_equivalent_to_rfam_or_cantera": False,
}

EXECUTION_STOP_CONDITIONS = [
    "any frozen-v3 authority path hash or byte count mismatch",
    "execution contract duplicate key unknown field nonfinite value or byte drift",
    "any holdout runner call before global training selection is frozen",
    "any nonfinite invalid or out-of-bound parameter runner residual sensitivity or covariance value",
    "optimizer protocol cache budget order bound or tie-break drift",
    "rank condition covariance or correlation protocol drift",
    "synthetic runner represented as mechanism experiment physical fit validation or production evidence",
    "attempt to enable HEM Fluent hot-start parameter publication or production capability",
]


class SyntheticExecutionError(RuntimeError):
    """A fail-closed contract, runner, or numerical protocol violation."""


def _strict_pairs(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise SyntheticExecutionError("duplicate-json-key:" + key)
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise SyntheticExecutionError("nonfinite-json-number:" + value)


def _load_json_bytes(raw: bytes, label: str) -> Dict[str, Any]:
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise SyntheticExecutionError("invalid-utf8:" + label) from exc
    try:
        value = json.loads(
            text,
            object_pairs_hook=_strict_pairs,
            parse_constant=_reject_constant,
        )
    except SyntheticExecutionError:
        raise
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise SyntheticExecutionError("invalid-json:" + label) from exc
    if not isinstance(value, dict):
        raise SyntheticExecutionError("json-root-not-object:" + label)
    return value


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _resolve_under(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise SyntheticExecutionError("invalid-relative-path:" + label)
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise SyntheticExecutionError("path-escapes-root:" + label) from exc
    return candidate


def _read_bound(root: Path, authority: Dict[str, Any], label: str) -> Tuple[Path, bytes]:
    path = _resolve_under(root, authority["path"], label + ".path")
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise SyntheticExecutionError("artifact-read-failed:" + label) from exc
    if len(raw) != authority["byte_count"]:
        raise SyntheticExecutionError("artifact-byte-count-mismatch:" + label)
    if _sha256(raw) != authority["sha256"]:
        raise SyntheticExecutionError("artifact-sha256-mismatch:" + label)
    return path, raw


def _require_object(value: Any, label: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise SyntheticExecutionError("expected-object:" + label)
    return value


def _require_list(value: Any, label: str) -> List[Any]:
    if not isinstance(value, list):
        raise SyntheticExecutionError("expected-list:" + label)
    return value


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SyntheticExecutionError("expected-number:" + label)
    result = float(value)
    if not math.isfinite(result):
        raise SyntheticExecutionError("nonfinite-number:" + label)
    return result


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise SyntheticExecutionError("contract-mismatch:" + label)
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise SyntheticExecutionError("contract-mismatch:" + label)
        for key in expected:
            _exact(actual[key], expected[key], label + "." + key)
        return
    if isinstance(expected, (list, tuple)):
        if len(actual) != len(expected):
            raise SyntheticExecutionError("contract-mismatch:" + label)
        for index, (actual_item, expected_item) in enumerate(zip(actual, expected)):
            _exact(actual_item, expected_item, label + "[" + str(index) + "]")
        return
    if actual != expected:
        raise SyntheticExecutionError("contract-mismatch:" + label)


def _exact_keys(value: Dict[str, Any], expected: Iterable[str], label: str) -> None:
    expected_set = set(expected)
    actual_set = set(value)
    missing = sorted(expected_set - actual_set)
    unknown = sorted(actual_set - expected_set)
    if missing:
        raise SyntheticExecutionError("missing-keys:" + label + ":" + ",".join(missing))
    if unknown:
        raise SyntheticExecutionError("unknown-keys:" + label + ":" + ",".join(unknown))


@dataclass(frozen=True)
class ParameterSpec:
    name: str
    transform: str
    lower: float
    upper: float
    transformed_lower: float
    transformed_upper: float

    def transformed_from_u(self, u: float) -> float:
        return self.transformed_lower + u * (self.transformed_upper - self.transformed_lower)

    def physical_from_u(self, u: float) -> float:
        transformed = self.transformed_from_u(u)
        if self.transform == "log10":
            return 10.0 ** transformed
        return transformed

    def u_from_physical(self, value: float) -> float:
        transformed = math.log10(value) if self.transform == "log10" else value
        return (transformed - self.transformed_lower) / (
            self.transformed_upper - self.transformed_lower
        )


@dataclass
class ObjectiveEvaluation:
    valid: bool
    objective: float
    residuals: Tuple[float, ...]
    case_objectives: Dict[str, float]
    reason: str = ""

    @classmethod
    def invalid(cls, reason: str) -> "ObjectiveEvaluation":
        return cls(False, math.inf, tuple(), {}, reason)


@dataclass
class SearchConfig:
    initial_step: float
    contraction: float
    minimum_step: float
    maximum_evaluations: int
    coordinate_count: int


@dataclass
class SearchResult:
    valid: bool
    start_index: int
    start_kind: str
    start_u: Tuple[float, ...]
    selected_u: Tuple[float, ...]
    evaluation: ObjectiveEvaluation
    evaluation_count: int
    cache_hits: int
    out_of_bound_skips: int
    completed_sweeps: int
    contractions: int
    accepted_coordinates: int
    termination: str
    trace_event_count: int
    trace_sha256: str


def radical_inverse(index: int, base: int) -> float:
    if isinstance(index, bool) or not isinstance(index, int) or index <= 0:
        raise SyntheticExecutionError("invalid-halton-index")
    if isinstance(base, bool) or not isinstance(base, int) or base <= 1:
        raise SyntheticExecutionError("invalid-halton-base")
    result = 0.0
    factor = 1.0 / float(base)
    remaining = index
    while remaining:
        result += float(remaining % base) * factor
        remaining //= base
        factor /= float(base)
    if not (0.0 <= result < 1.0) or not math.isfinite(result):
        raise SyntheticExecutionError("invalid-halton-result")
    return result


def _point_key(point: Sequence[float]) -> Tuple[str, ...]:
    return tuple(float(value).hex() for value in point)


def _trace_digest(events: Sequence[Dict[str, Any]]) -> str:
    raw = json.dumps(events, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return _sha256(raw)


def _transformed_key(point: Sequence[float], specs: Sequence[ParameterSpec]) -> Tuple[float, ...]:
    return tuple(spec.transformed_from_u(value) for spec, value in zip(specs, point))


def run_pattern_search(
    start_index: int,
    start_kind: str,
    start_u: Sequence[float],
    specs: Sequence[ParameterSpec],
    config: SearchConfig,
    evaluate: Callable[[Tuple[float, ...]], ObjectiveEvaluation],
) -> SearchResult:
    """Run one frozen bounded coordinate-pattern-search start."""

    if len(start_u) != config.coordinate_count or len(specs) != config.coordinate_count:
        raise SyntheticExecutionError("search-coordinate-count-mismatch")
    start = tuple(_finite(value, "search.start") for value in start_u)
    if any(value < 0.0 or value > 1.0 for value in start):
        raise SyntheticExecutionError("search-start-out-of-bound")

    cache: Dict[Tuple[str, ...], ObjectiveEvaluation] = {}
    events: List[Dict[str, Any]] = []
    evaluation_count = 0
    cache_hits = 0
    out_of_bound_skips = 0
    completed_sweeps = 0
    contractions = 0
    accepted_coordinates = 0

    def cached_evaluate(point: Tuple[float, ...]) -> ObjectiveEvaluation:
        nonlocal evaluation_count, cache_hits
        key = _point_key(point)
        if key in cache:
            cache_hits += 1
            events.append({"event": "cache-hit", "point": list(key)})
            return cache[key]
        if evaluation_count >= config.maximum_evaluations:
            raise SyntheticExecutionError("evaluation-budget-overrun")
        value = evaluate(point)
        if not isinstance(value, ObjectiveEvaluation):
            raise SyntheticExecutionError("evaluator-return-type-invalid")
        if value.valid and (not math.isfinite(value.objective) or value.objective < 0.0):
            value = ObjectiveEvaluation.invalid("invalid-objective")
        cache[key] = value
        evaluation_count += 1
        events.append({
            "event": "evaluate",
            "point": list(key),
            "valid": value.valid,
            "objective": float(value.objective).hex() if value.valid else "invalid",
        })
        return value

    incumbent_u = start
    incumbent = cached_evaluate(incumbent_u)
    if not incumbent.valid:
        events.append({"event": "terminate", "reason": "invalid-start"})
        return SearchResult(
            False, start_index, start_kind, start, start, incumbent,
            evaluation_count, cache_hits, out_of_bound_skips, completed_sweeps,
            contractions, accepted_coordinates, "invalid-start", len(events),
            _trace_digest(events),
        )

    step = config.initial_step
    termination = ""
    while not termination:
        accepted_in_sweep = False
        sweep_complete = True
        for coordinate in range(config.coordinate_count):
            entry_u = incumbent_u
            candidates: List[Tuple[str, Tuple[float, ...]]] = []
            for direction, delta in (("positive", step), ("negative", -step)):
                coordinate_value = entry_u[coordinate] + delta
                if coordinate_value < 0.0 or coordinate_value > 1.0:
                    out_of_bound_skips += 1
                    events.append({
                        "event": "out-of-bound-skip",
                        "coordinate": coordinate,
                        "direction": direction,
                        "step": float(step).hex(),
                        "entry": list(_point_key(entry_u)),
                    })
                    continue
                candidate = list(entry_u)
                candidate[coordinate] = coordinate_value
                candidates.append((direction, tuple(candidate)))

            uncached = {
                _point_key(candidate)
                for _, candidate in candidates
                if _point_key(candidate) not in cache
            }
            remaining_budget = config.maximum_evaluations - evaluation_count
            if len(uncached) > remaining_budget:
                events.append({
                    "event": "budget-stop-before-coordinate",
                    "coordinate": coordinate,
                    "needed": len(uncached),
                    "remaining": remaining_budget,
                })
                termination = "evaluation-limit-before-coordinate"
                sweep_complete = False
                break

            valid_candidates: List[Tuple[Tuple[float, ...], ObjectiveEvaluation]] = []
            for direction, candidate in candidates:
                value = cached_evaluate(candidate)
                events.append({
                    "event": "candidate-result",
                    "coordinate": coordinate,
                    "direction": direction,
                    "point": list(_point_key(candidate)),
                    "valid": value.valid,
                })
                if value.valid:
                    valid_candidates.append((candidate, value))

            if valid_candidates:
                candidate_u, candidate_value = min(
                    valid_candidates,
                    key=lambda item: (item[1].objective, _transformed_key(item[0], specs)),
                )
                if candidate_value.objective < incumbent.objective:
                    incumbent_u = candidate_u
                    incumbent = candidate_value
                    accepted_coordinates += 1
                    accepted_in_sweep = True
                    events.append({
                        "event": "accept",
                        "coordinate": coordinate,
                        "point": list(_point_key(candidate_u)),
                        "objective": float(candidate_value.objective).hex(),
                    })

        if not sweep_complete:
            break
        completed_sweeps += 1
        events.append({
            "event": "sweep-complete",
            "accepted": accepted_in_sweep,
            "step": float(step).hex(),
        })
        if accepted_in_sweep:
            continue
        next_step = step * config.contraction
        if next_step < config.minimum_step:
            termination = "minimum-step-reached"
            events.append({"event": "terminate", "reason": termination})
            break
        step = next_step
        contractions += 1
        events.append({"event": "contract", "step": float(step).hex()})

    if not termination:
        raise SyntheticExecutionError("search-termination-missing")
    return SearchResult(
        True, start_index, start_kind, start, incumbent_u, incumbent,
        evaluation_count, cache_hits, out_of_bound_skips, completed_sweeps,
        contractions, accepted_coordinates, termination, len(events),
        _trace_digest(events),
    )


class SyntheticAffineRunner:
    """Fixture-backed affine response with explicit split-scoped call accounting."""

    def __init__(
        self,
        objective_contract: Dict[str, Any],
        fixture: Dict[str, Any],
        execution: Dict[str, Any],
        parameter_specs: Sequence[ParameterSpec],
    ) -> None:
        self.objective_contract = objective_contract
        self.fixture = fixture
        self.parameter_specs = tuple(parameter_specs)
        self.parameter_index = {spec.name: index for index, spec in enumerate(parameter_specs)}
        trial = _require_list(fixture.get("synthetic_trial_parameter_vector"), "fixture.trial")
        if len(trial) != len(parameter_specs):
            raise SyntheticExecutionError("fixture-trial-parameter-count-mismatch")
        anchor_u = []
        for spec, value in zip(parameter_specs, trial):
            physical = _finite(value, "fixture.trial." + spec.name)
            u = spec.u_from_physical(physical)
            if u < 0.0 or u > 1.0 or not math.isfinite(u):
                raise SyntheticExecutionError("fixture-trial-out-of-bound:" + spec.name)
            anchor_u.append(u)
        self.anchor_u = tuple(anchor_u)
        runner_contract = _require_object(execution.get("synthetic_runner"), "execution.runner")
        entries = _require_list(
            runner_contract.get("training_response_entries"), "execution.training_response_entries"
        )
        self.response_entries: List[Dict[str, Any]] = []
        seen = set()
        for entry_value in entries:
            entry = _require_object(entry_value, "execution.response-entry")
            _exact_keys(
                entry,
                {"parameter", "case_id", "observable", "sample_index", "raw_slope_per_normalized_coordinate"},
                "execution.response-entry",
            )
            parameter = entry["parameter"]
            if parameter not in self.parameter_index:
                raise SyntheticExecutionError("response-parameter-unknown")
            case_id = entry["case_id"]
            observable = entry["observable"]
            sample_index = entry["sample_index"]
            if isinstance(sample_index, bool) or not isinstance(sample_index, int) or sample_index < 0:
                raise SyntheticExecutionError("response-sample-index-invalid")
            slope = _finite(entry["raw_slope_per_normalized_coordinate"], "response.slope")
            key = (parameter, case_id, observable, sample_index)
            if key in seen:
                raise SyntheticExecutionError("duplicate-response-entry")
            seen.add(key)
            self.response_entries.append({
                "parameter": parameter,
                "case_id": case_id,
                "observable": observable,
                "sample_index": sample_index,
                "slope": slope,
            })
        self._validate_case_inputs()
        self.calls: List[Dict[str, Any]] = []

    def _validate_case_inputs(self) -> None:
        mechanism = _require_object(
            self.objective_contract.get("mechanism_adapter"), "contract.mechanism-adapter"
        )
        domain = _require_object(self.objective_contract.get("domain"), "contract.domain")
        composition = _require_object(
            domain.get("initial_composition"), "contract.initial-composition"
        )
        observable_basis = _require_object(
            self.objective_contract.get("observables_and_basis"), "contract.observables"
        )
        energy_basis = _require_object(
            observable_basis.get("energy_basis"), "contract.energy-basis"
        )
        self.species_order = tuple(
            _require_list(self.fixture.get("species_order"), "fixture.species-order")
        )
        if list(self.species_order) != mechanism.get("species_order"):
            raise SyntheticExecutionError("fixture-mechanism-species-order-mismatch")
        if list(self.species_order) != composition.get("species_order"):
            raise SyntheticExecutionError("fixture-composition-species-order-mismatch")
        parameter_order = _require_list(self.fixture.get("parameter_order"), "fixture.parameter-order")
        if parameter_order != [spec.name for spec in self.parameter_specs]:
            raise SyntheticExecutionError("fixture-parameter-order-mismatch")
        self.energy_basis_id = self.fixture.get("energy_basis_id")
        if self.energy_basis_id != energy_basis.get("id"):
            raise SyntheticExecutionError("fixture-energy-basis-mismatch")
        required_inputs = _require_list(
            mechanism.get("required_runner_inputs"), "contract.required-runner-inputs"
        )
        self.required_input_keys = tuple(required_inputs)
        expected_input_keys = (
            "ordered-parameter-vector",
            "initial-temperature-k",
            "initial-pressure-pa",
            "equivalence-ratio",
            "diluent-mass-fraction",
            "ordered-initial-mass-fraction-vector",
            "ordered-sample-times-s",
            "energy-basis-id",
        )
        if self.required_input_keys != expected_input_keys:
            raise SyntheticExecutionError("runner-required-input-order-mismatch")
        partitions = _require_object(
            self.objective_contract.get("partitions"), "contract.partitions"
        )
        expected_case_ids = list(partitions["training_case_ids"]) + list(
            partitions["holdout_case_ids"]
        )
        cases = _require_list(self.fixture.get("cases"), "fixture.cases")
        if [case.get("case_id") for case in cases if isinstance(case, dict)] != expected_case_ids:
            raise SyntheticExecutionError("fixture-case-order-mismatch")
        self.case_request_counts = {"training": 0, "holdout": 0}
        for case_value in cases:
            case = _require_object(case_value, "fixture.case")
            condition = _require_object(case.get("condition"), "fixture.condition")
            self._validate_case_condition(condition, case)
            split = case.get("split")
            if split not in self.case_request_counts:
                raise SyntheticExecutionError("fixture-case-split-invalid")
            self.case_request_counts[split] += 1

    def _validate_case_condition(self, condition: Dict[str, Any], case: Dict[str, Any]) -> None:
        domain = self.objective_contract["domain"]
        temperature = _finite(condition.get("initial_temperature_k"), "condition.temperature")
        pressure = _finite(condition.get("initial_pressure_pa"), "condition.pressure")
        equivalence_ratio = _finite(condition.get("equivalence_ratio"), "condition.phi")
        dilution = _finite(condition.get("diluent_mass_fraction"), "condition.dilution")
        for value, limits, label in (
            (temperature, domain["initial_temperature_k"], "temperature"),
            (pressure, domain["initial_pressure_pa"], "pressure"),
            (equivalence_ratio, domain["equivalence_ratio"], "equivalence-ratio"),
            (dilution, domain["diluent_mass_fraction"], "diluent-mass-fraction"),
        ):
            lower = _finite(limits[0], "domain." + label + ".lower")
            upper = _finite(limits[1], "domain." + label + ".upper")
            if value < lower or value > upper:
                raise SyntheticExecutionError("fixture-case-out-of-domain:" + label)
        mass_fractions = _require_list(
            condition.get("initial_mass_fractions"), "condition.mass-fractions"
        )
        if len(mass_fractions) != len(self.species_order):
            raise SyntheticExecutionError("fixture-initial-mass-fraction-count-mismatch")
        y = [_finite(value, "condition.mass-fraction") for value in mass_fractions]
        if any(value < 0.0 or value > 1.0 for value in y):
            raise SyntheticExecutionError("fixture-initial-mass-fraction-out-of-range")
        composition = self.objective_contract["domain"]["initial_composition"]
        tolerance = _finite(composition["sum_tolerance"], "composition.sum-tolerance")
        if abs(math.fsum(y) - 1.0) > tolerance:
            raise SyntheticExecutionError("fixture-initial-mass-fraction-sum-mismatch")
        n2_index = self.species_order.index(composition["diluent_species"])
        if abs(y[n2_index] - dilution) > tolerance:
            raise SyntheticExecutionError("fixture-diluent-metadata-mismatch")
        fuel_index = self.species_order.index("C3H8")
        oxygen_index = self.species_order.index("O2")
        molecular_weights = composition["molecular_weights_kg_per_kmol"]
        stoichiometric_mass_ratio = _finite(
            molecular_weights["C3H8"], "composition.fuel-molecular-weight"
        ) / (3.5 * _finite(molecular_weights["O2"], "composition.oxygen-molecular-weight"))
        derived_phi = (y[fuel_index] / y[oxygen_index]) / stoichiometric_mass_ratio
        metadata_tolerance = _finite(
            composition["metadata_consistency_absolute_tolerance"],
            "composition.metadata-tolerance",
        )
        if abs(derived_phi - equivalence_ratio) > metadata_tolerance:
            raise SyntheticExecutionError("fixture-equivalence-ratio-metadata-mismatch")
        times = _require_list(case.get("sample_times_s"), "fixture.sample-times")
        parsed_times = [_finite(value, "fixture.sample-time") for value in times]
        if not parsed_times or any(
            parsed_times[index] <= parsed_times[index - 1]
            for index in range(1, len(parsed_times))
        ):
            raise SyntheticExecutionError("fixture-sample-times-not-strictly-increasing")
        time_domain = self.objective_contract["domain"]["sample_time_s"]
        if parsed_times[0] < _finite(time_domain[0], "domain.time.lower") or parsed_times[-1] > _finite(
            time_domain[1], "domain.time.upper"
        ):
            raise SyntheticExecutionError("fixture-sample-time-out-of-domain")

    def _case_request(
        self,
        case: Dict[str, Any],
        physical_parameters: Sequence[float],
    ) -> Dict[str, Any]:
        condition = _require_object(case.get("condition"), "fixture.condition")
        request = {
            "ordered-parameter-vector": list(physical_parameters),
            "initial-temperature-k": condition.get("initial_temperature_k"),
            "initial-pressure-pa": condition.get("initial_pressure_pa"),
            "equivalence-ratio": condition.get("equivalence_ratio"),
            "diluent-mass-fraction": condition.get("diluent_mass_fraction"),
            "ordered-initial-mass-fraction-vector": copy.deepcopy(
                condition.get("initial_mass_fractions")
            ),
            "ordered-sample-times-s": copy.deepcopy(case.get("sample_times_s")),
            "energy-basis-id": self.energy_basis_id,
        }
        if tuple(request) != self.required_input_keys:
            raise SyntheticExecutionError("runner-request-field-order-mismatch")
        self._consume_case_request(request)
        return request

    def _consume_case_request(self, request: Dict[str, Any]) -> None:
        _exact_keys(request, self.required_input_keys, "runner.request")
        parameters = _require_list(
            request["ordered-parameter-vector"], "runner.request.parameters"
        )
        if len(parameters) != len(self.parameter_specs) or any(
            not math.isfinite(_finite(value, "runner.request.parameter")) for value in parameters
        ):
            raise SyntheticExecutionError("runner-request-parameter-vector-invalid")
        if request["energy-basis-id"] != self.energy_basis_id:
            raise SyntheticExecutionError("runner-request-energy-basis-mismatch")
        y = _require_list(
            request["ordered-initial-mass-fraction-vector"], "runner.request.mass-fractions"
        )
        times = _require_list(request["ordered-sample-times-s"], "runner.request.sample-times")
        if len(y) != len(self.species_order) or not times:
            raise SyntheticExecutionError("runner-request-vector-shape-invalid")
        for key in (
            "initial-temperature-k", "initial-pressure-pa", "equivalence-ratio",
            "diluent-mass-fraction",
        ):
            _finite(request[key], "runner.request." + key)
        for value in y:
            _finite(value, "runner.request.mass-fraction")
        for value in times:
            _finite(value, "runner.request.sample-time")

    def physical_parameters(self, u: Sequence[float]) -> Tuple[float, ...]:
        self._validate_u(u)
        return tuple(spec.physical_from_u(value) for spec, value in zip(self.parameter_specs, u))

    def transformed_parameters(self, u: Sequence[float]) -> Tuple[float, ...]:
        self._validate_u(u)
        return tuple(spec.transformed_from_u(value) for spec, value in zip(self.parameter_specs, u))

    def _validate_u(self, u: Sequence[float]) -> None:
        if len(u) != len(self.parameter_specs):
            raise SyntheticExecutionError("runner-parameter-count-mismatch")
        for index, value in enumerate(u):
            number = _finite(value, "runner.u")
            if number < 0.0 or number > 1.0:
                raise SyntheticExecutionError("runner-parameter-out-of-bound:" + str(index))

    def run(self, u: Sequence[float], split: str, stage: str) -> List[Dict[str, Any]]:
        self._validate_u(u)
        if split not in ("training", "holdout"):
            raise SyntheticExecutionError("runner-split-invalid")
        point = tuple(float(value) for value in u)
        physical_parameters = self.physical_parameters(point)
        self.calls.append({
            "split": split,
            "stage": stage,
            "point": _point_key(point),
            "case_requests": self.case_request_counts[split],
        })
        output: List[Dict[str, Any]] = []
        for case_value in _require_list(self.fixture.get("cases"), "fixture.cases"):
            case = _require_object(case_value, "fixture.case")
            if case.get("split") != split:
                continue
            self._case_request(case, physical_parameters)
            prediction = copy.deepcopy(_require_object(case.get("prediction"), "fixture.prediction"))
            if split == "training":
                for entry in self.response_entries:
                    if entry["case_id"] != case.get("case_id"):
                        continue
                    delta = point[self.parameter_index[entry["parameter"]]] - self.anchor_u[
                        self.parameter_index[entry["parameter"]]
                    ]
                    observable = entry["observable"]
                    history = prediction.get(observable)
                    if not isinstance(history, list) or entry["sample_index"] >= len(history):
                        raise SyntheticExecutionError("response-target-invalid")
                    history[entry["sample_index"]] = _finite(
                        history[entry["sample_index"]], "response.anchor"
                    ) + entry["slope"] * delta
            output.append({
                "case_id": case.get("case_id"),
                "sample_times_s": copy.deepcopy(case.get("sample_times_s")),
                "prediction": prediction,
            })
        return output


class ObjectiveAdapter:
    def __init__(self, contract: Dict[str, Any], fixture: Dict[str, Any], runner: SyntheticAffineRunner):
        self.contract = contract
        self.fixture = fixture
        self.runner = runner
        self.case_by_id = {
            case["case_id"]: case for case in _require_list(fixture.get("cases"), "fixture.cases")
        }
        partitions = _require_object(contract.get("partitions"), "contract.partitions")
        self.split_case_ids = {
            "training": tuple(partitions["training_case_ids"]),
            "holdout": tuple(partitions["holdout_case_ids"]),
        }
        metrics = _require_object(
            contract.get("metrics_weights_tolerances"), "contract.metrics"
        )
        self.scales = _require_object(metrics.get("normalization_scales"), "contract.scales")
        self.weights = _require_object(metrics.get("component_weights"), "contract.weights")

    def evaluate(self, u: Sequence[float], split: str, stage: str) -> ObjectiveEvaluation:
        try:
            outputs = self.runner.run(u, split, stage)
            output_by_id = {output["case_id"]: output for output in outputs}
            expected_ids = self.split_case_ids[split]
            if tuple(output["case_id"] for output in outputs) != expected_ids:
                raise SyntheticExecutionError("runner-case-order-mismatch:" + split)
            case_objectives: Dict[str, float] = {}
            residuals: List[float] = []
            case_count = len(expected_ids)
            if case_count == 0:
                raise SyntheticExecutionError("runner-empty-split:" + split)
            for case_id in expected_ids:
                case = self.case_by_id[case_id]
                output = output_by_id[case_id]
                case_residuals, case_objective = self._case_residuals(
                    case, output, case_count
                )
                residuals.extend(case_residuals)
                case_objectives[case_id] = case_objective
            objective = math.fsum(case_objectives[case_id] for case_id in expected_ids) / case_count
            norm_squared = math.fsum(value * value for value in residuals)
            tolerance = _finite(
                self.contract["metrics_weights_tolerances"]["synthetic_regression_tolerances"]
                ["recomputed_objective_absolute"],
                "objective.identity-tolerance",
            )
            if abs(norm_squared - objective) > tolerance:
                raise SyntheticExecutionError("scaled-residual-objective-identity-failed")
            if not math.isfinite(objective) or objective < 0.0:
                raise SyntheticExecutionError("objective-invalid")
            return ObjectiveEvaluation(True, objective, tuple(residuals), case_objectives)
        except SyntheticExecutionError as exc:
            return ObjectiveEvaluation.invalid(str(exc))

    def _case_residuals(
        self,
        case: Dict[str, Any],
        output: Dict[str, Any],
        split_case_count: int,
    ) -> Tuple[List[float], float]:
        reference = _require_object(case.get("reference"), "case.reference")
        prediction = _require_object(output.get("prediction"), "runner.prediction")
        times = _require_list(case.get("sample_times_s"), "case.sample-times")
        if output.get("sample_times_s") != times or not times:
            raise SyntheticExecutionError("runner-sample-grid-mismatch")
        time_count = len(times)

        tau_reference = _finite(reference.get("ignition_delay_s"), "reference.ignition-delay")
        tau_prediction = _finite(prediction.get("ignition_delay_s"), "prediction.ignition-delay")
        if tau_reference <= 0.0 or tau_prediction <= 0.0:
            raise SyntheticExecutionError("runner-ignition-delay-not-positive")
        sigma_log = _finite(self.scales.get("ignition_delay_sigma_log"), "scale.idt")
        idt_unscaled = math.log(tau_prediction / tau_reference) / sigma_log
        idt_component = idt_unscaled * idt_unscaled
        residuals = [
            math.sqrt(_finite(self.weights["ignition_delay"], "weight.idt") / split_case_count)
            * idt_unscaled
        ]
        components: Dict[str, float] = {"ignition_delay": idt_component}

        source_keys = {
            "qdot": ("qdot_w_m3", "qdot_w_m3"),
            "temperature": ("temperature_k", "temperature_k"),
            "C3H8": ("mass_fractions", "C3H8_mass_fraction"),
            "O2": ("mass_fractions", "O2_mass_fraction"),
            "CO": ("mass_fractions", "CO_mass_fraction"),
            "CO2": ("mass_fractions", "CO2_mass_fraction"),
            "H2O": ("mass_fractions", "H2O_mass_fraction"),
        }
        for observable in HISTORY_OBSERVABLES:
            source_key, scale_key = source_keys[observable]
            if source_key == "mass_fractions":
                reference_values = _require_object(
                    reference.get("mass_fractions"), "reference.mass-fractions"
                ).get(observable)
                prediction_values = _require_object(
                    prediction.get("mass_fractions"), "prediction.mass-fractions"
                ).get(observable)
            else:
                reference_values = reference.get(source_key)
                prediction_values = prediction.get(source_key)
            reference_history = _require_list(reference_values, "reference." + observable)
            prediction_history = _require_list(prediction_values, "prediction." + observable)
            if len(reference_history) != time_count or len(prediction_history) != time_count:
                raise SyntheticExecutionError("runner-history-length-mismatch:" + observable)
            scale = _finite(self.scales.get(scale_key), "scale." + scale_key)
            weight = _finite(self.weights.get(observable), "weight." + observable)
            if scale <= 0.0 or weight < 0.0:
                raise SyntheticExecutionError("runner-scale-or-weight-invalid:" + observable)
            factor = math.sqrt(weight / (split_case_count * time_count)) / scale
            raw_squared: List[float] = []
            for index, (predicted_value, reference_value) in enumerate(
                zip(prediction_history, reference_history)
            ):
                predicted = _finite(predicted_value, "prediction." + observable)
                expected = _finite(reference_value, "reference." + observable)
                if observable == "temperature" and predicted <= 0.0:
                    raise SyntheticExecutionError("runner-temperature-not-positive")
                if observable in ("C3H8", "O2", "CO", "CO2", "H2O") and not (
                    0.0 <= predicted <= 1.0
                ):
                    raise SyntheticExecutionError("runner-mass-fraction-out-of-range:" + observable)
                difference = predicted - expected
                raw_squared.append((difference / scale) ** 2)
                residuals.append(factor * difference)
            components[observable] = math.fsum(raw_squared) / time_count

        case_objective = math.fsum(
            _finite(self.weights[name], "weight." + name) * components[name]
            for name in OBSERVABLE_ORDER
        )
        return residuals, case_objective


def _gram_matrix(matrix: Sequence[Sequence[float]], columns: int) -> List[List[float]]:
    if not matrix:
        raise SyntheticExecutionError("sensitivity-matrix-empty")
    gram = [[0.0 for _ in range(columns)] for _ in range(columns)]
    for row in matrix:
        if len(row) != columns:
            raise SyntheticExecutionError("sensitivity-row-width-mismatch")
    for i in range(columns):
        for j in range(i, columns):
            value = math.fsum(row[i] * row[j] for row in matrix)
            if not math.isfinite(value):
                raise SyntheticExecutionError("gram-nonfinite")
            gram[i][j] = value
            gram[j][i] = value
    return gram


def symmetric_jacobi(
    matrix: Sequence[Sequence[float]],
    relative_tolerance: float,
    maximum_sweeps: int,
) -> Tuple[List[float], List[List[float]], int]:
    """Deterministic cyclic Jacobi eigensolver for a small symmetric matrix."""

    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise SyntheticExecutionError("jacobi-matrix-not-square")
    a = [[_finite(value, "jacobi.matrix") for value in row] for row in matrix]
    for i in range(n):
        for j in range(i + 1, n):
            if a[i][j] != a[j][i]:
                raise SyntheticExecutionError("jacobi-matrix-not-symmetric")
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    sweeps = 0
    for sweep in range(maximum_sweeps + 1):
        diagonal_scale = max(abs(a[i][i]) for i in range(n))
        off_diagonal = max(
            [abs(a[i][j]) for i in range(n) for j in range(i + 1, n)] or [0.0]
        )
        threshold = relative_tolerance * diagonal_scale
        if off_diagonal <= threshold or off_diagonal == 0.0:
            sweeps = sweep
            break
        if sweep == maximum_sweeps:
            raise SyntheticExecutionError("jacobi-did-not-converge")
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) <= threshold:
                    continue
                app = a[p][p]
                aqq = a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                if tau >= 0.0:
                    t = 1.0 / (tau + math.sqrt(1.0 + tau * tau))
                else:
                    t = -1.0 / (-tau + math.sqrt(1.0 + tau * tau))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                a[p][p] = app - t * apq
                a[q][q] = aqq + t * apq
                a[p][q] = 0.0
                a[q][p] = 0.0
                for k in range(n):
                    if k == p or k == q:
                        continue
                    akp = a[k][p]
                    akq = a[k][q]
                    new_kp = c * akp - s * akq
                    new_kq = s * akp + c * akq
                    a[k][p] = new_kp
                    a[p][k] = new_kp
                    a[k][q] = new_kq
                    a[q][k] = new_kq
                for k in range(n):
                    vkp = v[k][p]
                    vkq = v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    else:
        raise SyntheticExecutionError("jacobi-did-not-converge")

    unsorted = [(a[index][index], index) for index in range(n)]
    order = sorted(range(n), key=lambda index: (-unsorted[index][0], unsorted[index][1]))
    eigenvalues = [unsorted[index][0] for index in order]
    eigenvectors = [[v[row][index] for index in order] for row in range(n)]
    scale = max([abs(value) for value in eigenvalues] or [0.0])
    for index, value in enumerate(eigenvalues):
        if value < 0.0 and abs(value) <= max(scale * relative_tolerance, 1e-300):
            eigenvalues[index] = 0.0
        elif value < 0.0:
            raise SyntheticExecutionError("jacobi-negative-eigenvalue")
    return eigenvalues, eigenvectors, sweeps


def analyze_identifiability(
    selected_u: Sequence[float],
    parameter_count: int,
    sensitivity_step: float,
    relative_rank_tolerance: float,
    maximum_condition_number: float,
    maximum_correlation: float,
    jacobi_tolerance: float,
    jacobi_sweeps: int,
    evaluate_residuals: Callable[[Tuple[float, ...]], ObjectiveEvaluation],
) -> Dict[str, Any]:
    point = tuple(_finite(value, "identifiability.point") for value in selected_u)
    if len(point) != parameter_count:
        raise SyntheticExecutionError("identifiability-parameter-count-mismatch")
    cache: Dict[Tuple[str, ...], ObjectiveEvaluation] = {}
    call_count = 0

    def evaluate(candidate: Tuple[float, ...]) -> ObjectiveEvaluation:
        nonlocal call_count
        key = _point_key(candidate)
        if key not in cache:
            value = evaluate_residuals(candidate)
            call_count += 1
            cache[key] = value
        return cache[key]

    base = evaluate(point)
    if not base.valid:
        return {"status": "failed-base", "reason": base.reason, "runner_calls": call_count}
    row_count = len(base.residuals)
    if row_count <= parameter_count:
        return {
            "status": "insufficient-residual-rows",
            "row_count": row_count,
            "parameter_count": parameter_count,
            "runner_calls": call_count,
        }

    columns: List[List[float]] = []
    difference_modes: List[str] = []
    for parameter in range(parameter_count):
        lower_available = point[parameter] - sensitivity_step >= 0.0
        upper_available = point[parameter] + sensitivity_step <= 1.0
        if lower_available and upper_available:
            lower = list(point)
            upper = list(point)
            lower[parameter] -= sensitivity_step
            upper[parameter] += sensitivity_step
            lower_value = evaluate(tuple(lower))
            upper_value = evaluate(tuple(upper))
            if not lower_value.valid or not upper_value.valid:
                return {
                    "status": "failed-perturbation",
                    "parameter_index": parameter,
                    "runner_calls": call_count,
                }
            column = [
                (upper_residual - lower_residual) / (2.0 * sensitivity_step)
                for upper_residual, lower_residual in zip(
                    upper_value.residuals, lower_value.residuals
                )
            ]
            difference_modes.append("central")
        elif upper_available:
            upper = list(point)
            upper[parameter] += sensitivity_step
            upper_value = evaluate(tuple(upper))
            if not upper_value.valid:
                return {
                    "status": "failed-perturbation",
                    "parameter_index": parameter,
                    "runner_calls": call_count,
                }
            column = [
                (upper_residual - base_residual) / sensitivity_step
                for upper_residual, base_residual in zip(upper_value.residuals, base.residuals)
            ]
            difference_modes.append("forward-one-sided")
        elif lower_available:
            lower = list(point)
            lower[parameter] -= sensitivity_step
            lower_value = evaluate(tuple(lower))
            if not lower_value.valid:
                return {
                    "status": "failed-perturbation",
                    "parameter_index": parameter,
                    "runner_calls": call_count,
                }
            column = [
                (base_residual - lower_residual) / sensitivity_step
                for base_residual, lower_residual in zip(base.residuals, lower_value.residuals)
            ]
            difference_modes.append("backward-one-sided")
        else:
            return {
                "status": "no-inward-perturbation",
                "parameter_index": parameter,
                "runner_calls": call_count,
            }
        if len(column) != row_count or any(not math.isfinite(value) for value in column):
            return {
                "status": "invalid-sensitivity-column",
                "parameter_index": parameter,
                "runner_calls": call_count,
            }
        columns.append(column)

    sensitivity = [
        [columns[column][row] for column in range(parameter_count)] for row in range(row_count)
    ]
    gram = _gram_matrix(sensitivity, parameter_count)
    eigenvalues, eigenvectors, jacobi_used = symmetric_jacobi(
        gram, jacobi_tolerance, jacobi_sweeps
    )
    singular_values = [math.sqrt(value) for value in eigenvalues]
    largest = singular_values[0]
    rank = sum(
        1 for value in singular_values if value > largest * relative_rank_tolerance
    )
    condition_number = (
        largest / singular_values[parameter_count - 1]
        if rank == parameter_count and singular_values[parameter_count - 1] > 0.0
        else math.inf
    )
    residual_sum_squares = math.fsum(value * value for value in base.residuals)
    residual_variance = residual_sum_squares / (row_count - parameter_count)

    result: Dict[str, Any] = {
        "status": "nonidentifiable",
        "coordinate_basis": "normalized-transformed-parameter-u",
        "row_count": row_count,
        "parameter_count": parameter_count,
        "difference_modes": difference_modes,
        "runner_calls": call_count,
        "sensitivity_matrix": sensitivity,
        "gram_matrix": gram,
        "jacobi_sweeps": jacobi_used,
        "singular_values": singular_values,
        "rank": rank,
        "required_rank": parameter_count,
        "condition_number": condition_number,
        "maximum_condition_number": maximum_condition_number,
        "residual_sum_squares": residual_sum_squares,
        "residual_variance": residual_variance,
        "covariance_matrix": None,
        "correlation_matrix": None,
        "maximum_absolute_offdiagonal_correlation": None,
        "rank_condition_correlation_pass": False,
        "uncertainty_outputs": {
            "status": "not-publishable-synthetic-covariance-diagnostic-only",
            "marginal_and_joint_parameter_intervals": "not-computed-no-frozen-coverage-distribution-rule",
            "training_and_holdout_predictive_intervals": "not-computed-no-frozen-observation-uncertainty-model",
            "domain_of_applicability": "not-established-by-synthetic-fixture",
            "physical_evidence": False,
        },
    }
    if rank != parameter_count or not math.isfinite(condition_number):
        return result

    covariance = [[0.0 for _ in range(parameter_count)] for _ in range(parameter_count)]
    for i in range(parameter_count):
        for j in range(parameter_count):
            covariance[i][j] = residual_variance * math.fsum(
                eigenvectors[i][k] * eigenvectors[j][k] / (singular_values[k] ** 2)
                for k in range(parameter_count)
            )
            if not math.isfinite(covariance[i][j]):
                raise SyntheticExecutionError("covariance-nonfinite")
    diagonal = [covariance[index][index] for index in range(parameter_count)]
    if any(value <= 0.0 for value in diagonal):
        return result
    correlation = [[0.0 for _ in range(parameter_count)] for _ in range(parameter_count)]
    maximum_absolute = 0.0
    for i in range(parameter_count):
        for j in range(parameter_count):
            value = covariance[i][j] / math.sqrt(diagonal[i] * diagonal[j])
            if not math.isfinite(value):
                raise SyntheticExecutionError("correlation-nonfinite")
            correlation[i][j] = value
            if i != j:
                maximum_absolute = max(maximum_absolute, abs(value))
    passes = condition_number <= maximum_condition_number and maximum_absolute <= maximum_correlation
    result.update({
        "status": "synthetic-full-rank" if passes else "synthetic-identifiability-threshold-failed",
        "covariance_matrix": covariance,
        "correlation_matrix": correlation,
        "maximum_absolute_offdiagonal_correlation": maximum_absolute,
        "rank_condition_correlation_pass": passes,
    })
    return result


def _verify_execution_contract(document: Dict[str, Any]) -> None:
    _exact_keys(
        document,
        {
            "schema", "contract_version", "status", "scope", "approvals",
            "frozen_v3_authorities", "synthetic_runner", "execution_protocol",
            "claims", "stop_conditions",
        },
        "execution.root",
    )
    _exact(document["schema"], "rfam.r2-joint-fit-synthetic-execution-contract/v1", "execution.schema")
    _exact(document["contract_version"], 1, "execution.contract_version")
    _exact(document["status"], "blocked-synthetic-optimizer-execution-only", "execution.status")
    _exact(
        document["scope"],
        "fixture-backed-affine-runner-and-deterministic-protocol-execution-no-physical-fit",
        "execution.scope",
    )
    approvals = _require_object(document["approvals"], "execution.approvals")
    _exact(approvals, EXECUTION_APPROVALS, "execution.approvals")
    authorities = _require_object(
        document["frozen_v3_authorities"], "execution.frozen_v3_authorities"
    )
    _exact(authorities, FROZEN_V3_AUTHORITIES, "execution.frozen_v3_authorities")

    synthetic_runner = _require_object(document["synthetic_runner"], "execution.synthetic_runner")
    _exact_keys(
        synthetic_runner,
        {
            "adapter_id", "input_class", "anchor_parameter_vector", "anchor_prediction",
            "parameter_coordinate", "runner_call_scope", "training_response_entries",
            "holdout_parameter_response", "initial_state_authority", "invalid_output_policy",
            "physical_interpretation",
        },
        "execution.synthetic_runner",
    )
    runner_scalars = {
        "adapter_id": "fixture-anchored-affine-response-v1",
        "input_class": "SYNTHETIC_NONPHYSICAL_TEST_ONLY",
        "anchor_parameter_vector": "synthetic_fixture.synthetic_trial_parameter_vector",
        "anchor_prediction": "synthetic_fixture.cases[].prediction",
        "parameter_coordinate": "optimizer-normalized-transformed-coordinate-u",
        "runner_call_scope": "one-parameter-vector-and-one-explicit-split-per-call",
        "holdout_parameter_response": "none-return-frozen-fixture-prediction",
        "initial_state_authority": "frozen-fixture-ordered-initial-mass-fraction-vector",
        "invalid_output_policy": "nonfinite-missing-misaligned-or-out-of-range-output-invalidates-objective-point",
        "physical_interpretation": "none-affine-response-exists-only-to-exercise-software-protocol",
    }
    for key, expected in runner_scalars.items():
        _exact(synthetic_runner[key], expected, "execution.synthetic_runner." + key)
    response_entries = _require_list(
        synthetic_runner["training_response_entries"],
        "execution.synthetic_runner.training_response_entries",
    )
    response_signatures = []
    for index, entry_value in enumerate(response_entries):
        entry = _require_object(entry_value, "execution.response-entry." + str(index))
        _exact_keys(
            entry,
            {
                "parameter", "case_id", "observable", "sample_index",
                "raw_slope_per_normalized_coordinate",
            },
            "execution.response-entry." + str(index),
        )
        sample_index = entry["sample_index"]
        if isinstance(sample_index, bool) or not isinstance(sample_index, int):
            raise SyntheticExecutionError("contract-mismatch:execution.response-entry.sample_index")
        slope = entry["raw_slope_per_normalized_coordinate"]
        _finite(slope, "execution.response-entry.slope")
        response_signatures.append((
            entry["parameter"], entry["case_id"], entry["observable"], sample_index, slope,
        ))
    _exact(tuple(response_signatures), EXPECTED_RESPONSE_ENTRIES, "execution.response-entries")

    protocol = _require_object(document["execution_protocol"], "execution.execution_protocol")
    _exact(protocol, EXECUTION_PROTOCOL, "execution.execution_protocol")
    claims = _require_object(document["claims"], "execution.claims")
    _exact(claims, EXECUTION_CLAIMS, "execution.claims")
    _exact(document["stop_conditions"], EXECUTION_STOP_CONDITIONS, "execution.stop_conditions")


def _load_authorities(root: Path) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    _, execution_raw = _read_bound(root, EXECUTION_AUTHORITY, "execution-contract")
    execution = _load_json_bytes(execution_raw, "execution-contract")
    _verify_execution_contract(execution)
    loaded: Dict[str, Dict[str, Any]] = {}
    for name, authority in FROZEN_V3_AUTHORITIES.items():
        _, raw = _read_bound(root, authority, name)
        if name in ("objective_contract", "synthetic_fixture"):
            loaded[name] = _load_json_bytes(raw, name)
    return execution, loaded["objective_contract"], loaded["synthetic_fixture"]


def _parameter_specs(contract: Dict[str, Any]) -> Tuple[ParameterSpec, ...]:
    bounds = _require_object(contract.get("parameter_bounds"), "contract.parameter_bounds")
    entries = _require_list(bounds.get("ordered_fit_parameters"), "contract.fit_parameters")
    specs = []
    for entry_value in entries:
        entry = _require_object(entry_value, "contract.fit-parameter")
        name = entry.get("name")
        transform = entry.get("transform")
        lower = _finite(entry.get("lower"), "parameter.lower")
        upper = _finite(entry.get("upper"), "parameter.upper")
        if not isinstance(name, str) or transform not in ("linear", "log10") or not lower < upper:
            raise SyntheticExecutionError("fit-parameter-invalid")
        if transform == "log10":
            if lower <= 0.0:
                raise SyntheticExecutionError("log10-parameter-lower-not-positive")
            transformed_lower = math.log10(lower)
            transformed_upper = math.log10(upper)
        else:
            transformed_lower = lower
            transformed_upper = upper
        specs.append(ParameterSpec(
            name, transform, lower, upper, transformed_lower, transformed_upper
        ))
    return tuple(specs)


def _start_points(contract: Dict[str, Any], parameter_count: int) -> List[Tuple[str, Tuple[float, ...]]]:
    optimizer = _require_object(contract.get("optimizer_protocol"), "contract.optimizer")
    sequence = _require_object(optimizer.get("start_sequence"), "contract.start_sequence")
    bases = _require_list(sequence.get("bases_in_parameter_order"), "contract.halton-bases")
    indices = _require_list(sequence.get("halton_indices"), "contract.halton-indices")
    if len(bases) != parameter_count or indices != [1, 16]:
        raise SyntheticExecutionError("start-sequence-invalid")
    starts = [("midpoint", tuple(0.5 for _ in range(parameter_count)))]
    for index in range(indices[0], indices[1] + 1):
        starts.append((
            "halton-" + str(index),
            tuple(radical_inverse(index, int(base)) for base in bases),
        ))
    return starts


def _round_float(value: float) -> float:
    if not math.isfinite(value):
        raise SyntheticExecutionError("result-nonfinite")
    rounded = float(format(value, ".15g"))
    return 0.0 if rounded == 0.0 else rounded


def _round_result(value: Any) -> Any:
    if isinstance(value, float):
        return _round_float(value)
    if isinstance(value, list):
        return [_round_result(item) for item in value]
    if isinstance(value, tuple):
        return [_round_result(item) for item in value]
    if isinstance(value, dict):
        return {key: _round_result(item) for key, item in value.items()}
    return value


def canonical_json_bytes(document: Dict[str, Any]) -> bytes:
    rounded = _round_result(document)
    return (json.dumps(rounded, indent=2, ensure_ascii=True) + "\n").encode("utf-8")


def execute(root: Path) -> Dict[str, Any]:
    root = root.resolve()
    if not root.is_dir():
        raise SyntheticExecutionError("root-not-directory")
    execution, contract, fixture = _load_authorities(root)
    specs = _parameter_specs(contract)
    runner = SyntheticAffineRunner(contract, fixture, execution, specs)
    adapter = ObjectiveAdapter(contract, fixture, runner)
    starts = _start_points(contract, len(specs))
    optimizer = _require_object(contract["optimizer_protocol"], "contract.optimizer")
    search = _require_object(optimizer["search"], "contract.search")
    config = SearchConfig(
        initial_step=_finite(search["initial_normalized_step"], "search.initial-step"),
        contraction=_finite(search["step_contraction"], "search.contraction"),
        minimum_step=_finite(search["minimum_normalized_step"], "search.minimum-step"),
        maximum_evaluations=int(search["maximum_objective_evaluations_per_start"]),
        coordinate_count=len(specs),
    )
    if config.maximum_evaluations <= 0:
        raise SyntheticExecutionError("search-budget-invalid")

    search_results: List[SearchResult] = []
    for start_index, (start_kind, start) in enumerate(starts):
        result = run_pattern_search(
            start_index,
            start_kind,
            start,
            specs,
            config,
            lambda point, start_index=start_index: adapter.evaluate(
                point, "training", "selection-start-" + str(start_index)
            ),
        )
        search_results.append(result)
    valid_results = [result for result in search_results if result.valid]
    if not valid_results:
        raise SyntheticExecutionError("all-starts-invalid")
    selected = min(
        valid_results,
        key=lambda result: (
            result.evaluation.objective,
            _transformed_key(result.selected_u, specs),
        ),
    )
    selection_call_count = len(runner.calls)
    if any(call["split"] != "training" for call in runner.calls):
        raise SyntheticExecutionError("holdout-called-during-selection")
    if selection_call_count != sum(result.evaluation_count for result in search_results):
        raise SyntheticExecutionError("selection-runner-accounting-mismatch")

    holdout = adapter.evaluate(selected.selected_u, "holdout", "post-selection-holdout")
    if not holdout.valid:
        raise SyntheticExecutionError("holdout-evaluation-invalid:" + holdout.reason)
    if sum(1 for call in runner.calls if call["split"] == "holdout") != 1:
        raise SyntheticExecutionError("holdout-runner-call-count-mismatch")

    ident = _require_object(
        contract["identifiability_uncertainty"], "contract.identifiability"
    )
    local = _require_object(ident["local_sensitivity"], "contract.local-sensitivity")
    requirements = _require_object(
        ident["identifiability_requirements"], "contract.identifiability-requirements"
    )
    execution_protocol = _require_object(
        execution["execution_protocol"], "execution.execution-protocol"
    )
    calls_before_identifiability = len(runner.calls)
    identifiability = analyze_identifiability(
        selected.selected_u,
        len(specs),
        _finite(local["step"], "identifiability.step"),
        _finite(requirements["svd_relative_rank_tolerance"], "identifiability.rank-tolerance"),
        _finite(requirements["maximum_condition_number"], "identifiability.max-condition"),
        _finite(
            requirements["maximum_absolute_parameter_correlation"],
            "identifiability.max-correlation",
        ),
        _finite(
            execution_protocol["jacobi_relative_offdiagonal_tolerance"],
            "execution.jacobi-tolerance",
        ),
        int(execution_protocol["jacobi_maximum_sweeps"]),
        lambda point: adapter.evaluate(point, "training", "identifiability"),
    )
    ident_calls = len(runner.calls) - calls_before_identifiability
    if ident_calls != identifiability.get("runner_calls"):
        raise SyntheticExecutionError("identifiability-runner-accounting-mismatch")

    per_start = []
    for result in search_results:
        per_start.append({
            "start_index": result.start_index,
            "start_kind": result.start_kind,
            "start_normalized_coordinates": list(result.start_u),
            "selected_normalized_coordinates": list(result.selected_u),
            "objective": result.evaluation.objective if result.valid else None,
            "valid": result.valid,
            "runner_evaluations": result.evaluation_count,
            "cache_hits": result.cache_hits,
            "out_of_bound_skips": result.out_of_bound_skips,
            "completed_sweeps": result.completed_sweeps,
            "step_contractions": result.contractions,
            "accepted_coordinates": result.accepted_coordinates,
            "termination": result.termination,
            "trace_event_count": result.trace_event_count,
            "trace_sha256": result.trace_sha256,
        })

    result_document = {
        "schema": "rfam.r2-joint-fit-synthetic-execution-result/v1",
        "status": "blocked-synthetic-optimizer-execution-only",
        "watermark": "SYNTHETIC_NONPHYSICAL_TEST_ONLY",
        "authorities": {
            "execution_contract": dict(EXECUTION_AUTHORITY),
            "frozen_v3": copy.deepcopy(FROZEN_V3_AUTHORITIES),
        },
        "runner": {
            "adapter_id": execution["synthetic_runner"]["adapter_id"],
            "affine_formula": "anchor_prediction+slope*(u-anchor_u)",
            "parameter_order": [spec.name for spec in specs],
            "required_case_request_fields": list(runner.required_input_keys),
            "initial_state_authority": "ordered-initial-mass-fraction-vector-no-reconstruction",
            "case_requests_per_training_call": runner.case_request_counts["training"],
            "case_requests_per_holdout_call": runner.case_request_counts["holdout"],
            "anchor_normalized_coordinates": list(runner.anchor_u),
            "selection_training_calls": selection_call_count,
            "holdout_calls_during_selection": 0,
            "holdout_calls_after_selection": 1,
            "identifiability_training_calls": ident_calls,
            "total_calls": len(runner.calls),
        },
        "optimizer": {
            "algorithm": search["algorithm"],
            "start_count": len(starts),
            "selected_start_index": selected.start_index,
            "selected_start_kind": selected.start_kind,
            "selected_normalized_coordinates": list(selected.selected_u),
            "selected_transformed_parameters": list(
                runner.transformed_parameters(selected.selected_u)
            ),
            "selected_parameter_vector": list(runner.physical_parameters(selected.selected_u)),
            "training_objective": selected.evaluation.objective,
            "training_case_objectives": selected.evaluation.case_objectives,
            "global_selection": search["global_start_selection"],
            "per_start": per_start,
        },
        "holdout": {
            "evaluation_count": 1,
            "objective": holdout.objective,
            "case_objectives": holdout.case_objectives,
            "used_during_selection": False,
        },
        "identifiability": identifiability,
        "claims": {
            "synthetic_optimizer_executed": True,
            "physical_optimizer_executed": False,
            "physical_fit_executed": False,
            "parameter_vector_publishable": False,
            "synthetic_identifiability_is_physical_evidence": False,
            "physical_uncertainty_claimed": False,
            "runner_equivalent_to_rfam_or_cantera": False,
            "hem_chemistry_enabled": False,
            "fluent_execution_approved": False,
            "reacting_hot_start_approved": False,
            "production_ready": False,
        },
        "scope_limits": [
            "synthetic affine fixture only; no RFAM or Cantera kinetics executed",
            "no detailed mechanism or experimental training or holdout data selected",
            "covariance is a deterministic software diagnostic, not physical uncertainty",
            "parameter and predictive intervals remain undefined and unpublished",
            "no HEM Fluent hot-start or production capability change",
        ],
    }
    return result_document


def _write_atomic(path: Path, raw: bytes) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(str(temporary), str(path))
    except Exception:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = execute(args.root)
        raw = canonical_json_bytes(result)
        if args.output is not None:
            _write_atomic(args.output, raw)
    except SyntheticExecutionError as exc:
        print(FAIL_MARKER + " " + str(exc), file=sys.stderr)
        return 2
    except Exception as exc:  # Fail closed without leaking a traceback.
        print(FAIL_MARKER + " unexpected:" + type(exc).__name__, file=sys.stderr)
        return 2
    print(
        BLOCKED_MARKER
        + " synthetic_only=1 starts=" + str(result["optimizer"]["start_count"])
        + " selected_start=" + str(result["optimizer"]["selected_start_index"])
        + " training_objective=" + format(result["optimizer"]["training_objective"], ".17g")
        + " holdout_objective=" + format(result["holdout"]["objective"], ".17g")
        + " selection_calls=" + str(result["runner"]["selection_training_calls"])
        + " holdout_calls=1 identifiability_rank=" + str(result["identifiability"].get("rank", 0))
        + " physical_fit_executed=0 production_ready=0"
        + " result_sha256=" + _sha256(raw)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
