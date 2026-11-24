#!/usr/bin/env python3
"""Evaluate one frozen RFAM two-step candidate on the WSL C3Mech matrix."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence


class CandidateError(RuntimeError):
    """Fail-closed candidate evaluation error."""


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPOSITORY_ROOT.parent
BASE_TOOL_PATH = REPOSITORY_ROOT / "tools/run_r2_c3mech_wsl_screen.py"
CANDIDATE_RESULT_REL = Path("tmp/a408-r2-c3mech-fit-screen-result-v1.json")
REFERENCE_REL = Path("tmp/r2-c3mech-wsl-reference-v1-final.json")
RUNNER_REL = Path("tmp/rfam_r2_0d_fit_runner")
CANDIDATE_YAML_REL = Path("tmp/rfam-twostep-c3mech-internal-candidate-v1.yaml")
EXPECTED_FILES = {
    "candidate-result": (
        CANDIDATE_RESULT_REL,
        10996,
        "08426f2dadce6b54139d455122858746d40c72ad68166b0b89bf9f9b6f413c04",
    ),
    "reference": (
        REFERENCE_REL,
        85975,
        "36437c95051f959c6499ebbd235604207558834cc42e8e2e78fffc20222337be",
    ),
    "runner": (
        RUNNER_REL,
        97816,
        "d080d0ed5a17fb029fc5e01fbad5dbe6f728f88f0686f0292606a8e81344f9ad",
    ),
}
PARAMETER_NAMES = (
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
)
PARAMETERS = (
    7498942.09332456,
    133333.333333333,
    1.06875,
    0.803571428571428,
    115856.513372443,
    298557.692307692,
)
OBJECTIVE_ORDER = (
    "ignition_delay",
    "qdot",
    "temperature",
    "C3H8",
    "O2",
    "CO",
    "CO2",
    "H2O",
)
WEIGHTS = {
    "ignition_delay": 0.30,
    "qdot": 0.15,
    "temperature": 0.20,
    "C3H8": 0.10,
    "O2": 0.10,
    "CO": 0.05,
    "CO2": 0.05,
    "H2O": 0.05,
}
SCHEMA = "rfam.r2-c3mech-wsl-candidate-evaluation/v1"
STATUS = "INTERNAL_C3MECH_CANDIDATE_EVALUATED_NOT_PHYSICALLY_VALIDATED"
MARKER = "R2_C3MECH_WSL_CANDIDATE_BLOCKED"
FAIL_MARKER = "R2_C3MECH_WSL_CANDIDATE_FAIL"


def load_base_tool():
    spec = importlib.util.spec_from_file_location("rfam_c3mech_wsl_base", BASE_TOOL_PATH)
    if spec is None or spec.loader is None:
        raise CandidateError("base-tool-import-spec-failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as error:
        sys.modules.pop(spec.name, None)
        raise CandidateError("base-tool-import-failed") from error
    return module


base = load_base_tool()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_pairs(pairs: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise CandidateError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise CandidateError("nonfinite-json-constant:" + value)


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        document = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except CandidateError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CandidateError("invalid-json:" + label) from error
    if not isinstance(document, dict):
        raise CandidateError("json-root-not-object:" + label)
    return document


def require_file(root: Path, artifact_id: str) -> Path:
    relative, expected_size, expected_sha = EXPECTED_FILES[artifact_id]
    declared = root / relative
    if declared.is_symlink() or not declared.is_file():
        raise CandidateError("missing-or-aliased-artifact:" + artifact_id)
    path = declared.resolve()
    if path != declared.absolute():
        raise CandidateError("artifact-path-mismatch:" + artifact_id)
    if path.stat().st_size != expected_size or sha256(path) != expected_sha:
        raise CandidateError("artifact-identity-mismatch:" + artifact_id)
    return path


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise CandidateError("not-a-number:" + label)
    try:
        result = float(value)
    except ValueError as error:
        raise CandidateError("not-a-number:" + label) from error
    if not math.isfinite(result):
        raise CandidateError("nonfinite:" + label)
    return result


def stable(value: float, label: str) -> float:
    if not math.isfinite(value):
        raise CandidateError("nonfinite:" + label)
    return float(format(float(value), ".12g"))


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (
        json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("ascii")


def publish_no_replace(path: Path, payload: bytes) -> None:
    if os.path.lexists(str(path)):
        raise CandidateError("output-path-already-exists")
    if not path.parent.is_dir():
        raise CandidateError("output-parent-missing")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix="." + path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise CandidateError("output-publication-race") from error
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def validate_candidate(document: dict[str, Any]) -> None:
    if document.get("schema") != "rfam.r2-c3mech-fit-screen-result/v1":
        raise CandidateError("candidate-schema-mismatch")
    if document.get("status") != "blocked-detailed-mechanism-candidate-fit-screen-only":
        raise CandidateError("candidate-status-mismatch")
    selection = document.get("selection")
    if not isinstance(selection, dict) or selection.get("split") != "training-only":
        raise CandidateError("candidate-selection-split-mismatch")
    values = selection.get("selected_parameter_vector")
    if not isinstance(values, list) or tuple(values) != PARAMETERS:
        raise CandidateError("candidate-parameter-vector-mismatch")
    accounting = document.get("runner_accounting")
    if not isinstance(accounting, dict):
        raise CandidateError("candidate-accounting-missing")
    if (
        accounting.get("holdout_access_attempts_during_selection") != 0
        or accounting.get("holdout_objective_calls_after_selection") != 1
    ):
        raise CandidateError("candidate-holdout-accounting-mismatch")
    approvals = document.get("approvals")
    if not isinstance(approvals, dict):
        raise CandidateError("candidate-approvals-missing")
    allowed_true = {"diagnostic_optimizer_execution_allowed"}
    for key, value in approvals.items():
        expected = key in allowed_true
        if value is not expected:
            raise CandidateError("candidate-approval-mismatch:" + key)
    claims = document.get("claims")
    if not isinstance(claims, dict) or any(
        claims.get(key) is not False
        for key in (
            "experimental_fit_executed",
            "parameter_vector_publishable",
            "physical_validation_approved",
            "hem_chemistry_enabled",
            "fluent_execution_approved",
            "production_ready",
        )
    ):
        raise CandidateError("candidate-claim-boundary-weakened")


def strip_fit_header(stdout: str, case_id: str) -> str:
    lines = stdout.splitlines(keepends=True)
    expected = ["fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n"]
    expected.extend(
        "fit_" + name + "=" + format(value, ".17g") + "\n"
        for name, value in zip(PARAMETER_NAMES, PARAMETERS)
    )
    if lines[: len(expected)] != expected:
        raise CandidateError("runner-fit-header-mismatch:" + case_id)
    return "".join(lines[len(expected) :])


def run_case(runner: Path, case: dict[str, Any]) -> tuple[dict[str, str], list[dict[str, Any]]]:
    initial = case["initial_mass_fractions"]
    command = [
        str(runner),
        *[format(value, ".17g") for value in PARAMETERS],
        case["id"],
        format(float(case["temperature_k"]), ".17g"),
        format(float(case["pressure_pa"]), ".17g"),
        str(len(base.SAMPLE_TIMES)),
        "0",
        *[format(float(initial[name]), ".17g") for name in base.SPECIES],
        *[format(value, ".17g") for value in base.SAMPLE_TIMES],
    ]
    result = subprocess.run(
        command, check=False, text=True, capture_output=True, timeout=60
    )
    if result.returncode != 0:
        raise CandidateError(
            "runner-rejected:" + case["id"] + ":rc=" + str(result.returncode)
        )
    trajectory = strip_fit_header(result.stdout, case["id"])
    try:
        return base.parse_runner_output(
            trajectory, result.stderr, case["id"], len(base.SAMPLE_TIMES)
        )
    except Exception as error:
        raise CandidateError("runner-output-invalid:" + case["id"]) from error


def rms(values: Sequence[float]) -> float:
    if not values:
        raise CandidateError("empty-rms-input")
    return math.sqrt(math.fsum(value * value for value in values) / len(values))


def score_case(
    case: dict[str, Any],
    headers: dict[str, str],
    predictions: list[dict[str, Any]],
) -> dict[str, Any]:
    paired = []
    for detailed, predicted in zip(case["samples"], predictions):
        if detailed["status"] != "eligible-within-mechanism-thermo-intersection":
            continue
        temperature = finite(predicted["temperature_k"], "prediction-temperature")
        if 300.0 <= temperature <= 3400.0:
            paired.append((detailed, predicted))
    if len(paired) < 2:
        raise CandidateError("insufficient-joint-prefix:" + case["id"])

    components: dict[str, float] = {}
    detailed_idt = case["ignition_delay"]
    predicted_idt = None
    idt_ratio = None
    if detailed_idt["status"] == "scored-within-mechanism-thermo-intersection":
        if headers.get("ignition_detected") != "1":
            raise CandidateError("prediction-ignition-not-detected:" + case["id"])
        predicted_idt = finite(headers["ignition_delay_s"], "prediction-idt")
        reference_idt = finite(detailed_idt["time_s"], "reference-idt")
        if predicted_idt <= 0.0 or reference_idt <= 0.0:
            raise CandidateError("ignition-delay-invalid:" + case["id"])
        idt_ratio = predicted_idt / reference_idt
        residual = math.log(idt_ratio) / 0.25
        components["ignition_delay"] = residual * residual

    qdot_reference = [
        finite(detailed["heat_release_rate_w_m3"], "reference-qdot")
        for detailed, _ in paired
    ]
    qdot_prediction = [
        finite(predicted["qdot_w_m3"], "prediction-qdot")
        for _, predicted in paired
    ]
    qdot_scale = max(max(abs(value) for value in qdot_reference), 1.0e8)
    qdot_errors = [
        (left - right) / qdot_scale
        for left, right in zip(qdot_prediction, qdot_reference)
    ]
    components["qdot"] = math.fsum(value * value for value in qdot_errors) / len(paired)

    temperature_errors = [
        finite(predicted["temperature_k"], "prediction-temperature")
        - finite(detailed["temperature_k"], "reference-temperature")
        for detailed, predicted in paired
    ]
    components["temperature"] = math.fsum(
        (value / 50.0) ** 2 for value in temperature_errors
    ) / len(paired)

    species_metrics: dict[str, Any] = {}
    for index, name in enumerate(base.SPECIES):
        errors = [
            finite(predicted["mass_fractions"][index], "prediction-Y")
            - finite(detailed["mass_fractions"][name], "reference-Y")
            for detailed, predicted in paired
        ]
        if name != "N2":
            components[name] = math.fsum((value / 0.02) ** 2 for value in errors) / len(paired)
        species_metrics[name] = {
            "rmse_mass_fraction": stable(rms(errors), "species-rmse"),
            "max_abs_mass_fraction": stable(
                max(abs(value) for value in errors), "species-max"
            ),
            "objective_weight": WEIGHTS.get(name, 0.0),
        }

    active_weight = math.fsum(WEIGHTS[name] for name in components)
    objective = math.fsum(
        WEIGHTS[name] * components[name] / active_weight for name in OBJECTIVE_ORDER
        if name in components
    )
    detailed_qdot0 = qdot_reference[0]
    predicted_qdot0 = qdot_prediction[0]
    return {
        "case_id": case["id"],
        "partition": case["partition"],
        "joint_prefix_sample_count": len(paired),
        "joint_prefix_last_time_s": stable(
            finite(paired[-1][0]["time_s"], "last-time"), "last-time"
        ),
        "objective": stable(objective, "case-objective"),
        "components": {
            name: stable(value, "component") for name, value in components.items()
        },
        "reference_ignition_delay": detailed_idt,
        "candidate_ignition_delay_s": (
            stable(predicted_idt, "candidate-idt") if predicted_idt is not None else None
        ),
        "candidate_to_reference_idt_ratio": (
            stable(idt_ratio, "idt-ratio") if idt_ratio is not None else None
        ),
        "initial_qdot_sign_matches": (
            (detailed_qdot0 == 0.0 and predicted_qdot0 == 0.0)
            or detailed_qdot0 * predicted_qdot0 > 0.0
        ),
        "initial_qdot_absolute_ratio": stable(
            abs(predicted_qdot0) / max(abs(detailed_qdot0), 1.0), "qdot0-ratio"
        ),
        "temperature_rmse_k": stable(rms(temperature_errors), "temperature-rmse"),
        "temperature_max_abs_error_k": stable(
            max(abs(value) for value in temperature_errors), "temperature-max"
        ),
        "qdot_normalized_rmse": stable(rms(qdot_errors), "qdot-rmse"),
        "species": species_metrics,
    }


def evaluate_split(
    runner: Path,
    cases: Sequence[dict[str, Any]],
    split: str,
) -> dict[str, Any]:
    expected = base.CASE_IDS["protocol-" + split]
    if tuple(case["id"] for case in cases) != expected:
        raise CandidateError("split-case-order-mismatch:" + split)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(run_case, runner, case) for case in cases]
        trajectories = [future.result() for future in futures]
    results = [
        score_case(case, headers, samples)
        for case, (headers, samples) in zip(cases, trajectories)
    ]
    return {
        "objective": stable(
            math.fsum(result["objective"] for result in results) / len(results),
            split + "-objective",
        ),
        "case_count": len(results),
        "case_process_attempts": len(results),
        "case_process_successes": len(results),
        "joint_prefix_sample_count": sum(
            result["joint_prefix_sample_count"] for result in results
        ),
        "scored_idt_count": sum(
            result["candidate_to_reference_idt_ratio"] is not None for result in results
        ),
        "initial_qdot_sign_mismatch_count": sum(
            not result["initial_qdot_sign_matches"] for result in results
        ),
        "maximum_temperature_error_k": stable(
            max(result["temperature_max_abs_error_k"] for result in results),
            split + "-temperature-max",
        ),
        "maximum_qdot_normalized_rmse": stable(
            max(result["qdot_normalized_rmse"] for result in results),
            split + "-qdot-max",
        ),
        "cases": results,
    }


def candidate_yaml_bytes() -> bytes:
    values = dict(zip(PARAMETER_NAMES, PARAMETERS))
    text = f"""description: |-
  Internal RFAM two-step candidate selected on C3Mech training cases only.
  This separate artifact does not replace RFAM defaults and is not an
  experimentally validated or production mechanism.

generator: rfam-r2-c3mech-wsl-candidate-v1
cantera-version: 3.1.0

units: {{length: m, time: s, quantity: mol, activation-energy: J/mol}}

elements:
- symbol: C
  atomic-weight: 12.0107
- symbol: H
  atomic-weight: 1.00794
- symbol: O
  atomic-weight: 15.9994
- symbol: N
  atomic-weight: 14.0067

phases:
- name: rfam-two-step-c3mech-internal-candidate
  thermo: ideal-gas
  kinetics: gas
  transport: none
  elements: [O, H, C, N]
  species:
  - ../realFluidAllMach/tests/reference/cantera_gri30_v3_0.yaml/species: [C3H8, O2, CO, CO2, H2O, N2]
  reactions: all
  state: {{T: 1200.0 K, P: 4.0 MPa}}

reactions:
- equation: C3H8 + 3.5 O2 => 3 CO + 4 H2O
  rate-constant: {{A: {format(values['R1_A_SI'], '.17g')}, b: 0.0, Ea: {format(values['R1_Ea_J_per_mol'], '.17g')}}}
  orders: {{C3H8: {format(values['R1_order_C3H8'], '.17g')}, O2: {format(values['R1_order_O2'], '.17g')}}}

- equation: CO + 0.5 O2 <=> CO2
  rate-constant: {{A: {format(values['R2_A_SI'], '.17g')}, b: 0.0, Ea: {format(values['R2_Ea_J_per_mol'], '.17g')}}}
"""
    return text.encode("ascii")


def evaluate(root: Path) -> tuple[dict[str, Any], bytes]:
    base.verify_authorities(root)
    candidate_path = require_file(root, "candidate-result")
    reference_path = require_file(root, "reference")
    runner = require_file(root, "runner")
    candidate = load_json(candidate_path, "candidate-result")
    validate_candidate(candidate)
    reference = load_json(reference_path, "reference")
    base.validate_reference(reference)

    training_cases = [
        case for case in reference["cases"] if case["partition"] == "protocol-training"
    ]
    holdout_cases = [
        case for case in reference["cases"] if case["partition"] == "protocol-holdout"
    ]
    access_order = ["candidate-selection-frozen", "training-evaluation"]
    training = evaluate_split(runner, training_cases, "training")
    access_order.extend(["holdout-execution-opened-after-training", "holdout-evaluation-once"])
    holdout = evaluate_split(runner, holdout_cases, "holdout")
    yaml_payload = candidate_yaml_bytes()
    temperatures = [float(case["temperature_k"]) for case in reference["cases"]]
    pressures = [float(case["pressure_pa"]) for case in reference["cases"]]
    phis = [float(case["equivalence_ratio"]) for case in reference["cases"]]
    return {
        "schema": SCHEMA,
        "status": STATUS,
        "source": {
            "training_selected_candidate_path": CANDIDATE_RESULT_REL.as_posix(),
            "training_selected_candidate_sha256": EXPECTED_FILES["candidate-result"][2],
            "expanded_reference_path": REFERENCE_REL.as_posix(),
            "expanded_reference_sha256": EXPECTED_FILES["reference"][2],
            "parameterized_runner_path": RUNNER_REL.as_posix(),
            "parameterized_runner_sha256": EXPECTED_FILES["runner"][2],
            "base_authority_composite_sha256": base.INPUT_COMPOSITE_SHA256,
            "evaluator_sha256": sha256(Path(__file__).resolve()),
        },
        "candidate": {
            "parameter_names": list(PARAMETER_NAMES),
            "parameter_vector": list(PARAMETERS),
            "selection_origin": "prior-C3Mech-training-only-pattern-search",
            "expanded_grid_refit_performed": False,
            "holdout_influenced_parameters": False,
            "candidate_yaml_logical_path": CANDIDATE_YAML_REL.as_posix(),
            "candidate_yaml_sha256": hashlib.sha256(yaml_payload).hexdigest(),
            "default_twostep_yaml_replaced": False,
            "runtime_installation_approved": False,
        },
        "objective": {
            "observable_order": list(OBJECTIVE_ORDER),
            "weights": WEIGHTS,
            "N2_policy": "reported-as-sixth-species-diagnostic-with-zero-fit-weight",
            "ignition_residual": "square(log(candidate/reference)/0.25)",
            "temperature_scale_k": 50.0,
            "species_mass_fraction_scale": 0.02,
            "qdot_scale": "max-absolute-reference-or-1e8-W-m3-floor",
            "case_aggregation": "arithmetic-mean-in-frozen-case-order",
        },
        "access_order": access_order,
        "runner_accounting": {
            "training_parameter_adjustment_calls": 0,
            "training_objective_calls": 1,
            "holdout_access_attempts_before_training_completed": 0,
            "holdout_objective_calls_after_selection": 1,
            "training_case_process_attempts": training["case_process_attempts"],
            "training_case_process_successes": training["case_process_successes"],
            "holdout_case_process_attempts": holdout["case_process_attempts"],
            "holdout_case_process_successes": holdout["case_process_successes"],
        },
        "training": training,
        "holdout": holdout,
        "applicability_envelope": {
            "reactor": "ideal-gas-adiabatic-constant-pressure-0D-screen",
            "initial_temperature_k": [min(temperatures), max(temperatures)],
            "pressure_pa": [min(pressures), max(pressures)],
            "equivalence_ratio": [min(phis), max(phis)],
            "sample_time_s": [base.SAMPLE_TIMES[0], base.SAMPLE_TIMES[-1]],
            "detailed_mechanism_thermo_intersection_k": [400.0, 2000.0],
            "scoring_policy": "joint-prefix-only-detailed-below-2000K-and-rfam-300-to-3400K",
            "phase_scope": "high-temperature-gas-only-no-VLE-no-HEM-no-solid",
            "extrapolation_allowed": False,
        },
        "rights_and_validation_boundary": {
            "c3mech_mechanism_code_license_candidate": "CC-BY-4.0",
            "experimental_training_rights_cleared": False,
            "experimental_holdout_rights_cleared": False,
            "experimental_values_ingested": False,
            "physical_validation_approved": False,
        },
        "approvals": {
            "candidate_parameter_publication_approved": False,
            "default_mechanism_replacement_approved": False,
            "phase_aware_chemistry_enabled": False,
            "fluent_execution_approved": False,
            "reacting_hot_start_approved": False,
            "production_ready": False,
        },
    }, yaml_payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=WORKSPACE_ROOT)
    parser.add_argument("--yaml-output", type=Path, required=True)
    parser.add_argument("--json-output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        root = args.workspace_root.resolve()
        yaml_output = args.yaml_output.resolve()
        json_output = args.json_output.resolve()
        if os.path.lexists(str(yaml_output)) or os.path.lexists(str(json_output)):
            raise CandidateError("output-path-already-exists")
        result, yaml_payload = evaluate(root)
        publish_no_replace(yaml_output, yaml_payload)
        publish_no_replace(json_output, canonical_bytes(result))
        print(
            MARKER
            + " training_objective="
            + format(result["training"]["objective"], ".12g")
            + " holdout_objective="
            + format(result["holdout"]["objective"], ".12g")
            + " holdout_calls=1 production_ready=0"
        )
        return 0
    except Exception as error:
        print(
            FAIL_MARKER + " " + type(error).__name__ + ":" + str(error),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
