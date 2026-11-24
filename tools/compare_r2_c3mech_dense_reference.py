#!/usr/bin/env python3
"""Evaluate the frozen coarse-grid R2 selection on dense C3Mech references.

No optimization occurs here.  The training reference is evaluated before the
holdout reference is opened, and the holdout result cannot alter the frozen
parameter vector.  The output is a diagnostic route-comparison artifact only.
"""

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
from typing import Any, Sequence


class DenseComparisonError(RuntimeError):
    """A fail-closed identity, schema, runner, or publication error."""


ARTIFACTS = {
    "training_reference": (
        Path("tmp/a408-r2-c3mech-44i0-dense-training-v1.json"),
        "ae629e1501531d803eec8a74db37dc1e840c5d6679528918e9f8f2edd3862247",
        55425,
    ),
    "holdout_reference": (
        Path("tmp/a408-r2-c3mech-44i0-dense-holdout-v1.json"),
        "e1cb70bee2ea084e7d1ed79ebe4481bc87c7c98b1df4a03e0b4fdd3cfa24b65e",
        54047,
    ),
    "selection": (
        Path("realFluidAllMach/tests/reference/r2_c3mech_coarse_selection_v1.json"),
        "2cd3842e2dac60cd82900abad760ac46a720ec051c3a5d32f50fccb6b5d6a562",
        883,
    ),
    "fit_module": (
        Path("realFluidAllMach/tools/run_r2_c3mech_fit_screen.py"),
        "140983b47636f3bfaec9eac8e8dcc008f5c5253b277c2dd769ee28b462e66c61",
        54647,
    ),
    "protocol_module": (
        Path("realFluidAllMach/tools/run_r2_0d_rfam_protocol.py"),
        "5a4642c37ce4b07ca97db82a6af218d39e039b8a73e417d20efb6094d43f58a8",
        20748,
    ),
    "runner": (
        Path("tmp/rfam_r2_0d_fit_runner"),
        "d080d0ed5a17fb029fc5e01fbad5dbe6f728f88f0686f0292606a8e81344f9ad",
        97816,
    ),
}
TRAINING_IDS = ("train_phi0p8_diluted", "train_phi1p0", "train_phi1p2")
HOLDOUT_IDS = ("holdout_phi0p6_diluted", "holdout_phi1p1", "holdout_phi1p4")
PARAMETER_NAMES = (
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
)
PARAMETER_VECTOR = (
    7498942.09332456,
    133333.333333333,
    1.06875,
    0.803571428571428,
    115856.513372443,
    298557.692307692,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DenseComparisonError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise DenseComparisonError("nonfinite-json-constant:" + value)


def require_file(root: Path, artifact_id: str) -> Path:
    relative, expected_sha, expected_size = ARTIFACTS[artifact_id]
    declared = root / relative
    if declared.is_symlink() or not declared.is_file():
        raise DenseComparisonError("missing-or-aliased-artifact:" + artifact_id)
    path = declared.resolve()
    if path.stat().st_size != expected_size or sha256(path) != expected_sha:
        raise DenseComparisonError("artifact-identity-mismatch:" + artifact_id)
    return path


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except DenseComparisonError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise DenseComparisonError("invalid-json:" + label) from error
    if not isinstance(value, dict):
        raise DenseComparisonError("json-root-not-object:" + label)
    return value


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise DenseComparisonError("module-import-spec-failed:" + name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as error:
        sys.modules.pop(name, None)
        raise DenseComparisonError("module-import-failed:" + name) from error
    return module


def load_selection(root: Path) -> dict[str, Any]:
    document = load_json(require_file(root, "selection"), "selection")
    expected_keys = {
        "schema", "status", "source_result", "selection_split",
        "holdout_influenced_selection", "parameter_names", "parameter_vector",
        "coarse_training_objective", "approvals",
    }
    if set(document) != expected_keys:
        raise DenseComparisonError("selection-schema-open")
    if (
        document["schema"] != "rfam.r2-c3mech-diagnostic-selection/v1"
        or document["status"] != "coarse-grid-training-only-selection-not-publishable"
        or document["selection_split"] != "training-only"
        or document["holdout_influenced_selection"] is not False
        or tuple(document["parameter_names"]) != PARAMETER_NAMES
        or tuple(document["parameter_vector"]) != PARAMETER_VECTOR
        or document["coarse_training_objective"] != 1.63488450935637
        or document["source_result"]
        != {
            "path": "tmp/a408-r2-c3mech-fit-screen-review-replacement-v2.json",
            "sha256": "08426f2dadce6b54139d455122858746d40c72ad68166b0b89bf9f9b6f413c04",
            "size": 10996,
        }
        or document["approvals"]
        != {
            "parameter_vector_publishable": False,
            "physical_validation_approved": False,
            "production_ready": False,
        }
    ):
        raise DenseComparisonError("selection-contract-mismatch")
    return document


def reference_view(document: dict[str, Any], expected_partition: str) -> dict[str, Any]:
    if document.get("schema") != "rfam.r2-c3mech-dense-calibration-reference/v1":
        raise DenseComparisonError("dense-reference-schema-mismatch:" + expected_partition)
    if document.get("partition") != expected_partition:
        raise DenseComparisonError("dense-reference-partition-mismatch:" + expected_partition)
    if document.get("sampling_policy") != {
        "candidate_or_reference_adaptive": False,
        "grid_kind": "fixed-predeclared-log-refinement",
        "grid_not_derived_from_rfam_or_fit": True,
        "grid_source": "checked-in-case-manifest",
        "holdout_grid_reused_after_selection": True,
        "sample_count": 42,
    }:
        raise DenseComparisonError("dense-sampling-policy-mismatch:" + expected_partition)
    source = document.get("source")
    if not isinstance(source, dict):
        raise DenseComparisonError("dense-source-invalid:" + expected_partition)
    expected_source_keys = {
        "mechanism_path", "mechanism_sha256", "mechanism_size", "license_path",
        "license_sha256", "license_size", "rights_report_path",
        "rights_report_sha256", "rights_report_size", "license_declaration",
        "base_generator_path", "base_generator_sha256", "case_manifest_path",
        "case_manifest_sha256", "case_manifest_size", "dense_generator_path",
        "dense_generator_sha256",
    }
    if set(source) != expected_source_keys:
        raise DenseComparisonError("dense-source-schema-open:" + expected_partition)
    if (
        source["mechanism_sha256"]
        != "d5cec179ffbb6221b1b91264e1439ec850658dd97a0298054c8dbc336f3ab148"
        or source["license_sha256"]
        != "7e7170e3cebf88a9f60c7b8421418323c09304da1af4d5e90f4da1dc1c8a2661"
        or source["rights_report_sha256"]
        != "6779ae447b8ae27d49f4b1011e91ab0747dc8c91435f44d5effcdf6cad5fd26f"
        or source["base_generator_sha256"]
        != "d8d4d3ad9cc4207039c49f9d2174d8e28ca8db601bea077a85d955d11d55bb5c"
        or source["dense_generator_sha256"]
        != "7b417d7d77f1a99690e560658dbf54cee890becf2bf5afc454535ab7e23da386"
        or source["license_declaration"]
        != "CC-BY-4.0-source-screen-candidate-not-fit-approval"
    ):
        raise DenseComparisonError("dense-source-authority-mismatch:" + expected_partition)
    view = {key: value for key, value in document.items() if key != "sampling_policy"}
    view["schema"] = "rfam.r2-c3mech-cantera-reference/v2"
    view["source"] = {
        "case_manifest_path": source["case_manifest_path"],
        "case_manifest_sha256": source["case_manifest_sha256"],
        "generator_sha256": source["base_generator_sha256"],
        "mechanism_path": source["mechanism_path"],
        "mechanism_sha256": source["mechanism_sha256"],
        "mechanism_size": source["mechanism_size"],
    }
    return view


def run_case(runner: Path, fit: Any, protocol: Any, parameters: Sequence[float], case: Any):
    command = [
        str(runner),
        *[format(float(value), ".17g") for value in parameters],
        case.case_id,
        format(case.temperature_k, ".17g"),
        format(case.pressure_pa, ".17g"),
        str(len(case.runner_times_s)),
        "0",
        *[format(value, ".17g") for value in case.initial_mass_fractions],
        *[format(value, ".17g") for value in case.runner_times_s],
    ]
    result = subprocess.run(command, check=False, text=True, capture_output=True, timeout=30)
    if result.returncode != 0:
        raise DenseComparisonError("runner-rejected:" + case.case_id + ":rc=" + str(result.returncode))
    trajectory = fit.strip_fit_parameter_header(result.stdout, parameters, case.case_id)
    headers, samples = protocol.parse_runner_output(
        trajectory, result.stderr, case.case_id, len(case.runner_times_s)
    )
    for index, (sample, expected_time) in enumerate(zip(samples, case.runner_times_s)):
        if abs(float(sample["time_s"]) - expected_time) > 2.0e-18:
            raise DenseComparisonError("runner-time-mismatch:" + case.case_id + ":" + str(index))
    return headers, samples


def evaluate_partition(runner: Path, fit: Any, protocol: Any, cases: Sequence[Any]) -> dict[str, Any]:
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(run_case, runner, fit, protocol, PARAMETER_VECTOR, case)
            for case in cases
        ]
        predictions = [future.result() for future in futures]
    residuals: list[float] = []
    case_objectives: dict[str, float] = {}
    case_components: dict[str, dict[str, float]] = {}
    weights = dict(fit.EXPECTED_OBJECTIVE["component_weights"])
    for case, (headers, samples) in zip(cases, predictions):
        objective, rows, components = fit.score_case(case, headers, samples, weights)
        case_objectives[case.case_id] = objective
        case_components[case.case_id] = components
        residuals.extend(value / math.sqrt(float(len(cases))) for value in rows)
    return {
        "objective": math.fsum(value * value for value in residuals),
        "case_objectives": case_objectives,
        "case_components": case_components,
        "eligible_sample_counts": {
            case.case_id: len(case.eligible_samples) for case in cases
        },
        "runner_sample_counts": {
            case.case_id: len(case.runner_times_s) for case in cases
        },
    }


def rounded(value: Any) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise DenseComparisonError("nonfinite-result")
        return float(format(value, ".15g"))
    if isinstance(value, list):
        return [rounded(item) for item in value]
    if isinstance(value, tuple):
        return [rounded(item) for item in value]
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    return value


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (json.dumps(rounded(document), indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")


def publish_no_replace(path: Path, raw: bytes) -> None:
    if os.path.lexists(path):
        raise DenseComparisonError("output-path-already-exists")
    if not path.parent.is_dir():
        raise DenseComparisonError("output-parent-missing")
    fd, temporary_name = tempfile.mkstemp(prefix=".r2-dense-compare-", suffix=".tmp", dir=str(path.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise DenseComparisonError("output-publication-race") from error
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def execute(root: Path) -> dict[str, Any]:
    selection = load_selection(root)
    fit_path = require_file(root, "fit_module")
    protocol_path = require_file(root, "protocol_module")
    runner = require_file(root, "runner")
    fit = load_module(fit_path, "rfam_bound_dense_fit")
    protocol = load_module(protocol_path, "rfam_bound_dense_protocol")

    access_order = ["selection", "training-reference"]
    training_document = load_json(require_file(root, "training_reference"), "training-reference")
    training_cases, _ = fit.validate_reference(
        reference_view(training_document, "protocol-training"),
        "protocol-training",
        TRAINING_IDS,
    )
    training = evaluate_partition(runner, fit, protocol, training_cases)

    # Holdout is intentionally first opened only after the frozen vector has
    # completed its dense training evaluation.  No selection step exists here.
    access_order.append("holdout-reference")
    holdout_document = load_json(require_file(root, "holdout_reference"), "holdout-reference")
    holdout_cases, _ = fit.validate_reference(
        reference_view(holdout_document, "protocol-holdout"),
        "protocol-holdout",
        HOLDOUT_IDS,
    )
    holdout = evaluate_partition(runner, fit, protocol, holdout_cases)
    return {
        "schema": "rfam.r2-c3mech-dense-comparison/v1",
        "status": "blocked-dense-detailed-reference-route-comparison-only",
        "source": {
            "comparator_sha256": sha256(Path(__file__).resolve()),
            "selection_sha256": ARTIFACTS["selection"][1],
            "training_reference_sha256": ARTIFACTS["training_reference"][1],
            "holdout_reference_sha256": ARTIFACTS["holdout_reference"][1],
            "runner_sha256": ARTIFACTS["runner"][1],
            "fit_module_sha256": ARTIFACTS["fit_module"][1],
            "protocol_module_sha256": ARTIFACTS["protocol_module"][1],
        },
        "selection": {
            "selection_split": selection["selection_split"],
            "holdout_influenced_selection": selection["holdout_influenced_selection"],
            "parameter_names": list(PARAMETER_NAMES),
            "parameter_vector": list(PARAMETER_VECTOR),
            "coarse_training_objective": selection["coarse_training_objective"],
            "parameter_vector_publishable": False,
        },
        "access_order": access_order,
        "training": training,
        "holdout": holdout,
        "approvals": {
            "mechanism_selected_for_fit": False,
            "experimental_training_selected": False,
            "experimental_holdout_selected": False,
            "parameter_vector_publishable": False,
            "physical_validation_approved": False,
            "production_ready": False,
        },
        "limitations": [
            "C3Mech is a detailed-mechanism candidate, not experimental truth",
            "the parameter vector was selected on the earlier coarse training grid",
            "dense holdout is evaluated once and cannot change the vector",
            "the 2000 K mechanism thermo boundary remains fail-closed",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        document = execute(args.workspace_root.resolve())
        publish_no_replace(args.output.resolve(), canonical_bytes(document))
        print(
            "R2_C3MECH_DENSE_COMPARISON_BLOCKED training_objective="
            + format(document["training"]["objective"], ".12g")
            + " holdout_objective="
            + format(document["holdout"]["objective"], ".12g")
        )
        return 0
    except Exception as error:
        print(
            "R2_C3MECH_DENSE_COMPARISON_FAIL "
            + type(error).__name__
            + ":"
            + str(error),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
