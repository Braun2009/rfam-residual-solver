#!/usr/bin/env python3
"""Generate a dense, fixed-grid C3Mech R2 calibration reference.

This is deliberately a separate route from the frozen six-case screen.  The
time grid is fixed in a checked-in manifest before Cantera is run; no time
point is selected from a RFAM result or from a holdout objective.  The tool
reuses the reviewed C3Mech trajectory implementation only after binding its
canonical source bytes and all external authorities.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Any


class DenseReferenceError(RuntimeError):
    """A fail-closed provenance, manifest, runtime, or generation error."""


ROOT = Path(__file__).resolve().parents[1]
BASE_GENERATOR_REL = Path("realFluidAllMach/scripts/generate_r2_c3mech_candidate_reference.py")
BASE_GENERATOR_SHA256 = "d8d4d3ad9cc4207039c49f9d2174d8e28ca8db601bea077a85d955d11d55bb5c"
BASE_GENERATOR_SIZE = 22366
MECHANISM_REL = Path(
    "tmp/c3mech-v4.0.1-rights-screen/C3MechV4.0.1_44I0_C0-C3-C4_LT-HT_Cantera.yaml"
)
MECHANISM_SHA256 = "d5cec179ffbb6221b1b91264e1439ec850658dd97a0298054c8dbc336f3ab148"
MECHANISM_SIZE = 2227356
LICENSE_REL = Path("tmp/c3mech-v4.0.1-rights-screen/LICENSE")
LICENSE_SHA256 = "7e7170e3cebf88a9f60c7b8421418323c09304da1af4d5e90f4da1dc1c8a2661"
LICENSE_SIZE = 18650
RIGHTS_REPORT_REL = Path("tmp/a408-c3mech-v4.0.1-propane-mechanism-rights-screen-v1-20260728.md")
RIGHTS_REPORT_SHA256 = "6779ae447b8ae27d49f4b1011e91ab0747dc8c91435f44d5effcdf6cad5fd26f"
RIGHTS_REPORT_SIZE = 9722
MANIFESTS = {
    "training": (
        Path("realFluidAllMach/tests/reference/r2_c3mech_dense_training_cases_v1.json"),
        "de56596277f253b748cf4fda168e64489ca47436f9045b1f7d22ad0f46649cdf",
        1213,
        "protocol-training",
        ("train_phi0p8_diluted", "train_phi1p0", "train_phi1p2"),
    ),
    "holdout": (
        Path("realFluidAllMach/tests/reference/r2_c3mech_dense_holdout_cases_v1.json"),
        "92496aeb3fdd38c6d5006919249ebf5698058b93595515c92a5edc204a2d7a6e",
        1217,
        "protocol-holdout",
        ("holdout_phi0p6_diluted", "holdout_phi1p1", "holdout_phi1p4"),
    ),
}

SAMPLE_TIMES = (
    0.0,
    1e-8,
    1.5e-8,
    2e-8,
    3e-8,
    4e-8,
    6e-8,
    8e-8,
    1e-7,
    1.5e-7,
    2e-7,
    3e-7,
    4e-7,
    6e-7,
    8e-7,
    1e-6,
    1.5e-6,
    2e-6,
    3e-6,
    4e-6,
    6e-6,
    8e-6,
    1e-5,
    1.5e-5,
    2e-5,
    3e-5,
    4e-5,
    6e-5,
    8e-5,
    1e-4,
    1.5e-4,
    2e-4,
    3e-4,
    4e-4,
    6e-4,
    8e-4,
    1e-3,
    1.5e-3,
    2e-3,
    3e-3,
    4e-3,
    5e-3,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fail_constant(value: str) -> None:
    raise DenseReferenceError("nonfinite-json-constant:" + value)


def strict_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DenseReferenceError("duplicate-json-key:" + key)
        result[key] = value
    return result


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=strict_pairs,
            parse_constant=fail_constant,
        )
    except DenseReferenceError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise DenseReferenceError("invalid-json:" + label) from error
    if not isinstance(value, dict):
        raise DenseReferenceError("json-root-not-object:" + label)
    return value


def require_file(root: Path, relative: Path, expected_sha: str, expected_size: int, label: str) -> Path:
    declared = root / relative
    if declared.is_symlink() or not declared.is_file():
        raise DenseReferenceError("missing-or-aliased-authority:" + label)
    path = declared.resolve()
    expected = (root / relative).resolve()
    if path != expected:
        raise DenseReferenceError("authority-path-mismatch:" + label)
    if path.stat().st_size != expected_size or sha256(path) != expected_sha:
        raise DenseReferenceError("authority-mismatch:" + label)
    return path


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DenseReferenceError("non-numeric:" + label)
    result = float(value)
    if not math.isfinite(result):
        raise DenseReferenceError("nonfinite:" + label)
    return result


def validate_manifest(root: Path, partition: str) -> tuple[Path, dict[str, Any]]:
    try:
        relative, expected_sha, expected_size, expected_partition, expected_ids = MANIFESTS[partition]
    except KeyError as error:
        raise DenseReferenceError("unknown-partition:" + partition) from error
    path = require_file(root, relative, expected_sha, expected_size, "case-manifest")
    document = load_json(path, "case-manifest")
    if set(document) != {"schema", "partition", "partition_locked", "sample_times_s", "cases"}:
        raise DenseReferenceError("case-manifest-schema-open")
    if document["schema"] != "rfam.r2-c3mech-case-partition/v1":
        raise DenseReferenceError("case-manifest-schema-mismatch")
    if document["partition"] != expected_partition or document["partition_locked"] is not True:
        raise DenseReferenceError("case-manifest-partition-mismatch")
    times = document["sample_times_s"]
    if not isinstance(times, list) or tuple(finite(value, "sample-time") for value in times) != SAMPLE_TIMES:
        raise DenseReferenceError("sample-grid-mismatch")
    cases = document["cases"]
    if not isinstance(cases, list) or tuple(case.get("id") for case in cases if isinstance(case, dict)) != expected_ids:
        raise DenseReferenceError("case-order-mismatch")
    allowed = {"id", "temperature_k", "pressure_pa", "equivalence_ratio", "fuel", "oxidizer"}
    for case in cases:
        if not isinstance(case, dict) or set(case) != allowed:
            raise DenseReferenceError("case-schema-mismatch")
        if finite(case["temperature_k"], "case-temperature") <= 0.0:
            raise DenseReferenceError("case-temperature-invalid")
        if finite(case["pressure_pa"], "case-pressure") <= 0.0:
            raise DenseReferenceError("case-pressure-invalid")
        if finite(case["equivalence_ratio"], "case-equivalence-ratio") <= 0.0:
            raise DenseReferenceError("case-equivalence-ratio-invalid")
        if not isinstance(case["fuel"], str) or not isinstance(case["oxidizer"], str):
            raise DenseReferenceError("case-composition-invalid")
    return path, document


def load_base_generator(root: Path):
    path = require_file(root, BASE_GENERATOR_REL, BASE_GENERATOR_SHA256, BASE_GENERATOR_SIZE, "base-generator")
    # The base module imports Cantera at import time.  Make the already bound
    # installed tree visible without allowing an arbitrary package shadow.
    package_root = Path(r"F:\tmp\cantera-3.2.0-python312")
    if str(package_root) not in sys.path:
        sys.path.insert(0, str(package_root))
    spec = importlib.util.spec_from_file_location("rfam_bound_c3mech_generator", path)
    if spec is None or spec.loader is None:
        raise DenseReferenceError("base-generator-import-spec-failed")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as error:
        raise DenseReferenceError("base-generator-import-failed") from error
    if sha256(path) != BASE_GENERATOR_SHA256:
        raise DenseReferenceError("base-generator-changed-after-import")
    return module


def generate(root: Path, partition: str) -> dict[str, Any]:
    manifest_path, manifest = validate_manifest(root, partition)
    mechanism = require_file(root, MECHANISM_REL, MECHANISM_SHA256, MECHANISM_SIZE, "mechanism")
    license_path = require_file(root, LICENSE_REL, LICENSE_SHA256, LICENSE_SIZE, "license")
    rights_report = require_file(root, RIGHTS_REPORT_REL, RIGHTS_REPORT_SHA256, RIGHTS_REPORT_SIZE, "rights-report")
    base = load_base_generator(root)
    try:
        runtime = base.verify_runtime()
        _, inventory = base.load_mechanism(mechanism, "44I0")
        if hasattr(base.ct, "suppress_thermo_warnings"):
            base.ct.suppress_thermo_warnings()
        cases = []
        for case in manifest["cases"]:
            # clone=False makes each reactor own the supplied phase.  Reload a
            # fresh, independently verified phase per case so a prior ReactorNet
            # lifetime cannot make trajectory results depend on case order.
            case_gas, case_inventory = base.load_mechanism(mechanism, "44I0")
            if case_inventory != inventory:
                raise DenseReferenceError("case-mechanism-inventory-drift")
            cases.append(
                base.run_case(
                    case_gas,
                    case,
                    list(SAMPLE_TIMES),
                    float(inventory["thermo_min_k"]),
                    float(inventory["thermo_max_k"]),
                )
            )
    except Exception as error:
        raise DenseReferenceError("cantera-generation-failed") from error
    any_exit = any(
        case["thermo_envelope_exit"]["status"]
        == "thermo-envelope-exit-detected-no-endpoint-state-published"
        for case in cases
    )
    return {
        "schema": "rfam.r2-c3mech-dense-calibration-reference/v1",
        "status": "OUT_OF_MECHANISM_THERMO_ENVELOPE" if any_exit else "WITHIN_MECHANISM_THERMO_ENVELOPE",
        "model": inventory,
        "partition": manifest["partition"],
        "runtime_identity": runtime,
        "source": {
            "mechanism_path": str(mechanism),
            "mechanism_sha256": sha256(mechanism),
            "mechanism_size": mechanism.stat().st_size,
            "license_path": str(license_path),
            "license_sha256": sha256(license_path),
            "license_size": license_path.stat().st_size,
            "rights_report_path": str(rights_report),
            "rights_report_sha256": sha256(rights_report),
            "rights_report_size": rights_report.stat().st_size,
            "license_declaration": "CC-BY-4.0-source-screen-candidate-not-fit-approval",
            "base_generator_path": str((root / BASE_GENERATOR_REL).resolve()),
            "base_generator_sha256": BASE_GENERATOR_SHA256,
            "case_manifest_path": str(manifest_path),
            "case_manifest_sha256": sha256(manifest_path),
            "case_manifest_size": manifest_path.stat().st_size,
            "dense_generator_path": str(Path(__file__).resolve()),
            "dense_generator_sha256": sha256(Path(__file__).resolve()),
        },
        "reactor": {
            "model": "IdealGasConstPressureReactor",
            "phase": "ideal-gas",
            "energy": "on",
            "rtol": 1.0e-10,
            "atol": 1.0e-18,
            "temperature_advance_limit_k": 10.0,
            "ignition_temperature_advance_limit_k": 1.0,
            "outside_envelope_policy": "stop-at-first-detected-exit; publish-no-state; leave-later-samples-unevaluated",
        },
        "sampling_policy": {
            "grid_kind": "fixed-predeclared-log-refinement",
            "grid_source": "checked-in-case-manifest",
            "sample_count": len(SAMPLE_TIMES),
            "candidate_or_reference_adaptive": False,
            "holdout_grid_reused_after_selection": True,
            "grid_not_derived_from_rfam_or_fit": True,
        },
        "sample_times_s": list(SAMPLE_TIMES),
        "major_species": list(base.SPECIES),
        "prefix_diagnostic": {
            "eligible_samples_only": True,
            "outside_endpoint_state_published": False,
            "later_samples_evaluated": False,
        },
        "cases": cases,
        "approvals": {
            "mechanism_selected_for_fit": False,
            "fit_executed": False,
            "holdout_used_for_selection": False,
            "physical_validation_approved": False,
            "production_ready": False,
        },
    }


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (json.dumps(document, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")


def publish_no_replace(path: Path, raw: bytes) -> None:
    if os.path.lexists(path):
        raise DenseReferenceError("output-path-already-exists")
    if not path.parent.is_dir():
        raise DenseReferenceError("output-parent-missing")
    fd, temporary_name = tempfile.mkstemp(prefix=".r2-dense-", suffix=".tmp", dir=str(path.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise DenseReferenceError("output-publication-race") from error
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=ROOT.parent)
    parser.add_argument("--partition", choices=sorted(MANIFESTS), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        root = args.workspace_root.resolve()
        output = args.output.resolve()
        document = generate(root, args.partition)
        publish_no_replace(output, canonical_bytes(document))
        print(
            "R2_C3MECH_DENSE_REFERENCE_BLOCKED partition="
            + args.partition
            + " cases="
            + str(len(document["cases"]))
            + " samples="
            + str(document["sampling_policy"]["sample_count"])
        )
        return 0
    except (DenseReferenceError, OSError, ValueError) as error:
        print("R2_C3MECH_DENSE_REFERENCE_FAIL " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
