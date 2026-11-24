#!/usr/bin/env python3
"""Fail-closed locator and regeneration gate for the internal UCSD R2 screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


class ProvenanceError(RuntimeError):
    pass


ARTIFACTS = [
    ("ucsd_chemkin_mechanism", "Vortex/tools/openfoam/sandiego20161214_mechCK.txt", "728aa95e372aa71c59528850ac35511bf4c5f457954b9803c6ef823996163c73", 24206, "official-chemkin-mechanism-local-copy"),
    ("ucsd_chemkin_thermo", "Vortex/tools/openfoam/sandiego20160815_therm.txt", "ba9518b7007d007353ca6ac927f881b05d6b302bbb56bc32c9e1ee80d3333d92", 34813, "official-chemkin-thermo-local-copy"),
    ("ucsd_chemkin_transport", "Vortex/tools/openfoam/sandiego20160815_trans.txt", "dc8134f2c62cd54c7999d2428ac1de22743e2e331caef2da8286763bedd5ce94", 75848, "official-chemkin-transport-local-copy"),
    ("ucsd_converted_yaml", "Vortex/tools/openfoam/sandiego2016_c3h8.yaml", "fbdbf6ce0bdbe85448229395f42096b59eb08828e92416d0f026bb6ae5040eaa", 66674, "cantera-converted-mechanism-candidate"),
    ("conversion_provenance_contract", "realFluidAllMach/tests/reference/r2_physics_candidate_v1.json", "cdd6046bab363ca5c9b0c637d341b9ff4ca3501b6e0f9dd142fdca1d17ba59d5", 14673, "blocked-conversion-provenance-contract"),
    ("cantera_reference_generator", "realFluidAllMach/scripts/generate_r2_0d_calibration_reference.py", "43ad404b9eb63a1fb6dd343023a6f726cb5efd4f3e84229c9a50336fb460ed49", 9428, "cantera-0d-reference-generator"),
    ("detailed_reference", "tmp/a408-r2-ucsd-internal-screen-v1.json", "ffcf5fb0c88b1cf2a4be465bfe2b53b9b306012154fc74c02b45a2b5a1c3984f", 28443, "generated-internal-detailed-trajectory"),
    ("comparison_tool", "realFluidAllMach/tools/screen_r2_detailed_candidate.py", "6e88e92ebe84db8ae378f95399de4764ca55f122ef9156ffb5c42e48874a14e7", 12960, "internal-screen-comparator"),
    ("comparison_protocol_module", "realFluidAllMach/tools/run_r2_0d_rfam_protocol.py", "5a4642c37ce4b07ca97db82a6af218d39e039b8a73e417d20efb6094d43f58a8", 20748, "internal-screen-transitive-protocol-module"),
    ("rfam_runner_source", "realFluidAllMach/tools/r2_0d_rfam_runner.cpp", "2f6acfe72015b4de91fe1ade806ea8f6d30bf8d9049fb7d89dc936e546eeabbb", 27822, "rfam-two-step-runner-source"),
    ("rfam_runner_executable", "realFluidAllMach/build-phase-foundation/rfam_r2_0d_runner", "4858ae43b17bec7cccf2f6ba3a75659f0e92ed092635be6d3b61b107e1d8a1d1", 93312, "frozen-linux-elf-behavior"),
    ("comparison_output", "tmp/a408-r2-ucsd-vs-rfam-screen-v1.json", "feaf8f280d4f92ef1d37d480a4d0e8bb2e32a47f48f2a2b201b5b9af55ac21b5", 5241, "compact-internal-screen-result"),
]

APPROVALS = {
    "mechanism_selected_for_fit": False,
    "experimental_data_selected": False,
    "r1_r2_fitting_approved": False,
    "physical_validation_approved": False,
    "production_ready": False,
    "reacting_hot_start_approved": False,
}

RUNTIME_IDENTITY = {
    "python_executable": "/usr/bin/python3.8",
    "python_sha256": "298a9e830ed52f36c299427565485d717d1ce0179c0597cc16560513eb780b06",
    "cantera_version": "3.1.0",
    "cantera_module": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/__init__.py",
    "cantera_module_sha256": "9a7ceae3bae1a7e41dc581e353d44e04af8447114748f5b7c75995e26f8819b9",
    "cantera_extension": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/_cantera.cpython-38-x86_64-linux-gnu.so",
    "cantera_extension_sha256": "c57a47c201b0cb8f965b38acc7ebfafbbe95db7464d80359499c64d50555fae8",
}

EXPECTED = {
    "schema": "rfam.r2-ucsd-internal-screen-provenance/v1",
    "contract_version": 1,
    "status": "blocked-internal-screening-only",
    "scope": "provenance-recovery-only-no-fit-no-physical-validation",
    "approvals": APPROVALS,
    "artifacts": [
        {"id": artifact_id, "path": path, "sha256": digest, "size": size, "kind": kind}
        for artifact_id, path, digest, size, kind in ARTIFACTS
    ],
    "runtime_identity": RUNTIME_IDENTITY,
    "reference_binding": {
        "mechanism_artifact_id": "ucsd_converted_yaml",
        "generator_artifact_id": "cantera_reference_generator",
        "output_artifact_id": "detailed_reference",
        "schema": "rfam.r2-0d-cantera-reference/v1",
        "status": "local-cantera-reference-not-fit-eligible",
        "case_count": 6,
        "major_species": ["C3H8", "O2", "CO", "CO2", "H2O", "N2"],
    },
    "comparison_binding": {
        "tool_artifact_id": "comparison_tool",
        "protocol_module_artifact_id": "comparison_protocol_module",
        "runner_source_artifact_id": "rfam_runner_source",
        "runner_artifact_id": "rfam_runner_executable",
        "reference_artifact_id": "detailed_reference",
        "output_artifact_id": "comparison_output",
        "schema": "rfam.r2-detailed-candidate-screen/v1",
        "status": "internal-screening-only-no-fit-no-physical-validation",
    },
    "conversion_rights": {
        "provenance_contract_artifact_id": "conversion_provenance_contract",
        "mechanism_status": "candidate-only",
        "license_status": "unverified",
        "fit_use_allowed": False,
        "redistribution_rights_clear": False,
    },
    "limitations": [
        "internal deterministic route screen only",
        "no mechanism selection or fit authority",
        "no experimental training or holdout authority",
        "no physical validation or production chemistry authority",
    ],
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ProvenanceError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=strict_object)
    except ProvenanceError:
        raise
    except Exception as exc:
        raise ProvenanceError(f"invalid-json:{label}:{exc}") from exc


def exact_typed(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise ProvenanceError(f"type-mismatch:{label}")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ProvenanceError(f"key-mismatch:{label}")
        for key, value in expected.items():
            exact_typed(actual[key], value, f"{label}.{key}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ProvenanceError(f"length-mismatch:{label}")
        for index, value in enumerate(expected):
            exact_typed(actual[index], value, f"{label}[{index}]")
    elif actual != expected:
        raise ProvenanceError(f"value-mismatch:{label}")


def resolve_artifact(root: Path, relative: str) -> Path:
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise ProvenanceError(f"invalid-artifact-path:{relative}")
    path = (root / relative_path).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ProvenanceError(f"artifact-path-escapes-root:{relative}") from exc
    return path


def verify_artifacts(root: Path) -> dict[str, Path]:
    resolved: dict[str, Path] = {}
    for artifact_id, relative, digest, size, _ in ARTIFACTS:
        path = resolve_artifact(root, relative)
        if not path.is_file():
            raise ProvenanceError(f"artifact-missing:{artifact_id}")
        if path.stat().st_size != size:
            raise ProvenanceError(f"artifact-size-mismatch:{artifact_id}")
        if sha256(path) != digest:
            raise ProvenanceError(f"artifact-sha256-mismatch:{artifact_id}")
        resolved[artifact_id] = path
    return resolved


def verify_embedded_identity(paths: dict[str, Path]) -> None:
    reference = load_json(paths["detailed_reference"], "detailed-reference")
    expected_reference = EXPECTED["reference_binding"]
    if reference.get("schema") != expected_reference["schema"]:
        raise ProvenanceError("reference-schema-mismatch")
    if reference.get("status") != expected_reference["status"]:
        raise ProvenanceError("reference-status-mismatch")
    if reference.get("generator_sha256") != ARTIFACTS[5][2]:
        raise ProvenanceError("reference-generator-sha256-mismatch")
    if reference.get("mechanism_sha256") != ARTIFACTS[3][2]:
        raise ProvenanceError("reference-mechanism-sha256-mismatch")
    if reference.get("cantera_version") != RUNTIME_IDENTITY["cantera_version"]:
        raise ProvenanceError("reference-cantera-version-mismatch")
    expected_runtime = dict(RUNTIME_IDENTITY)
    expected_runtime.pop("cantera_version")
    exact_typed(reference.get("runtime_identity"), expected_runtime, "reference.runtime_identity")
    if len(reference.get("cases", [])) != expected_reference["case_count"]:
        raise ProvenanceError("reference-case-count-mismatch")
    exact_typed(reference.get("major_species"), expected_reference["major_species"], "reference.major_species")

    comparison = load_json(paths["comparison_output"], "comparison-output")
    expected_comparison = EXPECTED["comparison_binding"]
    if comparison.get("schema") != expected_comparison["schema"]:
        raise ProvenanceError("comparison-schema-mismatch")
    if comparison.get("status") != expected_comparison["status"]:
        raise ProvenanceError("comparison-status-mismatch")
    exact_typed(comparison.get("approvals"), {
        "experimental_data_selected": False,
        "mechanism_selected_for_fit": False,
        "physical_validation_approved": False,
        "production_ready": False,
        "r1_r2_fitting_approved": False,
    }, "comparison.approvals")
    source = comparison.get("source", {})
    expected_source = {
        "detailed_reference_path": "/mnt/f/tmp/a408-r2-ucsd-internal-screen-v1.json",
        "detailed_reference_sha256": ARTIFACTS[6][2],
        "detailed_mechanism_sha256": ARTIFACTS[3][2],
        "detailed_generator_sha256": ARTIFACTS[5][2],
        "cantera_version": RUNTIME_IDENTITY["cantera_version"],
        "rfam_runner_path": "/mnt/f/realFluidAllMach/build-phase-foundation/rfam_r2_0d_runner",
        "rfam_runner_sha256": ARTIFACTS[10][2],
    }
    exact_typed(source, expected_source, "comparison.source")


def to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive
    if len(drive) != 2 or drive[1] != ":":
        raise ProvenanceError(f"unsupported-wsl-path:{resolved}")
    relative = resolved.relative_to(Path(drive + "\\"))
    return f"/mnt/{drive[0].lower()}/{relative.as_posix()}"


def run_checked(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as exc:
        raise ProvenanceError(f"{label}-launch-failed:{exc}") from exc
    if result.returncode != 0:
        detail = result.stderr.strip().replace("\n", " | ")
        raise ProvenanceError(f"{label}-failed:{result.returncode}:{detail}")
    return result


def verify_reference_regeneration(root: Path, paths: dict[str, Path]) -> None:
    temp_parent = root / "tmp"
    with tempfile.TemporaryDirectory(prefix="rfam-ucsd-reference-", dir=temp_parent) as directory:
        temp = Path(directory)
        generator = temp / "generate_reference.py"
        mechanism = temp / "mechanism.yaml"
        output = temp / "reference.json"
        shutil.copyfile(paths["cantera_reference_generator"], generator)
        shutil.copyfile(paths["ucsd_converted_yaml"], mechanism)
        run_checked([
            "wsl", "-e", "/usr/bin/python3.8", to_wsl(generator),
            "--mechanism", to_wsl(mechanism), "--output", to_wsl(output),
        ], "reference-regeneration")
        if output.read_bytes() != paths["detailed_reference"].read_bytes():
            raise ProvenanceError("reference-regeneration-byte-mismatch")


def verify_comparison_regeneration(root: Path, paths: dict[str, Path]) -> None:
    temp_parent = root / "tmp"
    with tempfile.TemporaryDirectory(prefix="rfam-ucsd-comparison-", dir=temp_parent) as directory:
        temp = Path(directory)
        comparator = temp / "screen.py"
        protocol_module = temp / "run_r2_0d_rfam_protocol.py"
        output = temp / "comparison.json"
        shutil.copyfile(paths["comparison_tool"], comparator)
        shutil.copyfile(paths["comparison_protocol_module"], protocol_module)
        run_checked([
            "wsl", "-e", "/usr/bin/python3.8", to_wsl(comparator),
            "--detailed-reference", "/mnt/f/tmp/a408-r2-ucsd-internal-screen-v1.json",
            "--rfam-runner", "/mnt/f/realFluidAllMach/build-phase-foundation/rfam_r2_0d_runner",
            "--output", to_wsl(output),
        ], "comparison-regeneration")
        if output.read_bytes() != paths["comparison_output"].read_bytes():
            raise ProvenanceError("comparison-regeneration-byte-mismatch")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--contract", type=Path, default=Path("realFluidAllMach/tests/reference/r2_ucsd_internal_screen_provenance_v1.json"))
    parser.add_argument("--verify-reference-regeneration", action="store_true")
    parser.add_argument("--verify-comparison-regeneration", action="store_true")
    args = parser.parse_args(argv)

    try:
        root = args.root.resolve()
        contract = args.contract if args.contract.is_absolute() else root / args.contract
        document = load_json(contract.resolve(), "contract")
        exact_typed(document, EXPECTED, "contract")
        paths = verify_artifacts(root)
        verify_embedded_identity(paths)
        checks = len(paths) + 3
        if args.verify_reference_regeneration:
            verify_reference_regeneration(root, paths)
            checks += 1
        if args.verify_comparison_regeneration:
            verify_comparison_regeneration(root, paths)
            checks += 1
    except Exception as exc:
        print(f"R2_UCSD_INTERNAL_SCREEN_PROVENANCE_FAIL reason={exc}", file=sys.stderr)
        return 2

    print(f"R2_UCSD_INTERNAL_SCREEN_PROVENANCE_BLOCKED checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
