#!/usr/bin/env python3
"""Focused regressions for the low-concentration regularization contract."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class RegressionFailure(RuntimeError):
    pass


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RegressionFailure(label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_validator(
    validator: Path,
    root: Path,
    contract: Path,
    executable: Path,
    cwd: Path,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(validator),
            "--root",
            str(root),
            "--contract",
            str(contract),
            "--executable",
            str(executable),
        ],
        cwd=cwd,
        check=False,
        text=True,
        capture_output=True,
    )


def require_rejected(result: subprocess.CompletedProcess[str], label: str) -> None:
    combined = result.stdout + result.stderr
    require(result.returncode == 2, f"{label}-return-code")
    require(
        "R2_LOW_CONCENTRATION_REGULARIZATION_FAIL" in combined,
        f"{label}-fail-marker",
    )
    require(
        "R2_LOW_CONCENTRATION_REGULARIZATION_BLOCKED" not in combined,
        f"{label}-blocked-marker",
    )
    require("Traceback" not in combined, f"{label}-traceback")


def copy_snapshot(source_root: Path, destination: Path) -> tuple[Path, Path]:
    contract_source = (
        source_root
        / "tests/reference/r2_low_concentration_regularization_contract_v1.json"
    )
    contract = json.loads(contract_source.read_text(encoding="utf-8"))
    for artifact in contract["artifacts"]:
        source = source_root / artifact["path"]
        target = destination / artifact["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    contract_target = destination / contract_source.relative_to(source_root)
    contract_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(contract_source, contract_target)
    executable = destination / contract["artifacts"][-1]["path"]
    return contract_target, executable


def update_artifact(document: dict, artifact_id: str, path: Path) -> None:
    for artifact in document["artifacts"]:
        if artifact["id"] == artifact_id:
            artifact["size"] = path.stat().st_size
            artifact["sha256"] = sha256(path)
            return
    raise RegressionFailure(f"artifact-not-found:{artifact_id}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, default=None,
        help="Directory containing built probe artifacts; "
             "defaults to <root>/build-phase-foundation")
    args = parser.parse_args()
    root = args.root.resolve()
    build_dir = args.build_dir.resolve() if args.build_dir else root / "build-phase-foundation"
    validator = root / "tools/validate_r2_low_concentration_regularization.py"
    contract = (
        root / "tests/reference/r2_low_concentration_regularization_contract_v1.json"
    )
    executable = (
        build_dir / "rfam_r2_low_concentration_regularization_probe"
    )
    checks = 0
    try:
        require(root.is_dir(), "root-missing")
        require(validator.is_file(), "validator-missing")
        require(contract.is_file(), "contract-missing")
        require(executable.is_file(), "executable-missing")
        checks += 4

        with tempfile.TemporaryDirectory(prefix="rfam-r2-low-concentration-") as directory:
            temporary = Path(directory)
            arbitrary_cwd = temporary / "arbitrary-cwd"
            arbitrary_cwd.mkdir()
            canonical = run_validator(
                validator, root, contract, executable, arbitrary_cwd)
            combined = canonical.stdout + canonical.stderr
            require(canonical.returncode == 0, "canonical-return-code")
            require(
                "R2_LOW_CONCENTRATION_REGULARIZATION_BLOCKED checks=148"
                in canonical.stdout,
                "canonical-blocked-marker",
            )
            require("production_ready=0" in canonical.stdout, "canonical-production")
            require("REGULARIZATION_FAIL" not in combined, "canonical-fail-marker")
            require("Traceback" not in combined, "canonical-traceback")
            checks += 5

            module_spec = importlib.util.spec_from_file_location(
                "r2_low_concentration_validator", validator)
            require(module_spec is not None and module_spec.loader is not None,
                    "validator-import-spec")
            module = importlib.util.module_from_spec(module_spec)
            module_spec.loader.exec_module(module)
            for payload, label in (
                (b'{"approval":true,"approval":false}', "duplicate-key"),
                (b'{"production_\\u0072eady":true,"production_ready":false}',
                 "escaped-duplicate-key"),
            ):
                rejected = False
                try:
                    module.load_json_bytes(payload, label)
                except module.ContractError as error:
                    rejected = str(error).startswith("duplicate-json-key:")
                require(rejected, label)
                checks += 1

            invalid_number = subprocess.run(
                [str(executable), "nan", "1"],
                check=False,
                text=True,
                capture_output=True,
            )
            require(invalid_number.returncode == 2, "probe-nan-return-code")
            require(
                "R2_LOW_CONCENTRATION_PROBE_FAIL" in invalid_number.stderr,
                "probe-nan-fail-marker",
            )
            require(
                "R2_LOW_CONCENTRATION_PROBE_OK" not in invalid_number.stdout,
                "probe-nan-ok-marker",
            )
            checks += 3

            missing_argument = subprocess.run(
                [str(executable), "3400"],
                check=False,
                text=True,
                capture_output=True,
            )
            require(missing_argument.returncode == 2, "probe-argc-return-code")
            require(
                "R2_LOW_CONCENTRATION_PROBE_FAIL" in missing_argument.stderr,
                "probe-argc-fail-marker",
            )
            checks += 2

            for relative_path, label in (
                ("src/combustion/FiniteRate/MassActionThermo.H", "mass-action"),
                ("src/combustion/FiniteRate/RhsCompositionPreparation.H",
                 "rhs-composition-helper"),
                ("src/combustion/FiniteRate/RealFluidCVODEWrapper.H",
                 "rhs-production-call-site"),
                ("src/thermo/EquilibriumThermoModel.H", "equilibrium-owner"),
                ("src/variables/StateContract.H", "energy-contract"),
            ):
                snapshot = temporary / f"live-tamper-{label}"
                snapshot_contract, snapshot_executable = copy_snapshot(
                    root, snapshot)
                artifact = snapshot / relative_path
                artifact.write_bytes(artifact.read_bytes() + b"\n")
                tampered = run_validator(
                    validator,
                    snapshot,
                    snapshot_contract,
                    snapshot_executable,
                    arbitrary_cwd,
                )
                require_rejected(tampered, f"live-artifact-tamper-{label}")
                checks += 4

            snapshot = temporary / "self-consistent-rebind"
            snapshot_contract, snapshot_executable = copy_snapshot(root, snapshot)
            document = json.loads(snapshot_contract.read_text(encoding="utf-8"))
            alternate = snapshot / "tools/alternate_probe.cpp"
            alternate.parent.mkdir(parents=True, exist_ok=True)
            alternate.write_text("int main(){return 0;}\n", encoding="ascii")
            for artifact in document["artifacts"]:
                if artifact["id"] == "probe_source":
                    artifact["path"] = "tools/alternate_probe.cpp"
                    artifact["size"] = alternate.stat().st_size
                    artifact["sha256"] = sha256(alternate)
            snapshot_contract.write_text(
                json.dumps(document, indent=2) + "\n", encoding="utf-8")
            rebind = run_validator(
                validator,
                snapshot,
                snapshot_contract,
                snapshot_executable,
                arbitrary_cwd,
            )
            require_rejected(rebind, "self-consistent-rebind")
            checks += 4

            for mutation, label in (
                (("approvals", "production_chemistry_ready", True), "approval"),
                (("acceptance", "transition_derivative_relative_tolerance", 1.0),
                 "tolerance"),
                (("regularization", "transition_concentration_mol_m3", 1.0e-3),
                 "transition"),
                (("regularization", "state_smoothing", {
                    "source": "unbound-helper",
                    "below_transition_formula": "c0*x^2*(2-x), x=c/c0",
                    "at_or_above_transition": "identity",
                    "nonpositive": "zero",
                }), "state-source"),
            ):
                snapshot = temporary / f"contract-{label}"
                snapshot_contract, snapshot_executable = copy_snapshot(root, snapshot)
                document = json.loads(snapshot_contract.read_text(encoding="utf-8"))
                section, key, value = mutation
                document[section][key] = value
                snapshot_contract.write_text(
                    json.dumps(document, indent=2) + "\n", encoding="utf-8")
                result = run_validator(
                    validator,
                    snapshot,
                    snapshot_contract,
                    snapshot_executable,
                    arbitrary_cwd,
                )
                require_rejected(result, f"contract-{label}")
                checks += 4

            wrapper = temporary / "substitute-executable.sh"
            marker = temporary / "substitute-executed"
            wrapper.write_text(
                "#!/bin/sh\n"
                f"touch {str(marker)!r}\n"
                "exit 0\n",
                encoding="ascii",
            )
            wrapper.chmod(wrapper.stat().st_mode | 0o111)
            substitute = run_validator(
                validator, root, contract, wrapper, arbitrary_cwd)
            require_rejected(substitute, "substitute-executable")
            require(not marker.exists(), "substitute-executable-launched")
            checks += 5

            rebuild = temporary / "rfam_r2_low_concentration_rebuild"
            contract_document = json.loads(contract.read_text(encoding="utf-8"))
            command = list(contract_document["toolchain"]["compile_command"])
            command[-1] = str(rebuild)
            built = subprocess.run(
                command,
                cwd=root,
                check=False,
                text=True,
                capture_output=True,
            )
            require(built.returncode == 0, "rebuild-return-code")
            require(built.stderr == "", "rebuild-stderr")
            require(rebuild.is_file(), "rebuild-output-missing")
            require(sha256(rebuild) == sha256(executable), "rebuild-byte-drift")
            checks += 4
    except (RegressionFailure, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"R2_LOW_CONCENTRATION_REGULARIZATION_TEST_FAIL {error}", file=sys.stderr)
        return 1
    except Exception as error:
        print(
            "R2_LOW_CONCENTRATION_REGULARIZATION_TEST_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 1

    print(f"R2_LOW_CONCENTRATION_REGULARIZATION_TEST_PASS checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
