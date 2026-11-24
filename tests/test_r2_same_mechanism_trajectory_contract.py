#!/usr/bin/env python3
"""Focused regressions for the R2 same-mechanism trajectory contract."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


class RegressionFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RegressionFailure(message)


def run_validator(root: Path, contract: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(root / "tools/validate_r2_same_mechanism_trajectory_contract.py"),
            "--root",
            str(root),
            "--contract",
            str(contract),
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def assert_rejected(result: subprocess.CompletedProcess[str], label: str) -> None:
    combined = result.stdout + result.stderr
    require(result.returncode == 2, f"{label}-return-code")
    require("R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_FAIL" in combined, f"{label}-fail-marker")
    require("R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_BLOCKED" not in combined, f"{label}-blocked-marker")
    require("Traceback" not in combined, f"{label}-traceback")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    canonical_path = root / "tests/reference/r2_same_mechanism_trajectory_contract_v1.json"
    contract = json.loads(canonical_path.read_text(encoding="utf-8"))
    checks = 0
    try:
        canonical = run_validator(root, canonical_path)
        require(canonical.returncode == 0, "canonical-return-code")
        require("R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_BLOCKED" in canonical.stdout, "canonical-marker")
        require("CONTRACT_FAIL" not in canonical.stdout + canonical.stderr, "canonical-fail-marker")
        checks += 3

        mutations: list[tuple[str, Callable[[dict[str, Any]], None]]] = [
            ("approval", lambda value: value["approvals"].__setitem__("production_ready", True)),
            ("guard", lambda value: value["thermo_envelope"].__setitem__("comparison_maximum_temperature_k", 3500.0)),
            ("reference-policy", lambda value: value["energy_reference"].__setitem__("rfam_reference_policy", "cantera-polynomial-extrapolation")),
            ("tail-species", lambda value: value["energy_reference"].__setitem__("lower_linear_tail_species", ["C3H8"])),
            ("tolerance", lambda value: value["acceptance_tolerances"].__setitem__("trajectory_temperature_absolute_k", 1.0)),
            ("count", lambda value: value["expected_counts"].__setitem__("eligible_samples", 732)),
            ("case", lambda value: value["cases"][0].__setitem__("temperature_k", 1201.0)),
            ("artifact", lambda value: value["artifacts"][0].__setitem__("path", "tests/reference/r2_instantaneous_reference_v1.json")),
        ]
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            for label, mutate in mutations:
                candidate = copy.deepcopy(contract)
                mutate(candidate)
                path = temporary / f"{label}.json"
                path.write_text(json.dumps(candidate), encoding="utf-8")
                assert_rejected(run_validator(root, path), label)
                checks += 4

            duplicate = temporary / "duplicate.json"
            raw = canonical_path.read_text(encoding="utf-8")
            duplicate.write_text(
                raw.replace(
                    '"production_ready": false,',
                    '"production_ready": true,\n    "production_ready": false,',
                    1,
                ),
                encoding="utf-8",
            )
            assert_rejected(run_validator(root, duplicate), "duplicate")
            checks += 4
    except (RegressionFailure, OSError, ValueError) as error:
        print(f"R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_TEST_FAIL {error}", file=sys.stderr)
        return 1
    except Exception as error:
        print(
            "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_TEST_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 1

    print(f"R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_TEST_PASS checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
