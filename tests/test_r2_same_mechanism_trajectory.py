#!/usr/bin/env python3
"""Focused regressions for the R2 same-mechanism trajectory comparator."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


class RegressionFailure(RuntimeError):
    pass


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RegressionFailure(label)


def run_comparator(
    root: Path,
    runner: Path,
    output: Path,
    comparator: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    comparator_path = comparator or (
        root / "tools/compare_r2_same_mechanism_trajectory.py")
    return subprocess.run(
        [
            sys.executable,
            str(comparator_path),
            "--root",
            str(root),
            "--contract",
            str(root / "tests/reference/r2_same_mechanism_trajectory_contract_v1.json"),
            "--runner",
            str(runner),
            "--output",
            str(output),
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def write_wrapper(path: Path, runner: Path, mutation: str) -> None:
    mutations = {
        "profile": (
            "stdout = stdout.replace("
            "'solver_tolerance_profile=strict-equivalence-v1', "
            "'solver_tolerance_profile=protocol-v1', 1)\n"
        ),
        "approval": (
            "stdout = stdout.replace('production_ready=0', "
            "'production_ready=1', 1)\n"
        ),
        "qdot": (
            "match = re.search(r'(qdot_w_m3=)([^|\\n]+)', stdout)\n"
            "if match:\n"
            "    changed = match.group(1) + format(float(match.group(2)) * 1.001, '.17g')\n"
            "    stdout = stdout[:match.start()] + changed + stdout[match.end():]\n"
        ),
        "truncate": (
            "lines = stdout.splitlines()\n"
            "for index in range(len(lines) - 1, -1, -1):\n"
            "    if lines[index].startswith('sample='):\n"
            "        del lines[index]\n"
            "        break\n"
            "stdout = '\\n'.join(lines) + '\\n'\n"
        ),
    }
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import re, subprocess, sys\n"
        f"runner = {str(runner)!r}\n"
        "result = subprocess.run([runner, *sys.argv[1:]], "
        "text=True, capture_output=True)\n"
        "stdout = result.stdout\n"
        + mutations[mutation]
        + "sys.stdout.write(stdout)\n"
        "sys.stderr.write(result.stderr)\n"
        "raise SystemExit(result.returncode)\n",
        encoding="ascii",
    )
    path.chmod(path.stat().st_mode | 0o111)


def write_comparator_mutation(
    path: Path,
    comparator: Path,
    mutation: str,
) -> None:
    replacements = {
        "trajectory-temperature-nan": (
            "float(reactor.T),\n"
            "            f\"cantera-trajectory:{case['id']}[{index}].temperature\",",
            "float(\"nan\"),\n"
            "            f\"cantera-trajectory:{case['id']}[{index}].temperature\",",
        ),
        "trajectory-composition-inf": (
            "float(reactor.thermo[name].Y[0]),\n"
            "                f\"cantera-trajectory:{case['id']}[{index}].Y[{name}]\",",
            "float(\"inf\"),\n"
            "                f\"cantera-trajectory:{case['id']}[{index}].Y[{name}]\",",
        ),
        "ignition-time-nan": (
            "float(network.step()), f\"cantera-ignition-time:{case['id']}\")",
            "float(\"nan\"), f\"cantera-ignition-time:{case['id']}\")",
        ),
        "derivative-qdot-inf": (
            "qdot = finite(rfam_qdot(gas, heats), \"cantera-derivative-qdot\")",
            "qdot = finite(float(\"inf\"), \"cantera-derivative-qdot\")",
        ),
        "derivative-species-nan": (
            "float(production[index]),\n"
            "                f\"cantera-net-production-rate:{species[index]}\",",
            "float(\"nan\"),\n"
            "                f\"cantera-net-production-rate:{species[index]}\",",
        ),
        "derivative-temperature-inf": (
            "temperature_rate = finite(\n"
            "        temperature_rate, \"cantera-temperature-derivative\")",
            "temperature_rate = finite(\n"
            "        float(\"inf\"), \"cantera-temperature-derivative\")",
        ),
    }
    source = comparator.read_text(encoding="utf-8")
    old, new = replacements[mutation]
    require(source.count(old) == 1, f"{mutation}-mutation-site")
    path.write_text(source.replace(old, new, 1), encoding="utf-8")


def require_rejected(
    result: subprocess.CompletedProcess[str],
    output: Path,
    label: str,
) -> None:
    combined = result.stdout + result.stderr
    require(result.returncode == 2, f"{label}-return-code")
    require("R2_SAME_MECHANISM_TRAJECTORY_FAIL" in combined, f"{label}-fail-marker")
    require("R2_SAME_MECHANISM_TRAJECTORY_PASS" not in combined, f"{label}-pass-marker")
    require("Traceback" not in combined, f"{label}-traceback")
    require(not output.exists(), f"{label}-output-published")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    runner = args.runner.resolve()
    checks = 0
    try:
        require(root.is_dir(), "root-missing")
        require(runner.is_file(), "runner-missing")
        checks += 2
        with tempfile.TemporaryDirectory(prefix="rfam-r2-trajectory-") as directory:
            temporary = Path(directory)
            output = temporary / "evidence.json"
            canonical = run_comparator(root, runner, output)
            require(canonical.returncode == 0, "canonical-return-code")
            require(
                "R2_SAME_MECHANISM_TRAJECTORY_PASS cases=6 samples=732 eligible=646"
                in canonical.stdout,
                "canonical-pass-marker",
            )
            require(
                "OUTSIDE_THERMO_ENVELOPE samples=86 scored=0 guard_samples=0"
                in canonical.stdout,
                "canonical-outside-marker",
            )
            require("R2_PHYSICS_GATE_BLOCKED" in canonical.stdout, "physics-block-marker")
            require("TRAJECTORY_FAIL" not in canonical.stdout + canonical.stderr,
                    "canonical-fail-marker")
            require(output.is_file(), "canonical-output-missing")
            checks += 6

            evidence = json.loads(output.read_text(encoding="ascii"))
            require(
                evidence["schema"]
                == "rfam.r2-same-mechanism-trajectory-evidence/v1",
                "evidence-schema",
            )
            require(
                evidence["status"]
                == "implementation-equivalence-pass-physical-validation-blocked",
                "evidence-status",
            )
            require(
                evidence["counts"]
                == {
                    "cases": 6,
                    "total_samples": 732,
                    "eligible_samples": 646,
                    "guard_samples": 0,
                    "outside_samples": 86,
                },
                "evidence-counts",
            )
            require(
                evidence["solver_profiles"]["rfam"]["name"]
                == "strict-equivalence-v1",
                "evidence-strict-profile",
            )
            require(
                evidence["energy_reference"]["lower_linear_tail_species"]
                == ["C3H8", "N2"],
                "evidence-lower-tail-species",
            )
            reference_construction = evidence["energy_reference_construction"]
            require(
                [item["species"] for item in reference_construction]
                == evidence["species_order"],
                "reference-construction-order",
            )
            require(
                [
                    item["species"]
                    for item in reference_construction
                    if item["branch"] == "lower-boundary-cp-linear-tail"
                ]
                == ["C3H8", "N2"],
                "reference-construction-lower-tail",
            )
            require(
                all(
                    item["rfam_minus_cantera_reference_j_per_kmol"] < 0.0
                    for item in reference_construction
                    if item["branch"] == "lower-boundary-cp-linear-tail"
                ),
                "reference-construction-offset-sign",
            )
            require(
                all(
                    item["rfam_minus_cantera_reference_j_per_kmol"] == 0.0
                    for item in reference_construction
                    if item["branch"] == "in-range-nasa7"
                ),
                "reference-construction-in-range-offset",
            )
            require(
                all(value is False for value in evidence["approvals"].values()),
                "evidence-approval-enabled",
            )
            metrics = evidence["maximum_metrics"]
            require(
                metrics["trajectory_temperature_absolute_k"] <= 2.0e-5,
                "temperature-tolerance",
            )
            require(
                metrics["trajectory_mass_fraction_absolute"] <= 2.0e-9,
                "composition-tolerance",
            )
            require(
                metrics["ignition_delay_absolute_s"] <= 2.0e-15,
                "ignition-tolerance",
            )
            require(
                metrics["qdot_scaled_error"] <= 2.0e-8
                and metrics["species_derivative_scaled_error"] <= 2.0e-8
                and metrics["temperature_derivative_scaled_error"] <= 2.0e-8,
                "derivative-tolerance",
            )
            require(
                all(math.isfinite(float(value)) for value in metrics.values()),
                "metric-nonfinite",
            )
            classifications = [
                sample
                for case in evidence["cases"]
                for sample in case["samples"]
            ]
            require(len(classifications) == 732, "audit-sample-count")
            require(
                sum(sample["classification"] == "ELIGIBLE"
                    for sample in classifications) == 646,
                "eligible-audit-count",
            )
            outside = [
                sample for sample in classifications
                if sample["classification"] == "OUTSIDE_THERMO_ENVELOPE"
            ]
            require(len(outside) == 86, "outside-audit-count")
            require(
                all(not sample["scored"] and sample["errors"] is None
                    for sample in outside),
                "outside-sample-scored",
            )
            require(
                all(
                    sample["scored"] and isinstance(sample["errors"], dict)
                    for sample in classifications
                    if sample["classification"] == "ELIGIBLE"
                ),
                "eligible-sample-not-scored",
            )
            checks += 19

            for mutation in ("profile", "approval", "qdot", "truncate"):
                wrapper = temporary / f"{mutation}_runner.py"
                write_wrapper(wrapper, runner, mutation)
                rejected_output = temporary / f"{mutation}.json"
                rejected = run_comparator(root, wrapper, rejected_output)
                require_rejected(rejected, rejected_output, mutation)
                checks += 5

            comparator = root / "tools/compare_r2_same_mechanism_trajectory.py"
            for mutation in (
                "trajectory-temperature-nan",
                "trajectory-composition-inf",
                "ignition-time-nan",
                "derivative-qdot-inf",
                "derivative-species-nan",
                "derivative-temperature-inf",
            ):
                mutated_comparator = temporary / f"{mutation}.py"
                write_comparator_mutation(
                    mutated_comparator, comparator, mutation)
                rejected_output = temporary / f"{mutation}.json"
                rejected = run_comparator(
                    root,
                    runner,
                    rejected_output,
                    comparator=mutated_comparator,
                )
                require_rejected(rejected, rejected_output, mutation)
                checks += 5

            occupied = temporary / "occupied.json"
            occupied.write_text("do-not-overwrite\n", encoding="ascii")
            occupied_result = run_comparator(root, runner, occupied)
            require_rejected(occupied_result, Path(str(occupied) + ".unused"), "occupied")
            require(
                occupied.read_text(encoding="ascii") == "do-not-overwrite\n",
                "occupied-output-overwritten",
            )
            checks += 1
    except (RegressionFailure, OSError, ValueError) as error:
        print(f"R2_SAME_MECHANISM_TRAJECTORY_TEST_FAIL {error}", file=sys.stderr)
        return 1
    except Exception as error:
        print(
            "R2_SAME_MECHANISM_TRAJECTORY_TEST_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 1

    print(f"R2_SAME_MECHANISM_TRAJECTORY_TEST_PASS checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
