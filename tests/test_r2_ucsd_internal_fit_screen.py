#!/usr/bin/env python3
"""Focused fail-close tests for the blocked UCSD internal fit screen."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
CONTRACT = (
    ROOT / "tests" / "reference" / "r2_ucsd_internal_fit_screen_contract_v1.json"
)
SCREEN_TOOL = TOOLS / "run_r2_ucsd_internal_fit_screen.py"

sys.path.insert(0, str(TOOLS))
import run_r2_ucsd_internal_fit_screen as screen  # noqa: E402
import split_r2_ucsd_internal_reference as splitter  # noqa: E402


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


class IdentitySpec:
    def physical_from_u(self, value: float) -> float:
        return float(value)


class R2UcsdInternalFitScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        screen.load_bound_helpers()

    def test_contract_hash_and_semantics_are_exact(self) -> None:
        raw = CONTRACT.read_bytes()
        self.assertEqual(sha256(raw), screen.EXPECTED_CONTRACT_SHA256)
        document, loaded_raw = screen.load_json(CONTRACT, "contract")
        self.assertEqual(loaded_raw, raw)
        specs, config, starts = screen.validate_contract(document)
        self.assertEqual([spec.name for spec in specs], screen.EXPECTED_PARAMETER_NAMES)
        self.assertEqual(config.maximum_evaluations, 800)
        self.assertEqual([kind for kind, _ in starts], [
            "midpoint", "halton-1", "halton-2", "halton-3", "halton-4",
        ])

        mutations = {
            "objective-policy": lambda value: value["objective"].__setitem__(
                "reduction", "arbitrary"
            ),
            "optimizer-budget": lambda value: value["optimizer"].__setitem__(
                "maximum_objective_evaluations_per_start", 1
            ),
            "parameter-bound": lambda value: value["parameters"][0].__setitem__(
                "upper", 1.0e19
            ),
            "runner-authority": lambda value: value["artifacts"][
                "parameterized_runner"
            ].__setitem__("sha256", "0" * 64),
            "interpretation": lambda value: value.__setitem__(
                "required_interpretation", ["weakened"]
            ),
            "authorization": lambda value: value["scope"].__setitem__(
                "production_ready", True
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                candidate = copy.deepcopy(document)
                mutate(candidate)
                with self.assertRaises(screen.ScreenError):
                    screen.validate_contract(candidate)

    def test_split_is_deterministic_and_partitions_are_isolated(self) -> None:
        training_ids = splitter.EXPECTED_CASES["protocol-training"]
        holdout_ids = splitter.EXPECTED_CASES["protocol-holdout"]
        document = {
            "schema": "synthetic-reference-for-split-test",
            "cases": [
                {
                    "id": case_id,
                    "partition": partition,
                    "partition_only_value": "only-" + case_id,
                }
                for partition, case_ids in splitter.EXPECTED_CASES.items()
                for case_id in case_ids
            ],
        }
        raw = splitter.canonical_bytes(document)
        with tempfile.TemporaryDirectory(prefix="rfam-r2-split-test-") as directory:
            source = Path(directory) / "source.json"
            source.write_bytes(raw)
            with mock.patch.object(splitter, "EXPECTED_SOURCE_SHA256", sha256(raw)), mock.patch.object(
                splitter, "validate_reference", return_value=([], document["cases"])
            ):
                first = splitter.split(source)
                second = splitter.split(source)

        self.assertEqual(first, second)
        training_raw = first["protocol-training"]
        holdout_raw = first["protocol-holdout"]
        training = json.loads(training_raw.decode("ascii"))
        holdout = json.loads(holdout_raw.decode("ascii"))
        self.assertEqual([case["id"] for case in training["cases"]], training_ids)
        self.assertEqual([case["id"] for case in holdout["cases"]], holdout_ids)
        for case_id in holdout_ids:
            self.assertNotIn(case_id.encode("ascii"), training_raw)
        for case_id in training_ids:
            self.assertNotIn(case_id.encode("ascii"), holdout_raw)
        self.assertEqual(training["partition_artifact"], "protocol-training")
        self.assertEqual(holdout["partition_artifact"], "protocol-holdout")

    def test_fit_parameter_header_is_stripped_only_when_exact(self) -> None:
        parameters = screen.DEFAULT_PARAMETERS
        header = "fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n"
        header += "".join(
            "fit_" + name + "=" + format(value, ".17g") + "\n"
            for name, value in zip(screen.EXPECTED_PARAMETER_NAMES, parameters)
        )
        trajectory = "schema=rfam.r2-0d-rfam-trajectory/v1\nR2_0D_RFAM_RUNNER_OK\n"
        self.assertEqual(
            screen.strip_fit_parameter_header(header + trajectory, parameters, "case"),
            trajectory,
        )

        rejected = {
            "missing-schema": "".join(header.splitlines(keepends=True)[1:]) + trajectory,
            "altered-value": (header + trajectory).replace("fit_R1_A_SI=10000000000", "fit_R1_A_SI=10000000001", 1),
            "truncated": "fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n",
            "crlf": (header + trajectory).replace("\n", "\r\n"),
        }
        for label, stdout in rejected.items():
            with self.subTest(label=label), self.assertRaises(screen.ScreenError):
                screen.strip_fit_parameter_header(stdout, parameters, "case")

    def test_holdout_cannot_be_installed_or_accessed_during_selection(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="ascii"))
        training_cases = [
            {"id": case_id}
            for case_id in contract["case_partition"]["training_case_ids"]
        ]
        holdout_cases = [
            {"id": case_id}
            for case_id in contract["case_partition"]["holdout_case_ids"]
        ]
        objective = screen.R2PhysicalObjective(
            Path("unused-runner"),
            [IdentitySpec() for _ in screen.EXPECTED_PARAMETER_NAMES],
            training_cases,
            [0.0, 1.0],
            contract,
        )
        self.assertNotIn("holdout", objective.case_ids)
        self.assertTrue(all(case["id"] not in objective.case_by_id for case in holdout_cases))
        with self.assertRaisesRegex(screen.ScreenError, "holdout-accessed-during-selection"):
            objective.evaluate([0.5] * 6, "holdout")
        self.assertEqual(objective.holdout_access_attempts_during_selection, 1)
        with self.assertRaisesRegex(screen.ScreenError, "holdout-installed-during-selection"):
            objective.install_holdout(holdout_cases, [0.0, 1.0])

        objective.selection_active = False
        objective.install_holdout(holdout_cases, [0.0, 1.0])
        self.assertEqual(objective.case_ids["holdout"], tuple(case["id"] for case in holdout_cases))
        with self.assertRaisesRegex(screen.ScreenError, "holdout-installation-state-invalid"):
            objective.install_holdout(holdout_cases, [0.0, 1.0])

    def test_repository_helper_drift_fails_preflight(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-helper-test-") as directory:
            root = Path(directory)
            helper = root / "tools" / "helper.py"
            helper.parent.mkdir(parents=True)
            helper.write_bytes(b"print('bound')\n")
            artifacts = {
                "repository_artifacts": [
                    {"path": "tools/helper.py", "sha256": sha256(helper.read_bytes())}
                ]
            }
            with mock.patch.object(screen, "EXPECTED_ARTIFACTS", artifacts):
                screen.verify_repository_artifacts(root)
                helper.write_bytes(b"print('drift')\n")
                with self.assertRaisesRegex(
                    screen.ScreenError, "repository-artifact-mismatch:tools/helper.py"
                ):
                    screen.verify_repository_artifacts(root)

    def test_duplicate_json_keys_are_rejected_by_both_tools(self) -> None:
        raw = b'{"scope":false,"scope":true}\n'
        with tempfile.TemporaryDirectory(prefix="rfam-r2-duplicate-test-") as directory:
            path = Path(directory) / "duplicate.json"
            path.write_bytes(raw)
            with self.assertRaisesRegex(screen.ScreenError, "duplicate-json-key:scope"):
                screen.load_json(path, "contract")
            with mock.patch.object(splitter, "EXPECTED_SOURCE_SHA256", sha256(raw)):
                with self.assertRaisesRegex(splitter.SplitError, "duplicate-json-key:scope"):
                    splitter.split(path)

    def test_failure_removes_stale_output_without_success_or_traceback(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-fail-close-test-") as directory:
            temporary = Path(directory)
            output = temporary / "stale-result.json"
            output.write_text('{"status":"stale-success"}\n', encoding="ascii")
            missing = temporary / "missing"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCREEN_TOOL),
                    "--contract",
                    str(CONTRACT),
                    "--training-reference",
                    str(missing / "training.json"),
                    "--holdout-reference",
                    str(missing / "holdout.json"),
                    "--runner",
                    str(missing / "runner"),
                    "--output",
                    str(output),
                ],
                check=False,
                text=True,
                capture_output=True,
                timeout=10,
            )
            combined = result.stdout + result.stderr
            self.assertEqual(result.returncode, 2, combined)
            self.assertIn("R2_UCSD_INTERNAL_FIT_SCREEN_FAIL", result.stderr)
            self.assertNotIn("R2_UCSD_INTERNAL_FIT_SCREEN_BLOCKED", combined)
            self.assertNotIn("Traceback", combined)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
