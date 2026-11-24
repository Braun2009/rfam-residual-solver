#!/usr/bin/env python3
"""Focused regressions for the frozen UCSD internal fit-screen result."""

from __future__ import annotations

import copy
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
VALIDATOR = TOOLS / "validate_r2_ucsd_internal_fit_screen_result.py"
TMP = Path("F:/tmp") if Path("F:/tmp").is_dir() else Path("/mnt/f/tmp")
RESULT = TMP / "a408-r2-ucsd-internal-fit-screen-result-v4.json"
TRAINING = TMP / "a408-r2-ucsd-internal-training-v1.json"
HOLDOUT = TMP / "a408-r2-ucsd-internal-holdout-v1.json"
RUNNER = TMP / "rfam_r2_0d_fit_runner"
BASE_RAW = RESULT.read_bytes()
BASE = json.loads(BASE_RAW.decode("ascii"))

sys.path.insert(0, str(TOOLS))
import validate_r2_ucsd_internal_fit_screen_result as validator  # noqa: E402


class R2UcsdInternalFitScreenResultTests(unittest.TestCase):
    def invoke(
        self,
        result: Path = RESULT,
        root: Path = ROOT,
        training: Path = TRAINING,
        holdout: Path = HOLDOUT,
        runner: Path = RUNNER,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VALIDATOR),
                "--root",
                str(root),
                "--result",
                str(result),
                "--training-reference",
                str(training),
                "--holdout-reference",
                str(holdout),
                "--runner",
                str(runner),
            ],
            check=False,
            text=True,
            capture_output=True,
            timeout=10,
        )

    def require_failure(
        self,
        result: subprocess.CompletedProcess[str],
        reason: str,
    ) -> None:
        combined = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, combined)
        self.assertIn("R2_UCSD_INTERNAL_FIT_SCREEN_RESULT_FAIL", result.stderr)
        self.assertIn(reason, result.stderr)
        self.assertNotIn("R2_UCSD_INTERNAL_FIT_SCREEN_RESULT_BLOCKED", combined)
        self.assertNotIn("Traceback", combined)

    def mutated(
        self,
        mutation: Callable[[dict[str, Any]], None],
        canonical: bool = True,
    ) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(BASE)
        mutation(document)
        with tempfile.TemporaryDirectory(prefix="rfam-r2-result-mutation-") as directory:
            path = Path(directory) / "result.json"
            if canonical:
                path.write_bytes(validator.canonical_bytes(document))
            else:
                path.write_text(json.dumps(document, indent=2) + "\n", encoding="ascii")
            return self.invoke(result=path)

    def test_canonical_package_is_blocked(self) -> None:
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(
            "R2_UCSD_INTERNAL_FIT_SCREEN_RESULT_BLOCKED checks=106 "
            "physical_validation=0 production_ready=0",
            result.stdout,
        )
        self.assertNotIn("FAIL", result.stdout + result.stderr)

    def test_recursive_unknown_missing_and_type_alias_fields_fail(self) -> None:
        mutations = {
            "unknown-nested": (
                lambda value: value["baseline"].__setitem__("extra", 0),
                "key-schema-mismatch:baseline",
            ),
            "missing-nested": (
                lambda value: value["starts"][0].pop("trace_sha256"),
                "key-schema-mismatch:starts[0]",
            ),
            "bool-for-int": (
                lambda value: value["runner_accounting"].__setitem__(
                    "training_objective_calls", True
                ),
                "type-mismatch:runner_accounting.training_objective_calls",
            ),
            "int-for-float": (
                lambda value: value["selection"].__setitem__(
                    "training_objective", 172707
                ),
                "type-mismatch:selection.training_objective",
            ),
            "wrong-list-shape": (
                lambda value: value["selection"][
                    "selected_normalized_coordinates"
                ].pop(),
                "selected-coordinate-shape-mismatch",
            ),
        }
        for label, (mutation, reason) in mutations.items():
            with self.subTest(label=label):
                self.require_failure(self.mutated(mutation), reason)

    def test_every_claim_has_exact_boolean_value(self) -> None:
        for key, expected in validator.EXPECTED_CLAIMS.items():
            with self.subTest(key=key):
                self.require_failure(
                    self.mutated(
                        lambda value, key=key, expected=expected: value["claims"].__setitem__(
                            key, not expected
                        )
                    ),
                    "value-mismatch:claims." + key,
                )

    def test_optimizer_and_accounting_relations_fail_closed(self) -> None:
        mutations = {
            "wrong-selected-start": (
                lambda value: value["selection"].__setitem__("selected_start_index", 1),
                "selected-start-kind-mismatch",
            ),
            "evaluation-total": (
                lambda value: value["starts"][4].__setitem__("evaluation_count", 111),
                "training-objective-accounting-mismatch",
            ),
            "training-calls": (
                lambda value: value["runner_accounting"].__setitem__(
                    "training_objective_calls", 648
                ),
                "training-objective-accounting-mismatch",
            ),
            "training-attempts": (
                lambda value: value["runner_accounting"].__setitem__(
                    "training_case_process_attempts", 1942
                ),
                "training-attempt-accounting-mismatch",
            ),
            "holdout-access": (
                lambda value: value["runner_accounting"].__setitem__(
                    "holdout_access_attempts_during_selection", 1
                ),
                "holdout-accounting-mismatch",
            ),
            "valid-count": (
                lambda value: value["selection"].__setitem__("valid_start_count", 4),
                "valid-start-count-mismatch",
            ),
            "selected-start-parameters": (
                lambda value: (
                    value["starts"][0].__setitem__(
                        "selected_normalized_coordinates",
                        list(value["starts"][1]["selected_normalized_coordinates"]),
                    ),
                    value["starts"][0].__setitem__(
                        "selected_parameter_vector",
                        list(value["starts"][1]["selected_parameter_vector"]),
                    ),
                ),
                "selected-start-parameter-mismatch",
            ),
        }
        for label, (mutation, reason) in mutations.items():
            with self.subTest(label=label):
                self.require_failure(self.mutated(mutation), reason)

    def test_source_authority_mutations_fail_semantically(self) -> None:
        mutations = {
            "contract": lambda value: value["source"].__setitem__(
                "contract_sha256", "0" * 64
            ),
            "training": lambda value: value["source"]["training_reference"].__setitem__(
                "sha256", "0" * 64
            ),
            "runner": lambda value: value["source"]["parameterized_runner"].__setitem__(
                "sha256", "0" * 64
            ),
            "helper": lambda value: value["source"]["repository_artifacts"][4].__setitem__(
                "sha256", "0" * 64
            ),
        }
        for label, mutation in mutations.items():
            with self.subTest(label=label):
                self.require_failure(self.mutated(mutation), "value-mismatch:source")

    def test_duplicate_keys_fail_without_success_marker(self) -> None:
        needle = '"production_ready":false'
        self.assertEqual(BASE_RAW.decode("ascii").count(needle), 1)
        duplicate = BASE_RAW.decode("ascii").replace(
            needle,
            '"production_ready":true,"production_ready":false',
            1,
        )
        with tempfile.TemporaryDirectory(prefix="rfam-r2-result-duplicate-") as directory:
            path = Path(directory) / "result.json"
            path.write_text(duplicate, encoding="ascii")
            self.require_failure(self.invoke(result=path), "duplicate-json-key:production_ready")

    def test_noncanonical_and_hash_drift_are_distinguished(self) -> None:
        self.require_failure(
            self.mutated(lambda value: None, canonical=False),
            "result-not-canonical",
        )
        self.require_failure(
            self.mutated(
                lambda value: value["runner_accounting"][
                    "invalid_candidate_reasons"
                ].__setitem__(
                    "ignition-not-detected:train_phi0p8_diluted", 75
                )
            ),
            "result-sha256-mismatch",
        )
        with tempfile.TemporaryDirectory(prefix="rfam-r2-result-copy-") as directory:
            copy_path = Path(directory) / RESULT.name
            copy_path.write_bytes(BASE_RAW)
            self.require_failure(self.invoke(result=copy_path), "result-path-mismatch")

    def test_selected_coordinates_bind_the_parameter_vector(self) -> None:
        self.require_failure(
            self.mutated(
                lambda value: value["selection"]["selected_parameter_vector"].__setitem__(
                    0, value["selection"]["selected_parameter_vector"][0] * 1.01
                )
            ),
            "selected-parameter-transform-mismatch:0",
        )

    def test_repository_and_external_drift_are_rejected_pre_execution(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-result-root-") as directory:
            root = Path(directory)
            for relative in validator.EXPECTED_REPOSITORY:
                source = ROOT / relative
                target = root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
            drift = root / "tools" / "run_r2_joint_fit_synthetic.py"
            drift.write_bytes(drift.read_bytes() + b"\n")
            self.require_failure(
                self.invoke(root=root),
                "repository-artifact-mismatch:tools/run_r2_joint_fit_synthetic.py",
            )

        with tempfile.TemporaryDirectory(prefix="rfam-r2-result-runner-") as directory:
            temporary = Path(directory)
            sentinel = temporary / "executed.txt"
            substitute = temporary / "runner.py"
            substitute.write_text(
                "from pathlib import Path\n"
                + "Path(" + repr(str(sentinel)) + ").write_text('executed')\n",
                encoding="ascii",
            )
            self.require_failure(
                self.invoke(runner=substitute),
                "external-artifact-mismatch:parameterized_runner",
            )
            self.assertFalse(sentinel.exists())


if __name__ == "__main__":
    unittest.main()
