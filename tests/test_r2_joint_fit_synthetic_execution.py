#!/usr/bin/env python3
"""Focused and adversarial tests for synthetic R2 optimizer execution."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple


TEST_FILE = Path(__file__).resolve()
ROOT = TEST_FILE.parents[1]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

import run_r2_joint_fit_synthetic as runner  # noqa: E402
import validate_r2_joint_fit_synthetic_execution as package_gate  # noqa: E402


TOOL = TOOLS / "run_r2_joint_fit_synthetic.py"
GATE = TOOLS / "validate_r2_joint_fit_synthetic_execution.py"
EXECUTION_CONTRACT = (
    ROOT / "tests" / "reference" / "r2_joint_fit_synthetic_execution_contract_v1.json"
)
EXPECTED_RESULT = ROOT / "tests" / "reference" / "r2_joint_fit_synthetic_result_v1.json"
DOCUMENTATION = ROOT / "docs" / "r2_joint_fit_synthetic_execution.md"

EXPECTED_EXECUTION_HASH = "4c3db0f2711a00d93d8bdff4b08a72df41107db1a930cd072cc10fa7d4b0e268"
EXPECTED_RESULT_HASH = "582a00acc9779f0de4d1be2ecbe9efb91b5720fd464cb84e8d7c743d5aaf70ca"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def all_floats_finite(value: Any) -> bool:
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, list):
        return all(all_floats_finite(item) for item in value)
    if isinstance(value, dict):
        return all(all_floats_finite(item) for item in value.values())
    return True


class SyntheticR2ExecutionTests(unittest.TestCase):
    checks = 0

    def _invoke(
        self,
        root: Path,
        output: Path,
        env: Optional[Dict[str, str]] = None,
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(TOOL), "--root", str(root), "--output", str(output)],
            check=False,
            text=True,
            capture_output=True,
            env=env,
        )
        self.__class__.checks += 1
        return result

    def _invoke_gate(self, root: Path) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(GATE), "--root", str(root)],
            check=False,
            text=True,
            capture_output=True,
        )
        self.__class__.checks += 1
        return result

    def _minimal_runner_root(self, directory: Path) -> Path:
        paths = [
            EXECUTION_CONTRACT,
            *[
                ROOT / authority["path"]
                for authority in runner.FROZEN_V3_AUTHORITIES.values()
            ],
        ]
        for source in paths:
            target = directory / source.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        return directory

    def _minimal_gate_root(self, directory: Path) -> Path:
        self._minimal_runner_root(directory)
        for source in (TOOL, TEST_FILE, EXPECTED_RESULT, DOCUMENTATION):
            target = directory / source.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        return directory

    def _require_failure(
        self,
        result: subprocess.CompletedProcess[str],
        marker: str,
        success_marker: str,
    ) -> None:
        combined = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, combined)
        self.assertIn(marker, result.stderr)
        self.assertNotIn(success_marker, combined)
        self.assertNotIn("Traceback", combined)

    def test_bound_execution_contract_identity(self) -> None:
        self.assertEqual(sha256(EXECUTION_CONTRACT), EXPECTED_EXECUTION_HASH)
        self.assertEqual(EXECUTION_CONTRACT.stat().st_size, 6028)
        for name, authority in runner.FROZEN_V3_AUTHORITIES.items():
            path = ROOT / authority["path"]
            with self.subTest(name=name):
                self.assertEqual(sha256(path), authority["sha256"])
                self.assertEqual(path.stat().st_size, authority["byte_count"])

    def test_execution_contract_semantics_reject_synchronized_rehash(self) -> None:
        document = json.loads(EXECUTION_CONTRACT.read_text(encoding="utf-8"))
        cases = [
            (
                "scope",
                lambda value: value.__setitem__("scope", "physical-fit"),
                "contract-mismatch:execution.scope",
            ),
            (
                "missing-approval",
                lambda value: value["approvals"].pop("production_ready"),
                "contract-mismatch:execution.approvals",
            ),
            (
                "unknown-approval",
                lambda value: value["approvals"].__setitem__("new_approval", False),
                "contract-mismatch:execution.approvals",
            ),
            (
                "enabled-approval",
                lambda value: value["approvals"].__setitem__("production_ready", True),
                "contract-mismatch:execution.approvals",
            ),
            (
                "approval-int-zero",
                lambda value: value.__setitem__(
                    "approvals", {key: 0 for key in runner.EXECUTION_APPROVALS}
                ),
                "contract-mismatch:execution.approvals",
            ),
            (
                "adapter",
                lambda value: value["synthetic_runner"].__setitem__(
                    "adapter_id", "physical-runner"
                ),
                "contract-mismatch:execution.synthetic_runner.adapter_id",
            ),
            (
                "slope",
                lambda value: value["synthetic_runner"]["training_response_entries"][0]
                .__setitem__("raw_slope_per_normalized_coordinate", 1.0),
                "contract-mismatch:execution.response-entries",
            ),
            (
                "integer-slope-alias",
                lambda value: value["synthetic_runner"]["training_response_entries"][0]
                .__setitem__("raw_slope_per_normalized_coordinate", 100000000),
                "contract-mismatch:execution.response-entries",
            ),
            (
                "holdout-response",
                lambda value: value["synthetic_runner"].__setitem__(
                    "holdout_parameter_response", "affine"
                ),
                "contract-mismatch:execution.synthetic_runner.holdout_parameter_response",
            ),
            (
                "jacobi-tolerance",
                lambda value: value["execution_protocol"].__setitem__(
                    "jacobi_relative_offdiagonal_tolerance", 1e-6
                ),
                "contract-mismatch:execution.execution_protocol",
            ),
            (
                "unknown-protocol-field",
                lambda value: value["execution_protocol"].__setitem__("parallel", True),
                "contract-mismatch:execution.execution_protocol",
            ),
            (
                "protocol-bool-count-aliases",
                lambda value: value["execution_protocol"].update(
                    {
                        "holdout_runner_calls_during_selection": False,
                        "holdout_runner_calls_after_selection": True,
                    }
                ),
                "contract-mismatch:execution.execution_protocol",
            ),
            (
                "protocol-float-sweep-alias",
                lambda value: value["execution_protocol"].__setitem__(
                    "jacobi_maximum_sweeps", 100.0
                ),
                "contract-mismatch:execution.execution_protocol",
            ),
            (
                "unknown-claim",
                lambda value: value["claims"].__setitem__("production_authorized", False),
                "contract-mismatch:execution.claims",
            ),
            (
                "enabled-claim",
                lambda value: value["claims"].__setitem__(
                    "physical_optimizer_executed", True
                ),
                "contract-mismatch:execution.claims",
            ),
            (
                "claim-int-aliases",
                lambda value: value.__setitem__(
                    "claims",
                    {
                        key: (1 if expected else 0)
                        for key, expected in runner.EXECUTION_CLAIMS.items()
                    },
                ),
                "contract-mismatch:execution.claims",
            ),
            (
                "authority-float-byte-count-alias",
                lambda value: value["frozen_v3_authorities"]["objective_contract"]
                .__setitem__("byte_count", 15830.0),
                "contract-mismatch:execution.frozen_v3_authorities",
            ),
            (
                "stop-condition",
                lambda value: value["stop_conditions"].pop(),
                "contract-mismatch:execution.stop_conditions",
            ),
        ]
        for label, mutate, needle in cases:
            with self.subTest(label=label):
                candidate = copy.deepcopy(document)
                mutate(candidate)
                with self.assertRaises(runner.SyntheticExecutionError) as caught:
                    runner._verify_execution_contract(candidate)
                self.assertIn(needle, str(caught.exception))

    def test_result_schema_rejects_authorization_and_type_aliases(self) -> None:
        document = json.loads(EXPECTED_RESULT.read_text(encoding="utf-8"))
        cases = [
            (
                "missing-negative-claim",
                lambda value: value["claims"].pop("production_ready"),
            ),
            (
                "unknown-false-claim",
                lambda value: value["claims"].__setitem__("new_approval", False),
            ),
            (
                "numeric-claim-alias",
                lambda value: value["claims"].__setitem__("production_ready", 0),
            ),
            (
                "unknown-root-authorization",
                lambda value: value.__setitem__("unknown_root", {"production_ready": True}),
            ),
            (
                "missing-scope-limit",
                lambda value: value["scope_limits"].pop(),
            ),
            (
                "authority-float-byte-count-alias",
                lambda value: value["authorities"]["frozen_v3"]["objective_contract"]
                .__setitem__("byte_count", 15830.0),
            ),
            (
                "runner-field-order",
                lambda value: value["runner"]["required_case_request_fields"].reverse(),
            ),
            (
                "runner-unknown-field",
                lambda value: value["runner"].__setitem__("physical_runner", True),
            ),
            (
                "optimizer-unknown-field",
                lambda value: value["optimizer"].__setitem__("production_ready", False),
            ),
            (
                "per-start-unknown-field",
                lambda value: value["optimizer"]["per_start"][0]
                .__setitem__("accepted", True),
            ),
            (
                "uncertainty-output-drift",
                lambda value: value["identifiability"]["uncertainty_outputs"]
                .__setitem__("physical_evidence", 0),
            ),
            (
                "holdout-bool-int-alias",
                lambda value: value["holdout"].__setitem__("used_during_selection", 0),
            ),
        ]
        for label, mutate in cases:
            with self.subTest(label=label):
                candidate = copy.deepcopy(document)
                mutate(candidate)
                with self.assertRaises(package_gate.PackageError):
                    package_gate._verify_result(candidate)

    def test_canonical_result_is_byte_deterministic(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-synthetic-run-") as directory:
            first = Path(directory) / "first.json"
            second = Path(directory) / "second.json"
            first_result = self._invoke(ROOT, first)
            environment = os.environ.copy()
            environment["PYTHONHASHSEED"] = "937"
            second_result = self._invoke(ROOT, second, environment)
            self.assertEqual(first_result.returncode, 0, first_result.stdout + first_result.stderr)
            self.assertEqual(second_result.returncode, 0, second_result.stdout + second_result.stderr)
            self.assertIn(runner.BLOCKED_MARKER, first_result.stdout)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(first.read_bytes(), EXPECTED_RESULT.read_bytes())
            self.assertEqual(sha256(first), EXPECTED_RESULT_HASH)

    def test_result_scope_accounting_selection_and_identifiability(self) -> None:
        result = json.loads(EXPECTED_RESULT.read_text(encoding="utf-8"))
        self.assertEqual(result["status"], "blocked-synthetic-optimizer-execution-only")
        self.assertEqual(result["watermark"], "SYNTHETIC_NONPHYSICAL_TEST_ONLY")
        self.assertEqual(result["optimizer"]["start_count"], 17)
        self.assertEqual(result["optimizer"]["selected_start_index"], 12)
        self.assertEqual(result["optimizer"]["selected_start_kind"], "halton-12")
        self.assertAlmostEqual(result["optimizer"]["training_objective"], 0.170813797398293, 15)
        self.assertAlmostEqual(result["holdout"]["objective"], 0.317243682724575, 15)
        self.assertEqual(result["runner"]["selection_training_calls"], 3284)
        self.assertEqual(result["runner"]["case_requests_per_training_call"], 2)
        self.assertEqual(result["runner"]["case_requests_per_holdout_call"], 2)
        self.assertEqual(
            result["runner"]["initial_state_authority"],
            "ordered-initial-mass-fraction-vector-no-reconstruction",
        )
        self.assertEqual(len(result["runner"]["required_case_request_fields"]), 8)
        self.assertEqual(
            result["runner"]["selection_training_calls"],
            sum(start["runner_evaluations"] for start in result["optimizer"]["per_start"]),
        )
        self.assertEqual(result["runner"]["holdout_calls_during_selection"], 0)
        self.assertEqual(result["runner"]["holdout_calls_after_selection"], 1)
        self.assertEqual(result["holdout"]["evaluation_count"], 1)
        self.assertFalse(result["holdout"]["used_during_selection"])
        self.assertEqual(result["runner"]["identifiability_training_calls"], 13)
        self.assertEqual(result["runner"]["total_calls"], 3298)
        self.assertEqual(result["identifiability"]["row_count"], 58)
        self.assertEqual(result["identifiability"]["rank"], 6)
        self.assertEqual(result["identifiability"]["status"], "synthetic-full-rank")
        self.assertTrue(result["identifiability"]["rank_condition_correlation_pass"])
        self.assertLess(result["identifiability"]["condition_number"], 2.0)
        self.assertEqual(
            result["identifiability"]["maximum_absolute_offdiagonal_correlation"], 0.0
        )
        self.assertEqual(
            result["identifiability"]["uncertainty_outputs"]["status"],
            "not-publishable-synthetic-covariance-diagnostic-only",
        )
        self.assertTrue(all_floats_finite(result))
        for key, value in result["claims"].items():
            if key == "synthetic_optimizer_executed":
                self.assertTrue(value)
            else:
                self.assertFalse(value, key)

    def test_selected_point_matches_affine_optimum_within_minimum_step(self) -> None:
        result = json.loads(EXPECTED_RESULT.read_text(encoding="utf-8"))
        anchor = result["runner"]["anchor_normalized_coordinates"]
        expected_delta = [-0.1, 0.2, -0.1, -0.25, 0.5, 0.5]
        optimum = [value + delta for value, delta in zip(anchor, expected_delta)]
        selected = result["optimizer"]["selected_normalized_coordinates"]
        for actual, expected in zip(selected, optimum):
            self.assertLessEqual(abs(actual - expected), 0.0009765625)

    def test_halton_definition_and_start_order(self) -> None:
        self.assertEqual(runner.radical_inverse(1, 2), 0.5)
        self.assertEqual(runner.radical_inverse(2, 2), 0.25)
        self.assertEqual(runner.radical_inverse(3, 2), 0.75)
        self.assertAlmostEqual(runner.radical_inverse(2, 3), 2.0 / 3.0, 16)
        contract = json.loads(
            (ROOT / runner.FROZEN_V3_AUTHORITIES["objective_contract"]["path"])
            .read_text(encoding="utf-8")
        )
        specs = runner._parameter_specs(contract)
        starts = runner._start_points(contract, len(specs))
        self.assertEqual(starts[0][0], "midpoint")
        self.assertEqual(starts[1][0], "halton-1")
        self.assertEqual(starts[-1][0], "halton-16")
        self.assertEqual(starts[1][1][0], 0.5)

    def test_coordinate_candidates_share_entry_and_tie_break_lexicographically(self) -> None:
        calls = []
        spec = (runner.ParameterSpec("x", "linear", 0.0, 1.0, 0.0, 1.0),)

        def evaluate(point: Tuple[float, ...]) -> runner.ObjectiveEvaluation:
            calls.append(point)
            objective = 0.0 if point[0] in (0.25, 0.75) else 1.0
            return runner.ObjectiveEvaluation(True, objective, (math.sqrt(objective),), {})

        result = runner.run_pattern_search(
            0,
            "test",
            (0.5,),
            spec,
            runner.SearchConfig(0.25, 0.5, 0.25, 3, 1),
            evaluate,
        )
        self.assertEqual(calls, [(0.5,), (0.75,), (0.25,)])
        self.assertEqual(result.selected_u, (0.25,))
        self.assertEqual(result.evaluation_count, 3)
        self.assertEqual(result.termination, "evaluation-limit-before-coordinate")

    def test_per_start_cache_and_budget_preflight(self) -> None:
        spec = (runner.ParameterSpec("x", "linear", 0.0, 1.0, 0.0, 1.0),)
        cache_calls = []

        def cache_evaluate(point: Tuple[float, ...]) -> runner.ObjectiveEvaluation:
            cache_calls.append(point)
            objective = (point[0] - 0.75) ** 2
            return runner.ObjectiveEvaluation(True, objective, (math.sqrt(objective),), {})

        cached = runner.run_pattern_search(
            0,
            "cache",
            (0.5,),
            spec,
            runner.SearchConfig(0.25, 0.5, 0.25, 6, 1),
            cache_evaluate,
        )
        self.assertEqual(cached.selected_u, (0.75,))
        self.assertEqual(cached.cache_hits, 1)
        self.assertEqual(cache_calls.count((0.5,)), 1)

        second_calls = []
        runner.run_pattern_search(
            1,
            "cache-new-start",
            (0.5,),
            spec,
            runner.SearchConfig(0.25, 0.5, 0.25, 6, 1),
            lambda point: (
                second_calls.append(point)
                or runner.ObjectiveEvaluation(True, (point[0] - 0.75) ** 2, (0.0,), {})
            ),
        )
        self.assertEqual(second_calls.count((0.5,)), 1)

        budget_calls = []
        budget = runner.run_pattern_search(
            2,
            "budget",
            (0.5,),
            spec,
            runner.SearchConfig(0.25, 0.5, 0.25, 2, 1),
            lambda point: (
                budget_calls.append(point)
                or runner.ObjectiveEvaluation(True, point[0], (point[0],), {})
            ),
        )
        self.assertEqual(budget_calls, [(0.5,)])
        self.assertEqual(budget.evaluation_count, 1)
        self.assertEqual(budget.termination, "evaluation-limit-before-coordinate")

    def test_invalid_start_and_candidate_fail_closed(self) -> None:
        spec = (runner.ParameterSpec("x", "linear", 0.0, 1.0, 0.0, 1.0),)
        invalid_start = runner.run_pattern_search(
            0,
            "invalid-start",
            (0.5,),
            spec,
            runner.SearchConfig(0.25, 0.5, 0.25, 10, 1),
            lambda point: runner.ObjectiveEvaluation.invalid("synthetic-failure"),
        )
        self.assertFalse(invalid_start.valid)
        self.assertEqual(invalid_start.evaluation_count, 1)
        self.assertEqual(invalid_start.termination, "invalid-start")

        def invalid_positive(point: Tuple[float, ...]) -> runner.ObjectiveEvaluation:
            if point[0] > 0.5:
                return runner.ObjectiveEvaluation.invalid("synthetic-failure")
            return runner.ObjectiveEvaluation(True, 0.5 - point[0], (0.0,), {})

        result = runner.run_pattern_search(
            1,
            "invalid-candidate",
            (0.5,),
            spec,
            runner.SearchConfig(0.25, 0.5, 0.25, 6, 1),
            invalid_positive,
        )
        self.assertTrue(result.valid)
        self.assertEqual(result.selected_u, (0.5,))

    def test_jacobi_and_rank_deficient_identifiability(self) -> None:
        eigenvalues, eigenvectors, _ = runner.symmetric_jacobi(
            [[2.0, 1.0], [1.0, 2.0]], 1e-15, 100
        )
        self.assertAlmostEqual(eigenvalues[0], 3.0, 14)
        self.assertAlmostEqual(eigenvalues[1], 1.0, 14)
        self.assertAlmostEqual(
            sum(value * value for value in (eigenvectors[0][0], eigenvectors[1][0])),
            1.0,
            14,
        )

        def flat_residuals(point: Tuple[float, ...]) -> runner.ObjectiveEvaluation:
            del point
            residuals = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            return runner.ObjectiveEvaluation(True, 1.0, residuals, {})

        result = runner.analyze_identifiability(
            (0.5, 0.5), 2, 1e-4, 1e-8, 1e8, 0.98, 1e-15, 100,
            flat_residuals,
        )
        self.assertEqual(result["status"], "nonidentifiable")
        self.assertEqual(result["rank"], 0)
        self.assertIsNone(result["covariance_matrix"])
        self.assertFalse(result["rank_condition_correlation_pass"])

    def test_execution_or_v3_authority_tamper_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-synthetic-tamper-") as directory:
            root = self._minimal_runner_root(Path(directory))
            execution = root / EXECUTION_CONTRACT.relative_to(ROOT)
            execution.write_bytes(execution.read_bytes() + b" ")
            output = root / "result.json"
            result = self._invoke(root, output)
            self._require_failure(result, runner.FAIL_MARKER, runner.BLOCKED_MARKER)
            self.assertFalse(output.exists())

        with tempfile.TemporaryDirectory(prefix="rfam-r2-v3-tamper-") as directory:
            root = self._minimal_runner_root(Path(directory))
            fixture = root / runner.FROZEN_V3_AUTHORITIES["synthetic_fixture"]["path"]
            fixture.write_bytes(fixture.read_bytes() + b" ")
            output = root / "result.json"
            result = self._invoke(root, output)
            self._require_failure(result, runner.FAIL_MARKER, runner.BLOCKED_MARKER)
            self.assertFalse(output.exists())

    def test_self_consistent_execution_rebinding_and_duplicate_key_fail(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-execution-rebind-") as directory:
            root = self._minimal_runner_root(Path(directory))
            execution = root / EXECUTION_CONTRACT.relative_to(ROOT)
            document = json.loads(execution.read_text(encoding="utf-8"))
            source = root / runner.FROZEN_V3_AUTHORITIES["objective_contract"]["path"]
            substitute = root / "tests" / "reference" / "substitute-objective.json"
            shutil.copy2(source, substitute)
            document["frozen_v3_authorities"]["objective_contract"]["path"] = str(
                substitute.relative_to(root)
            ).replace("\\", "/")
            document["frozen_v3_authorities"]["objective_contract"]["sha256"] = sha256(substitute)
            document["frozen_v3_authorities"]["objective_contract"]["byte_count"] = substitute.stat().st_size
            execution.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            output = root / "result.json"
            result = self._invoke(root, output)
            self._require_failure(result, runner.FAIL_MARKER, runner.BLOCKED_MARKER)
            self.assertFalse(output.exists())

        with tempfile.TemporaryDirectory(prefix="rfam-r2-execution-duplicate-") as directory:
            root = self._minimal_runner_root(Path(directory))
            execution = root / EXECUTION_CONTRACT.relative_to(ROOT)
            text = execution.read_text(encoding="utf-8")
            needle = '  "status": "blocked-synthetic-optimizer-execution-only",'
            execution.write_text(
                text.replace(needle, '  "status": "production-ready",\n' + needle, 1),
                encoding="utf-8",
            )
            output = root / "result.json"
            result = self._invoke(root, output)
            self._require_failure(result, runner.FAIL_MARKER, runner.BLOCKED_MARKER)
            self.assertFalse(output.exists())

    def test_gate_accepts_only_frozen_execution_package(self) -> None:
        result = self._invoke_gate(ROOT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_JOINT_FIT_SYNTHETIC_PACKAGE_BLOCKED", result.stdout)
        self.assertIn("physical_fit_executed=0 production_ready=0", result.stdout)

    def test_gate_rejects_self_consistently_rehashed_substitute_runner_before_execution(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-substitute-runner-") as directory:
            root = self._minimal_gate_root(Path(directory))
            substitute = root / TOOL.relative_to(ROOT)
            sentinel = root / "substitute-executed.txt"
            substitute.write_text(
                "from pathlib import Path\n"
                + "Path(" + repr(str(sentinel)) + ").write_text('executed', encoding='utf-8')\n",
                encoding="utf-8",
            )
            execution = root / EXECUTION_CONTRACT.relative_to(ROOT)
            document = json.loads(execution.read_text(encoding="utf-8"))
            rebound = document["frozen_v3_authorities"]["objective_validator"]
            rebound["path"] = str(substitute.relative_to(root)).replace("\\", "/")
            rebound["sha256"] = sha256(substitute)
            rebound["byte_count"] = substitute.stat().st_size
            execution.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            result = self._invoke_gate(root)
            self._require_failure(
                result,
                "R2_JOINT_FIT_SYNTHETIC_PACKAGE_FAIL",
                "R2_JOINT_FIT_SYNTHETIC_PACKAGE_BLOCKED",
            )
            self.assertFalse(sentinel.exists())

    def test_gate_rejects_self_consistently_replaced_result_test_or_docs(self) -> None:
        for source in (EXPECTED_RESULT, TEST_FILE, DOCUMENTATION):
            with self.subTest(source=source.name):
                with tempfile.TemporaryDirectory(prefix="rfam-r2-package-tamper-") as directory:
                    root = self._minimal_gate_root(Path(directory))
                    target = root / source.relative_to(ROOT)
                    target.write_bytes(target.read_bytes() + b" ")
                    result = self._invoke_gate(root)
                    self._require_failure(
                        result,
                        "R2_JOINT_FIT_SYNTHETIC_PACKAGE_FAIL",
                        "R2_JOINT_FIT_SYNTHETIC_PACKAGE_BLOCKED",
                    )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(SyntheticR2ExecutionTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(
            "R2_JOINT_FIT_SYNTHETIC_EXECUTION_TEST_PASS checks="
            + str(SyntheticR2ExecutionTests.checks)
        )
    raise SystemExit(0 if result.wasSuccessful() else 1)
