#!/usr/bin/env python3
"""Adversarial coverage for the synthetic-only joint R1/R2 objective gate."""

from __future__ import annotations

import copy
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "validate_r2_joint_fit_objective.py"
CONTRACT = ROOT / "tests" / "reference" / "r2_joint_fit_objective_contract_v1.json"
FIXTURE = ROOT / "tests" / "reference" / "r2_joint_fit_synthetic_fixture_v1.json"
CONTRACT_TEXT = CONTRACT.read_text(encoding="utf-8")
FIXTURE_TEXT = FIXTURE.read_text(encoding="utf-8")
CONTRACT_DOCUMENT = json.loads(CONTRACT_TEXT)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class JointR2FitObjectiveContractTests(unittest.TestCase):
    checks = 0

    def _invoke(self, root: Path, contract: Path) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(TOOL), "--root", str(root), "--contract", str(contract)],
            check=False,
            text=True,
            capture_output=True,
        )
        self.__class__.checks += 1
        return result

    def _run_document(self, document: dict[str, Any]) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-contract-") as directory:
            contract = Path(directory) / "contract.json"
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            return self._invoke(ROOT, contract)

    def _run_raw_contract(self, text: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-contract-raw-") as directory:
            contract = Path(directory) / "contract.json"
            contract.write_text(text, encoding="utf-8")
            return self._invoke(ROOT, contract)

    def _mutated_contract(self, mutation: Callable[[dict[str, Any]], None]) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(CONTRACT_DOCUMENT)
        mutation(document)
        return self._run_document(document)

    def _minimal_root(self, directory: Path) -> Path:
        for source in (CONTRACT, FIXTURE):
            relative = source.relative_to(ROOT)
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        return directory

    def _require_failure(self, result: subprocess.CompletedProcess[str], needle: str) -> None:
        combined = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, combined)
        self.assertIn("R2_JOINT_FIT_OBJECTIVE_FAIL", result.stderr)
        self.assertIn(needle, result.stderr)
        self.assertNotIn("R2_JOINT_FIT_OBJECTIVE_SYNTHETIC_ONLY_BLOCKED", combined)
        self.assertNotIn("Traceback", combined)

    def test_canonical_contract_recomputes_separate_split_objectives_and_remains_blocked(self) -> None:
        result = self._invoke(ROOT, CONTRACT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_JOINT_FIT_OBJECTIVE_SYNTHETIC_ONLY_BLOCKED", result.stdout)
        self.assertIn("training_cases=2 holdout_cases=2", result.stdout)
        self.assertIn("training_objective=0.1828606650975268", result.stdout)
        self.assertIn("holdout_objective=0.31724368272457459", result.stdout)
        self.assertIn(
            "fit_executed=0 physical_validation_approved=0 production_ready=0",
            result.stdout,
        )
        self.assertNotIn("R2_JOINT_FIT_OBJECTIVE_FAIL", result.stdout + result.stderr)

    def test_every_approval_flag_fails_closed(self) -> None:
        for name in CONTRACT_DOCUMENT["approvals"]:
            with self.subTest(name=name):
                result = self._mutated_contract(
                    lambda document, name=name: document["approvals"].__setitem__(name, True)
                )
                self._require_failure(result, "approval-must-remain-false:" + name)

    def test_duplicate_root_and_nested_keys_fail_before_last_value_wins(self) -> None:
        root_needle = '  "scope": "mechanism-agnostic-joint-r1-r2-objective-schema-no-fit-no-physical-validation",'
        self.assertEqual(CONTRACT_TEXT.count(root_needle), 1)
        self._require_failure(
            self._run_raw_contract(
                CONTRACT_TEXT.replace(
                    root_needle,
                    '  "scope": "production-fit",\n' + root_needle,
                    1,
                )
            ),
            "duplicate-json-key:scope",
        )
        approval_needle = '    "production_ready": false'
        self.assertEqual(CONTRACT_TEXT.count(approval_needle), 1)
        self._require_failure(
            self._run_raw_contract(
                CONTRACT_TEXT.replace(
                    approval_needle,
                    '    "production_ready": true,\n' + approval_needle,
                    1,
                )
            ),
            "duplicate-json-key:production_ready",
        )

    def test_duplicate_fixture_key_is_rejected_before_identity_comparison(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-fixture-duplicate-") as directory:
            root = self._minimal_root(Path(directory))
            fixture = root / FIXTURE.relative_to(ROOT)
            needle = '  "watermark": "SYNTHETIC_NONPHYSICAL_TEST_ONLY",'
            self.assertEqual(FIXTURE_TEXT.count(needle), 1)
            fixture.write_text(
                FIXTURE_TEXT.replace(
                    needle,
                    '  "watermark": "PHYSICAL_REFERENCE",\n' + needle,
                    1,
                ),
                encoding="utf-8",
            )
            result = self._invoke(root, root / CONTRACT.relative_to(ROOT))
        self._require_failure(result, "duplicate-json-key:watermark")

    def test_nonfinite_json_in_contract_or_fixture_is_structured_failure(self) -> None:
        needle = '      "ignition_delay_sigma_log": 0.1,'
        self.assertEqual(CONTRACT_TEXT.count(needle), 1)
        self._require_failure(
            self._run_raw_contract(CONTRACT_TEXT.replace(needle, needle.replace("0.1", "NaN"), 1)),
            "nonfinite-json-number:NaN",
        )
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-fixture-nan-") as directory:
            root = self._minimal_root(Path(directory))
            fixture = root / FIXTURE.relative_to(ROOT)
            needle = '        "ignition_delay_s": 0.0003,'
            self.assertEqual(FIXTURE_TEXT.count(needle), 1)
            fixture.write_text(FIXTURE_TEXT.replace(needle, needle.replace("0.0003", "NaN"), 1),
                               encoding="utf-8")
            result = self._invoke(root, root / CONTRACT.relative_to(ROOT))
        self._require_failure(result, "nonfinite-json-number:NaN")

    def test_unknown_authorization_field_fails_closed(self) -> None:
        self._require_failure(
            self._mutated_contract(
                lambda document: document.__setitem__("unreviewed_fit_authorization", True)
            ),
            "unknown-keys:root:unreviewed_fit_authorization",
        )

    def test_contract_version_bool_fails(self) -> None:
        self._require_failure(
            self._mutated_contract(lambda document: document.__setitem__("contract_version", True)),
            "invalid-int:root.contract_version",
        )

    def test_scientific_contract_weakening_changes_frozen_authority(self) -> None:
        mutations: list[tuple[str, Callable[[dict[str, Any]], None]]] = [
            ("holdout-use", lambda d: d["partitions"].__setitem__(
                "holdout_objective_calls_during_parameter_selection", 1)),
            ("partition", lambda d: d["partitions"]["holdout_case_ids"].__setitem__(
                0, "SYNTHETIC_TRAIN_001")),
            ("observable", lambda d: d["observables_and_basis"]["observable_order"].pop()),
            ("enthalpy", lambda d: d["observables_and_basis"]["energy_basis"].__setitem__(
                "basis_offset_fitting_allowed", True)),
            ("weight", lambda d: d["metrics_weights_tolerances"]["component_weights"].__setitem__(
                "qdot", 0.0)),
            ("physical-tolerance", lambda d: d["metrics_weights_tolerances"]
             ["physical_acceptance_tolerances"].__setitem__("status", "approved")),
            ("parameter-bound", lambda d: d["parameter_bounds"]["ordered_fit_parameters"][5]
             .__setitem__("upper", 10.0)),
            ("concentration-unit", lambda d: d["parameter_bounds"].__setitem__(
                "concentration_unit", "kmol/m3")),
            ("activation-energy-unit", lambda d: d["parameter_bounds"]
             ["ordered_fit_parameters"][1].__setitem__("units", "J/kmol")),
            ("reversible-r2-order", lambda d: d["parameter_bounds"]["fixed_parameters"]
             .__setitem__("R2_order_O2", 0.25)),
            ("optimizer", lambda d: d["optimizer_protocol"].__setitem__(
                "selection_split", "training-and-holdout")),
            ("optimizer-halton-origin", lambda d: d["optimizer_protocol"]["start_sequence"]
             .__setitem__("halton_radical_inverse_index_origin", 0)),
            ("optimizer-cache-scope", lambda d: d["optimizer_protocol"]["search"]
             .__setitem__("cache_scope", "global-across-starts")),
            ("identifiability", lambda d: d["identifiability_uncertainty"].__setitem__(
                "status", "passed")),
            ("identifiability-scaling", lambda d: d["identifiability_uncertainty"]
             ["scaled_residual_definition"].__setitem__("holdout_rows_excluded", False)),
            ("identifiability-covariance", lambda d: d["identifiability_uncertainty"]
             ["uncertainty_requirements"].__setitem__(
                 "covariance_definition", "unspecified")),
            ("rights", lambda d: d["provenance_rights"].__setitem__(
                "physical_fit_use_permission_bound", True)),
            ("domain", lambda d: d["domain"].__setitem__(
                "phase_scope", "gas-and-two-phase")),
            ("initial-state", lambda d: d["domain"]["initial_composition"]
             .__setitem__("diluent_species", "unspecified")),
        ]
        for label, mutation in mutations:
            with self.subTest(label=label):
                self._require_failure(
                    self._mutated_contract(mutation),
                    "contract-sha256-contract-mismatch",
                )

    def test_self_consistently_rehashed_fixture_substitution_fails(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-substitute-") as directory:
            root = self._minimal_root(Path(directory))
            fixture = root / FIXTURE.relative_to(ROOT)
            payload = json.loads(fixture.read_text(encoding="utf-8"))
            payload["cases"][0]["prediction"]["temperature_k"][1] = 1220.0
            fixture.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            contract = root / CONTRACT.relative_to(ROOT)
            document = json.loads(contract.read_text(encoding="utf-8"))
            document["artifact"]["sha256"] = sha256(fixture)
            document["artifact"]["byte_count"] = fixture.stat().st_size
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            result = self._invoke(root, contract)
        self._require_failure(result, "contract-sha256-contract-mismatch")

    def test_self_consistent_fixture_path_rebinding_fails(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-rebind-") as directory:
            root = self._minimal_root(Path(directory))
            source = root / FIXTURE.relative_to(ROOT)
            replacement = root / "tests" / "reference" / "substitute_fixture.json"
            shutil.copy2(source, replacement)
            contract = root / CONTRACT.relative_to(ROOT)
            document = json.loads(contract.read_text(encoding="utf-8"))
            document["artifact"]["path"] = "tests/reference/substitute_fixture.json"
            document["artifact"]["sha256"] = sha256(replacement)
            document["artifact"]["byte_count"] = replacement.stat().st_size
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            result = self._invoke(root, contract)
        self._require_failure(result, "contract-sha256-contract-mismatch")

    def test_fixture_tamper_with_canonical_contract_fails_without_objective_marker(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-joint-tamper-") as directory:
            root = self._minimal_root(Path(directory))
            fixture = root / FIXTURE.relative_to(ROOT)
            fixture.write_text(fixture.read_text(encoding="utf-8") + " ", encoding="utf-8")
            result = self._invoke(root, root / CONTRACT.relative_to(ROOT))
        self._require_failure(result, "artifact-byte-count-mismatch:synthetic_objective_fixture")


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        JointR2FitObjectiveContractTests
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(
            "R2_JOINT_FIT_OBJECTIVE_CONTRACT_TEST_PASS checks="
            + str(JointR2FitObjectiveContractTests.checks)
        )
    raise SystemExit(0 if result.wasSuccessful() else 1)
