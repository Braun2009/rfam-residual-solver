#!/usr/bin/env python3
"""Negative-contract coverage for the blocked local Cantera R2 0D protocol."""

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
TOOL = ROOT / "tools" / "validate_r2_0d_calibration.py"
MANIFEST = ROOT / "tests" / "reference" / "r2_0d_calibration_contract_v1.json"
BASE_TEXT = MANIFEST.read_text(encoding="utf-8")
BASE = json.loads(BASE_TEXT)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class R2ZeroDCalibrationContractTests(unittest.TestCase):
    checks = 0

    def _run(self, document: dict[str, Any], root: Path = ROOT) -> subprocess.CompletedProcess[str]:
        contract = root / "tests" / "reference" / "r2_0d_calibration_contract_v1.json"
        contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
        result = subprocess.run(
            [sys.executable, str(TOOL), "--root", str(root), "--contract", str(contract)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.__class__.checks += 1
        return result

    def _run_raw(self, text: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-0d-raw-") as directory:
            contract = Path(directory) / "contract.json"
            contract.write_text(text, encoding="utf-8")
            result = subprocess.run(
                [sys.executable, str(TOOL), "--root", str(ROOT), "--contract", str(contract)],
                capture_output=True,
                text=True,
                check=False,
            )
        self.__class__.checks += 1
        return result

    def _mutated(self, mutation: Callable[[dict[str, Any]], None]) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(BASE)
        mutation(document)
        with tempfile.TemporaryDirectory(prefix="rfam-r2-0d-mutated-") as directory:
            contract = Path(directory) / "contract.json"
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            result = subprocess.run(
                [sys.executable, str(TOOL), "--root", str(ROOT), "--contract", str(contract)],
                capture_output=True,
                text=True,
                check=False,
            )
        self.__class__.checks += 1
        return result

    def _require_failure(self, result: subprocess.CompletedProcess[str], needle: str) -> None:
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("R2_0D_CALIBRATION_CONTRACT_FAIL", result.stderr)
        self.assertIn(needle, result.stderr)
        self.assertNotIn("R2_0D_CALIBRATION_CONTRACT_BLOCKED", result.stdout + result.stderr)
        self.assertNotIn("Traceback", result.stdout + result.stderr)

    def _copy_minimal_root(self, directory: Path) -> Path:
        relative_paths = [
            "tests/reference/r2_0d_calibration_contract_v1.json",
            "tests/reference/r2_0d_cantera_reference_v1.json",
            "tests/reference/cantera_gri30_v3_0.yaml",
            "tests/reference/r2_instantaneous_reference_v1.json",
            "tests/reference/r2_physics_candidate_v1.json",
            "scripts/generate_r2_0d_calibration_reference.py",
        ]
        for relative in relative_paths:
            source = ROOT / relative
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        return directory

    def test_base_contract_is_blocked(self) -> None:
        result = subprocess.run(
            [sys.executable, str(TOOL), "--root", str(ROOT), "--contract", str(MANIFEST)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_0D_CALIBRATION_CONTRACT_BLOCKED", result.stdout)

    def test_all_approval_flags_are_fail_closed(self) -> None:
        for name in BASE["approvals"]:
            with self.subTest(name=name):
                self._require_failure(
                    self._mutated(lambda document, name=name: document["approvals"].__setitem__(name, True)),
                    "approval-must-remain-false:" + name,
                )

    def test_contract_version_bool_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document.__setitem__("contract_version", True)),
            "invalid-int:root.contract_version",
        )

    def test_local_reference_cannot_be_selected_for_fit(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document["local_reference"]["mechanism"].__setitem__("status", "selected-for-fit")),
            "local_reference.mechanism.status",
        )

    def test_local_reference_fit_use_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document["local_reference"]["mechanism"].__setitem__("fit_use_allowed", True)),
            "local-reference-must-not-be-fit-eligible",
        )

    def test_network_or_fluent_cannot_be_claimed(self) -> None:
        for name in ("network_access", "fluent_access"):
            with self.subTest(name=name):
                self._require_failure(
                    self._mutated(lambda document, name=name: document["local_reference"]["runtime"].__setitem__(name, True)),
                    "local-reference-" + name.replace("_", "-") + "-must-remain-false",
                )

    def test_runtime_identity_cannot_be_rebound(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document["local_reference"]["runtime"].__setitem__("cantera_extension_sha256", "f" * 64)),
            "local_reference.runtime.cantera_extension_sha256",
        )

    def test_detailed_mechanism_and_experiment_gates_cannot_be_advanced(self) -> None:
        for key, status in (("detailed_mechanism", "selected"), ("experimental_training_dataset", "selected"), ("experimental_holdout_dataset", "selected"), ("metrics_and_tolerances", "frozen")):
            with self.subTest(key=key):
                self._require_failure(
                    self._mutated(lambda document, key=key, status=status: document["fit_and_experiment_gates"][key].__setitem__("status", status)),
                    "fit_and_experiment_gates." + key + ".status",
                )

    def test_protocol_partition_overlap_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document["protocol"].__setitem__("protocol_holdout_case_ids", ["train_phi0p8_diluted", "holdout_phi1p1", "holdout_phi1p4"])),
            "protocol.protocol_holdout_case_ids",
        )

    def test_protocol_observable_or_common_basis_weakening_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document["protocol"].__setitem__("common_basis_status", "established")),
            "protocol.common_basis_status",
        )
        self._require_failure(
            self._mutated(lambda document: document["protocol"]["declared_observables"].pop()),
            "protocol.declared_observables",
        )

    def test_artifact_path_and_hash_swap_fails_self_consistently(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first, second = document["artifacts"][0], document["artifacts"][1]
            for key in ("path", "sha256", "kind"):
                first[key], second[key] = second[key], first[key]

        self._require_failure(self._mutated(mutate), "artifact-path-contract-mismatch:local_cantera_mechanism")

    def test_reference_hash_rebind_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document["artifacts"][2].__setitem__("sha256", "a" * 64)),
            "artifact-sha256-contract-mismatch:reference_output",
        )

    def test_duplicate_root_scope_fails_before_last_value_wins(self) -> None:
        needle = '  "scope": "local-cantera-0d-protocol-only-no-fitting-or-physical-validation",'
        replacement = '  "scope": "physical-validation-approved",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(self._run_raw(BASE_TEXT.replace(needle, replacement, 1)), "duplicate-json-key:scope")

    def test_duplicate_nested_approval_fails_before_false_override(self) -> None:
        needle = '    "production_ready": false,'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(self._run_raw(BASE_TEXT.replace(needle, replacement, 1)), "duplicate-json-key:production_ready")

    def test_unknown_root_authorization_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda document: document.__setitem__("unreviewed_authorization", "allow-fit")),
            "unknown-keys:root:unreviewed_authorization",
        )

    def test_reference_regeneration_requires_the_frozen_wsl_entrypoint(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(TOOL),
                "--root",
                str(ROOT),
                "--contract",
                str(MANIFEST),
                "--verify-reference",
                "--cantera-python",
                "python3",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.__class__.checks += 1
        self._require_failure(result, "cantera-python-command-contract-mismatch")

    def test_self_consistently_rehashed_source_and_reference_still_fail(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-0d-rehash-") as directory:
            root = self._copy_minimal_root(Path(directory))
            generator = root / "scripts" / "generate_r2_0d_calibration_reference.py"
            generator.write_text(generator.read_text(encoding="utf-8") + "\n# self-consistent replacement\n", encoding="utf-8")
            output = root / "tests" / "reference" / "r2_0d_cantera_reference_v1.json"
            reference = json.loads(output.read_text(encoding="utf-8"))
            reference["generator_sha256"] = sha256(generator)
            output.write_text(json.dumps(reference, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
            contract = json.loads((root / "tests" / "reference" / "r2_0d_calibration_contract_v1.json").read_text(encoding="utf-8"))
            contract["artifacts"][1]["sha256"] = sha256(generator)
            contract["artifacts"][2]["sha256"] = sha256(output)
            result = self._run(contract, root)
        self._require_failure(result, "artifact-sha256-contract-mismatch:reference_generator")

    def test_self_consistent_reference_payload_rebind_fails(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-0d-output-") as directory:
            root = self._copy_minimal_root(Path(directory))
            output = root / "tests" / "reference" / "r2_0d_cantera_reference_v1.json"
            reference = json.loads(output.read_text(encoding="utf-8"))
            reference["cases"][0]["ignition_delay_s"] = 1.0e-6
            output.write_text(json.dumps(reference, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
            contract = json.loads((root / "tests" / "reference" / "r2_0d_calibration_contract_v1.json").read_text(encoding="utf-8"))
            contract["artifacts"][2]["sha256"] = sha256(output)
            result = self._run(contract, root)
        self._require_failure(result, "artifact-sha256-contract-mismatch:reference_output")


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(R2ZeroDCalibrationContractTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print("R2_0D_CALIBRATION_CONTRACT_TEST_PASS checks=" + str(R2ZeroDCalibrationContractTests.checks))
    raise SystemExit(0 if result.wasSuccessful() else 1)
