#!/usr/bin/env python3
"""Negative and positive tests for the WP2 provenance gate."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "validate_r2_physics_candidate.py"
MANIFEST = ROOT / "tests" / "reference" / "r2_physics_candidate_v1.json"
BASE = json.loads(MANIFEST.read_text(encoding="utf-8"))
BASE_TEXT = MANIFEST.read_text(encoding="utf-8")


class R2PhysicsCandidateTests(unittest.TestCase):
    checks = 0

    def _run(self, document: dict[str, Any]) -> subprocess.CompletedProcess[str]:
        return self._run_raw(json.dumps(document, indent=2) + "\n")

    def _run_raw(self, payload: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-contract-") as temp:
            path = Path(temp) / "candidate.json"
            path.write_text(payload, encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOL),
                    "--root",
                    str(ROOT),
                    "--contract",
                    str(path),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
        self.__class__.checks += 1
        return result

    def _mutated(self, mutate: Callable[[dict[str, Any]], None]) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(BASE)
        mutate(document)
        return self._run(document)

    def _require_failure(self, result: subprocess.CompletedProcess[str], fragment: str) -> None:
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("R2_PHYSICS_CANDIDATE_FAIL", result.stdout)
        self.assertIn(fragment, result.stdout)
        self.assertNotIn("R2_PHYSICS_CANDIDATE_BLOCKED", result.stdout)
        self.assertNotIn("Traceback", result.stdout + result.stderr)

    def test_frozen_candidate_is_blocked(self) -> None:
        result = self._run(copy.deepcopy(BASE))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_PHYSICS_CANDIDATE_BLOCKED", result.stdout)

    def test_artifact_path_swap_fails_even_when_hashes_follow(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first, second = document["artifacts"][0], document["artifacts"][1]
            for key in ("path", "url", "sha256", "kind"):
                first[key], second[key] = second[key], first[key]

        self._require_failure(self._mutated(mutate), "artifact-path-contract-mismatch")

    def test_duplicate_artifact_path_key_fails_before_last_value_wins(self) -> None:
        needle = '      "path": "Vortex/tools/openfoam/sandiego20161214_mechCK.txt",'
        replacement = '      "path": "Vortex/tools/openfoam/unrelated.txt",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run_raw(BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:path",
        )

    def test_artifact_order_change_fails(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["artifacts"][0], document["artifacts"][1] = document["artifacts"][1], document["artifacts"][0]

        self._require_failure(self._mutated(mutate), "artifact-order-contract-mismatch")

    def test_artifact_url_swap_fails(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first, second = document["artifacts"][0], document["artifacts"][1]
            first["url"], second["url"] = second["url"], first["url"]

        self._require_failure(self._mutated(mutate), "artifact[ucsd_mechanism].url")

    def test_artifact_hash_tamper_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["artifacts"][0].__setitem__("sha256", "a" * 64)),
            "artifact-sha-contract-mismatch",
        )

    def test_upstream_path_and_hash_swap_fails(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first, second = document["upstream"][0], document["upstream"][1]
            first["path"], second["path"] = second["path"], first["path"]
            first["sha256"], second["sha256"] = second["sha256"], first["sha256"]

        self._require_failure(self._mutated(mutate), "upstream-path-contract-mismatch")

    def test_generator_version_change_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["conversion"]["generator"].__setitem__("version", "3.1.0")),
            "conversion.generator.version",
        )

    def test_generator_wheel_link_change_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["conversion"]["generator"].__setitem__("wheel_artifact_id", "converted_yaml")),
            "conversion.generator.wheel_artifact_id",
        )

    def test_runtime_artifact_path_change_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["conversion"]["generator"]["runtime_artifacts"][0].__setitem__("path", "C:/other/python.exe")),
            "runtime-artifact-path-contract-mismatch",
        )

    def test_normalized_yaml_fingerprint_change_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["conversion"]["output"].__setitem__("normalized_sha256", "b" * 64)),
            "conversion-normalized-sha-contract-mismatch",
        )

    def test_runtime_species_fingerprint_change_fails(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["runtime_probe"]["observed"].__setitem__("species_order_sha256", "c" * 64)),
            "runtime_probe.observed.species_order_sha256",
        )

    def test_runtime_warning_record_change_fails(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["runtime_probe"]["observed"]["warnings"]["records"][0]["species"] = "O2"

        self._require_failure(self._mutated(mutate), "runtime_probe.observed.warnings.records")

    def test_duplicate_nested_runtime_version_fails(self) -> None:
        needle = '    "python": {\n      "version": "3.8.10",'
        replacement = '    "python": {\n      "version": "0.0.0",\n      "version": "3.8.10",'
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run_raw(BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:version",
        )

    def test_license_gate_cannot_be_deleted(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            del document["mechanism_candidate"]["license"]["required_before_selection"]

        self._require_failure(self._mutated(mutate), "missing-keys:mechanism_candidate.license")

    def test_candidate_cannot_be_marked_selected(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["mechanism_candidate"].__setitem__("status", "selected")),
            "mechanism_candidate.status",
        )

    def test_training_dataset_cannot_be_silently_filled(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["experimental_evidence"]["training_dataset"].__setitem__("status", "selected")),
            "experimental_evidence.training_dataset",
        )

    def test_duplicate_nested_experimental_status_fails(self) -> None:
        needle = '    "training_dataset": {\n      "status": "unselected",'
        replacement = '    "training_dataset": {\n      "status": "selected",\n      "status": "unselected",'
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run_raw(BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:status",
        )

    def test_holdout_partition_is_required(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            del document["experimental_evidence"]["holdout_dataset"]

        self._require_failure(self._mutated(mutate), "missing-keys:experimental_evidence")

    def test_metrics_must_be_frozen_before_fit(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["validation_protocol"]["metrics"].__setitem__("status", "frozen")),
            "validation_protocol.metrics",
        )

    def test_approval_flag_is_not_an_authorization_path(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["approvals"].__setitem__("physical_validation_approved", True)),
            "approval-must-remain-false",
        )

    def test_duplicate_approval_key_fails_before_false_override(self) -> None:
        needle = '    "production_ready": false,'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run_raw(BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:production_ready",
        )

    def test_cross_scope_quality_owner_cannot_move_to_uds(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d["cross_scope"]["hem"].__setitem__("quality_uds_authority", True)),
            "cross_scope.hem",
        )

    def test_experiment_observable_object_fails_structured(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["experimental_evidence"]["training_dataset"]["observables"] = [{"name": "ignition_delay"}]

        self._require_failure(self._mutated(mutate), "experimental_evidence.training_dataset")

    def test_contract_version_bool_is_not_integer(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d.__setitem__("contract_version", True)),
            "invalid-int:root.contract_version",
        )

    def test_unknown_root_key_fails_without_success_marker(self) -> None:
        self._require_failure(
            self._mutated(lambda d: d.__setitem__("unreviewed_authorization", True)),
            "unknown-keys:root:unreviewed_authorization",
        )

    def test_duplicate_root_scope_fails_before_frozen_override(self) -> None:
        needle = '  "scope": "provenance-screening-only-no-fitting-or-physical-validation",'
        replacement = '  "scope": "physical-validation-approved",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run_raw(BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:scope",
        )

    def test_optional_runtime_and_regeneration_evidence(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(TOOL),
                "--root",
                str(ROOT),
                "--contract",
                str(MANIFEST),
                "--verify-runtime-probe",
                "--verify-regeneration",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_PHYSICS_CANDIDATE_BLOCKED", result.stdout)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(R2PhysicsCandidateTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"R2_PHYSICS_CANDIDATE_TEST_PASS checks={R2PhysicsCandidateTests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
