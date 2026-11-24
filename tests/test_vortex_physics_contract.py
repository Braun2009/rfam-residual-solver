#!/usr/bin/env python3
"""Negative-case matrix for the WP0 Vortex physics contract."""

from __future__ import annotations

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
VALIDATOR = ROOT / "tools" / "validate_vortex_physics_contract.py"
CONTRACT_RELATIVE = Path("tests/reference/vortex_physics_validation_contract_v1.json")
CONTRACT = ROOT / CONTRACT_RELATIVE


class PhysicsContractTests(unittest.TestCase):
    checks = 0

    def _check(self, condition: bool, message: str) -> None:
        type(self).checks += 1
        self.assertTrue(condition, message)

    def _run(self, root: Path, contract: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VALIDATOR),
                "--root",
                str(root),
                "--contract",
                str(contract),
            ],
            check=False,
            text=True,
            capture_output=True,
        )

    def _copy_contract_tree(self, destination: Path) -> Path:
        document = json.loads(CONTRACT.read_text(encoding="utf-8"))
        for artifact in document["artifacts"]:
            relative = Path(artifact["path"])
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, target)
        contract_target = destination / CONTRACT_RELATIVE
        contract_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(CONTRACT, contract_target)
        return contract_target

    @staticmethod
    def _read_contract(path: Path) -> dict[str, Any]:
        return json.loads(path.read_text(encoding="utf-8"))

    @staticmethod
    def _write_contract(path: Path, document: dict[str, Any]) -> None:
        path.write_text(
            json.dumps(document, indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )

    def _mutated_result(
        self,
        callback: Callable[[dict[str, Any], Path], None],
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            contract = self._copy_contract_tree(root)
            document = self._read_contract(contract)
            callback(document, root)
            self._write_contract(contract, document)
            return self._run(root, contract)

    def _require_failure(
        self,
        result: subprocess.CompletedProcess[str],
        expected: str,
    ) -> None:
        self._check(result.returncode == 2, result.stdout + result.stderr)
        self._check(expected in result.stderr, result.stderr)
        self._check("VORTEX_PHYSICS_CONTRACT_FAIL" in result.stderr, result.stderr)
        self._check("VORTEX_PHYSICS_CONTRACT_BLOCKED" not in result.stdout, result.stdout)
        self._check("Traceback" not in result.stderr, result.stderr)

    @staticmethod
    def _artifact(document: dict[str, Any], artifact_id: str) -> dict[str, Any]:
        return next(item for item in document["artifacts"] if item["id"] == artifact_id)

    def test_frozen_contract_is_structurally_valid_but_blocked(self) -> None:
        result = self._run(ROOT, CONTRACT)
        self._check(result.returncode == 0, result.stderr)
        self._check("VORTEX_PHYSICS_CONTRACT_BLOCKED" in result.stdout, result.stdout)
        self._check("production_ready=0" in result.stdout, result.stdout)
        self._check("VORTEX_PHYSICS_CONTRACT_FAIL" not in result.stdout, result.stdout)

    def test_missing_artifact_fails(self) -> None:
        def mutate(document: dict[str, Any], root: Path) -> None:
            (root / "scripts" / "cantera_reference.py").unlink()

        self._require_failure(
            self._mutated_result(mutate), "artifact-missing:cantera_generator"
        )

    def test_modified_artifact_fails(self) -> None:
        def mutate(document: dict[str, Any], root: Path) -> None:
            target = root / "tests/reference/cantera_gri30_propane_oxygen_v1.txt"
            target.write_bytes(target.read_bytes() + b"stale=1\n")

        self._require_failure(
            self._mutated_result(mutate), "artifact-sha256-mismatch:cantera_reference"
        )

    def test_artifact_ids_are_bound_to_canonical_paths(self) -> None:
        def mutate(document: dict[str, Any], root: Path) -> None:
            source = self._artifact(document, "two_step_source")
            generator = self._artifact(document, "cantera_generator")
            source["path"], generator["path"] = generator["path"], source["path"]
            source["sha256"], generator["sha256"] = (
                generator["sha256"],
                source["sha256"],
            )

        self._require_failure(
            self._mutated_result(mutate),
            "artifact-path-contract-mismatch:cantera_generator",
        )

    def test_contract_version_bool_fails(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document.__setitem__("contract_version", True)
            ),
            "invalid-int:root.contract_version",
        )

    def test_production_flag_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document.__setitem__("production_ready", True)
            ),
            "production_ready-must-be-false",
        )

    def test_path_traversal_fails_before_read(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document["artifacts"][0].__setitem__(
                    "path", "../outside"
                )
            ),
            "path-escapes-root",
        )

    def test_energy_basis_is_exact(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document["energy_basis"].__setitem__(
                    "flow_energy", "anything-nonempty"
                )
            ),
            "energy-basis-contract-mismatch",
        )

    def test_envelope_required_fields_cannot_be_removed(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document["target_envelope"].pop(
                    "required_before_fit"
                )
            ),
            "missing-keys:target_envelope:required_before_fit",
        )

    def test_envelope_anchor_must_be_positive(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document["target_envelope"]["known_anchors"].__setitem__(
                    "pressure_pa", [-1.0, 4000000.0]
                )
            ),
            "target-envelope-anchor-nonpositive:pressure_pa",
        )

    def test_envelope_anchor_must_be_unique_and_sorted(self) -> None:
        for values in ([4000000.0, 101325.0], [101325.0, 101325.0]):
            with self.subTest(values=values):
                self._require_failure(
                    self._mutated_result(
                        lambda document, root, values=values: document["target_envelope"][
                            "known_anchors"
                        ].__setitem__("pressure_pa", values)
                    ),
                    "target-envelope-anchor-order-or-duplicate:pressure_pa",
                )

    def test_r2_contract_is_closed_to_deletion_and_extension(self) -> None:
        mutations = [
            (lambda document, root: document["r2"].pop("detailed_mechanism"),
             "missing-keys:r2:detailed_mechanism"),
            (lambda document, root: document["r2"].pop("experimental_dataset"),
             "missing-keys:r2:experimental_dataset"),
            (lambda document, root: document["r2"].pop("required_observables"),
             "missing-keys:r2:required_observables"),
            (lambda document, root: document["r2"].__setitem__("authorization", True),
             "unknown-keys:r2:authorization"),
        ]
        for mutation, expected in mutations:
            with self.subTest(expected=expected):
                self._require_failure(self._mutated_result(mutation), expected)

    def test_kij_contract_fields_cannot_be_removed(self) -> None:
        mutations = [
            (lambda document, root: document["kij"]["explicit_zero"].pop("required_fields"),
             "missing-keys:kij.explicit_zero:required_fields"),
            (lambda document, root: document["kij"].pop("sensitivity_outputs"),
             "missing-keys:kij:sensitivity_outputs"),
            (lambda document, root: document["kij"].__setitem__(
                "missing_pair_policy", "allow-if-insensitive"),
             "kij-missing-policy-mismatch"),
        ]
        for mutation, expected in mutations:
            with self.subTest(expected=expected):
                self._require_failure(self._mutated_result(mutation), expected)

    def test_hem_dependencies_cannot_be_removed(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document["hem"].pop("box_dependencies")
            ),
            "missing-keys:hem:box_dependencies",
        )

    def test_soot_physics_requirements_cannot_be_removed(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document["condensed_carbon"].pop(
                    "physics_required"
                )
            ),
            "missing-keys:condensed_carbon:physics_required",
        )

    def test_comparison_contract_cannot_be_weakened(self) -> None:
        mutations = [
            (lambda document, root: document["comparison"].__setitem__(
                "quantitative_single_phase_only", False),
             "comparison-must-remain-single-phase-only"),
            (lambda document, root: document["comparison"].pop(
                "quantitative_single_phase_only"),
             "missing-keys:comparison:quantitative_single_phase_only"),
            (lambda document, root: document["comparison"]["excluded_quantities"].append({}),
             "invalid-string:comparison.excluded_quantities"),
        ]
        for mutation, expected in mutations:
            with self.subTest(expected=expected):
                self._require_failure(self._mutated_result(mutation), expected)

    def test_external_origin_rejects_zero_hash_and_unknown_snapshot(self) -> None:
        mutations = [
            (lambda document, root: document["external_origins"][0].__setitem__(
                "source_sha256", "0" * 64),
             "zero-sha256-forbidden:external_origins[0].source_sha256"),
            (lambda document, root: document["external_origins"][0].__setitem__(
                "snapshot_artifact", "missing-snapshot"),
             "external-origin-snapshot-artifact-missing:missing-snapshot"),
        ]
        for mutation, expected in mutations:
            with self.subTest(expected=expected):
                self._require_failure(self._mutated_result(mutation), expected)

    def test_snapshot_provenance_is_bound_even_with_matching_artifact_hash(self) -> None:
        def mutate(document: dict[str, Any], root: Path) -> None:
            artifact = self._artifact(document, "injection_snapshot")
            path = root / artifact["path"]
            snapshot = json.loads(path.read_text(encoding="utf-8"))
            snapshot["source_sha256"] = "1" * 64
            path.write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
            artifact["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()

        self._require_failure(
            self._mutated_result(mutate), "external-snapshot-source-sha256-mismatch"
        )

    def test_external_origin_is_bound_to_canonical_source_path(self) -> None:
        def mutate(document: dict[str, Any], root: Path) -> None:
            live_source = root / "unrelated-live-source.json"
            payload = {"unrelated": True}
            live_source.write_text(
                json.dumps(payload, indent=2) + "\n", encoding="utf-8"
            )
            source_path = str(live_source.resolve())
            source_sha = hashlib.sha256(live_source.read_bytes()).hexdigest()

            origin = document["external_origins"][0]
            origin["source_path"] = source_path
            origin["source_sha256"] = source_sha

            artifact = self._artifact(document, "injection_snapshot")
            snapshot_path = root / artifact["path"]
            snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
            snapshot["source_path"] = source_path
            snapshot["source_sha256"] = source_sha
            snapshot["payload"] = payload
            snapshot_path.write_text(
                json.dumps(snapshot, indent=2) + "\n", encoding="utf-8"
            )
            artifact["sha256"] = hashlib.sha256(
                snapshot_path.read_bytes()
            ).hexdigest()

        self._require_failure(
            self._mutated_result(mutate),
            "external-origin-source-path-contract-mismatch",
        )

    def test_unknown_root_key_fails(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document.__setitem__("unreviewed_override", True)
            ),
            "unknown-keys:root:unreviewed_override",
        )

    def test_ready_status_is_not_an_authorization_path(self) -> None:
        self._require_failure(
            self._mutated_result(
                lambda document, root: document.__setitem__("status", "ready")
            ),
            "status-must-remain-blocked",
        )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PhysicsContractTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"VORTEX_PHYSICS_CONTRACT_TEST_PASS checks={PhysicsContractTests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
