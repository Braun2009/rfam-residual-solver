#!/usr/bin/env python3
"""Adversarial coverage for the offline phase-caloric feasibility schema."""

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
TOOL = ROOT / "tools" / "validate_phase_caloric_feasibility_contract.py"
CONTRACT = ROOT / "tests" / "reference" / "phase_caloric_feasibility_contract_v1.json"
TEXT = CONTRACT.read_text(encoding="utf-8")
DOCUMENT = json.loads(TEXT)


class PhaseCaloricFeasibilityContractTests(unittest.TestCase):
    checks = 0

    def invoke(self, contract: Path) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(TOOL), "--root", str(ROOT), "--contract", str(contract)],
            check=False, text=True, capture_output=True,
        )
        self.__class__.checks += 1
        return result

    def require_failure(self, result: subprocess.CompletedProcess[str], needle: str) -> None:
        combined = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, combined)
        self.assertIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_FAIL", result.stderr)
        self.assertIn(needle, result.stderr)
        self.assertNotIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_BLOCKED", combined)
        self.assertNotIn("Traceback", combined)

    def mutated(self, mutation: Callable[[dict[str, Any]], None]) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(DOCUMENT)
        mutation(document)
        with tempfile.TemporaryDirectory(prefix="rfam-phase-caloric-contract-") as directory:
            path = Path(directory) / "contract.json"
            path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            return self.invoke(path)

    def raw(self, text: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-phase-caloric-contract-raw-") as directory:
            path = Path(directory) / "contract.json"
            path.write_text(text, encoding="utf-8")
            return self.invoke(path)

    def test_canonical_contract_is_explicitly_blocked(self) -> None:
        result = self.invoke(CONTRACT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_BLOCKED", result.stdout)
        self.assertIn("runtime_implementation=0", result.stdout)
        self.assertIn("two_phase_chemistry_enabled=0", result.stdout)
        self.assertIn("production_ready=0", result.stdout)

    def test_all_approval_flags_fail_closed(self) -> None:
        for name in DOCUMENT["approvals"]:
            with self.subTest(name=name):
                self.require_failure(
                    self.mutated(lambda d, name=name: d["approvals"].__setitem__(name, True)),
                    "contract-mismatch:approval-must-remain-false:" + name,
                )

    def test_duplicate_and_nonfinite_json_fail_before_authority_use(self) -> None:
        scope = '  "scope": "phase-local-vapor-or-supercritical-caloric-feasibility-no-runtime-implementation-no-two-phase-chemistry-success",'
        self.assertEqual(TEXT.count(scope), 1)
        self.require_failure(self.raw(TEXT.replace(scope, '  "scope": "runtime",\n' + scope, 1)),
                             "duplicate-json-key:scope")
        tolerance = '      "pressure_relative_residual_max": 1e-8,'
        self.assertEqual(TEXT.count(tolerance), 1)
        self.require_failure(self.raw(TEXT.replace(tolerance, tolerance.replace("1e-8", "NaN"), 1)),
                             "nonfinite-json-number:NaN")

    def test_phase_result_cannot_be_promoted_to_equilibrium_state(self) -> None:
        self.require_failure(
            self.mutated(lambda d: d["phase_caloric_schema"].__setitem__(
                "conversion_to_equilibrium_state", "allowed")),
            "contract-sha256-mismatch",
        )

    def test_seed_root_and_branch_weakening_fail_closed(self) -> None:
        mutations = [
            lambda d: d["immutable_request_binding"]["continuation"].__setitem__(
                "seed_role", "branch-authority"),
            lambda d: d["immutable_request_binding"].__setitem__(
                "automatic_vapor_to_supercritical", "allowed"),
            lambda d: d["immutable_request_binding"].__setitem__(
                "multi_root_fallback", "allowed"),
            lambda d: d["phase_local_validity"].__setitem__(
                "near_critical_ambiguity", "accept"),
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.require_failure(self.mutated(mutation), "contract-sha256-mismatch")

    def test_parent_token_and_restart_weakening_fail_closed(self) -> None:
        mutations = [
            lambda d: d["immutable_request_binding"]["continuation"].__setitem__(
                "parent_generation_mismatch", "accept"),
            lambda d: d["immutable_request_binding"]["continuation"].__setitem__(
                "root_continuity_required", False),
            lambda d: d["checkpoint_restart"].__setitem__(
                "token_reuse_across_rollback_restart_or_final_reflash", "accept"),
            lambda d: d["checkpoint_restart"].__setitem__(
                "generic_equilibrium_state_serialization", "allowed"),
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.require_failure(self.mutated(mutation), "contract-sha256-mismatch")

    def test_constant_pressure_and_publication_bypasses_fail_closed(self) -> None:
        mutations = [
            lambda d: d["chemistry_substep"].__setitem__(
                "constant_pressure_chemistry", "allowed"),
            lambda d: d["chemistry_substep"].__setitem__(
                "phase_density_relative_displacement_max", 0.01),
            lambda d: d["conservative_writeback_and_publication"].__setitem__(
                "only_publisher", "phase-caloric-result"),
            lambda d: d["conservative_writeback_and_publication"].__setitem__(
                "three_phase_or_unsupported_topology", "coerce-two-phase"),
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.require_failure(self.mutated(mutation), "contract-sha256-mismatch")

    def test_negative_matrix_and_source_scope_are_immutable(self) -> None:
        self.require_failure(
            self.mutated(lambda d: d["negative_cases_required"].pop()),
            "contract-sha256-mismatch",
        )
        self.require_failure(
            self.mutated(lambda d: d["evidence_boundary"]["source_files_allowed"].append(
                "src/thermo/EquilibriumThermoModel.H")),
            "contract-sha256-mismatch",
        )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        PhaseCaloricFeasibilityContractTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print("PHASE_CALORIC_FEASIBILITY_CONTRACT_TEST_PASS checks="
              + str(PhaseCaloricFeasibilityContractTests.checks))
    raise SystemExit(0 if result.wasSuccessful() else 1)
