#!/usr/bin/env python3
"""Adversarial checks for the blocked phase-local caloric v3 schema."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "validate_phase_caloric_feasibility_contract_v3.py"
CONTRACT = ROOT / "tests" / "reference" / "phase_caloric_feasibility_contract_v3.json"
RAW = CONTRACT.read_bytes()
DOCUMENT = json.loads(RAW.decode("utf-8"))

spec = importlib.util.spec_from_file_location("phase_caloric_v3_validator", TOOL)
assert spec is not None and spec.loader is not None
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


class PhaseCaloricFeasibilityContractV3Tests(unittest.TestCase):
    checks = 0

    def invoke(self, path: Path) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(TOOL), "--root", str(ROOT), "--contract", str(path)],
            check=False, text=True, capture_output=True,
        )
        self.__class__.checks += 1
        return result

    def require_failure(self, result: subprocess.CompletedProcess[str], needle: str) -> None:
        combined = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, combined)
        self.assertIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_V3_FAIL", result.stderr)
        self.assertIn(needle, result.stderr)
        self.assertNotIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_V3_BLOCKED", combined)
        self.assertNotIn("Traceback", combined)

    def raw_failure(self, raw: bytes, needle: str) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-phase-caloric-v3-raw-") as directory:
            path = Path(directory) / "contract.json"
            path.write_bytes(raw)
            self.require_failure(self.invoke(path), needle)

    def semantic_failure(self, mutation: Callable[[dict[str, Any]], None], needle: str) -> None:
        document = copy.deepcopy(DOCUMENT)
        mutation(document)
        raw = (json.dumps(document, indent=2) + "\n").encode("utf-8")
        previous = validator.EXPECTED_CONTRACT_SHA256
        validator.EXPECTED_CONTRACT_SHA256 = hashlib.sha256(raw).hexdigest()
        try:
            with self.assertRaises(validator.ContractError) as caught:
                validator._verify(document, raw)
        finally:
            validator.EXPECTED_CONTRACT_SHA256 = previous
        self.assertIn(needle, str(caught.exception))
        self.__class__.checks += 1

    def test_canonical_contract_remains_explicitly_blocked(self) -> None:
        result = self.invoke(CONTRACT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_V3_BLOCKED", result.stdout)
        self.assertIn("runtime_implementation=0", result.stdout)
        self.assertIn("two_phase_chemistry_enabled=0", result.stdout)
        self.assertIn("production_ready=0", result.stdout)

    def test_digest_and_json_syntax_fail_closed(self) -> None:
        mutated = copy.deepcopy(DOCUMENT)
        mutated["scope"] = "runtime"
        self.raw_failure((json.dumps(mutated, indent=2) + "\n").encode("utf-8"),
                         "contract-mismatch:scope")
        duplicate = RAW.replace(
            b'  "scope": "phase-local-vapor-or-supercritical-caloric-feasibility-no-runtime-implementation-no-two-phase-chemistry-success",\n',
            b'  "scope": "runtime",\n  "scope": "phase-local-vapor-or-supercritical-caloric-feasibility-no-runtime-implementation-no-two-phase-chemistry-success",\n',
            1,
        )
        self.raw_failure(duplicate, "duplicate-json-key:scope")
        nonfinite = RAW.replace(
            b'"pressure_relative": 1e-8', b'"pressure_relative": NaN', 1)
        self.raw_failure(nonfinite, "nonfinite-json-number:NaN")

    def test_every_approval_remains_false(self) -> None:
        for name in DOCUMENT["approvals"]:
            with self.subTest(name=name):
                self.semantic_failure(
                    lambda d, name=name: d["approvals"].__setitem__(name, True),
                    "contract-mismatch:approval:" + name,
                )

    def test_parent_generation_epoch_and_token_are_closed_schemas(self) -> None:
        cases = [
            (lambda d: d["parent_binding_schema"]["required_parent_fields"].remove(
                "parent_phase_payload_fingerprint"), "contract-mismatch:parent-fields"),
            (lambda d: d["parent_binding_schema"]["generation_identity"].__setitem__(
                "type", "string"), "contract-mismatch:generation-type"),
            (lambda d: d["parent_binding_schema"]["generation_identity"].__setitem__(
                "restore_policy", "restore-live-token"), "contract-mismatch:generation-restore"),
            (lambda d: d["parent_binding_schema"]["runtime_epoch"].__setitem__(
                "type", "optional"), "contract-mismatch:epoch-type"),
            (lambda d: d["root_token_schema"]["token_binds"].remove(
                "runtime_epoch"), "contract-mismatch:token-binds"),
            (lambda d: d["root_token_schema"].__setitem__(
                "rollback_policy", "restore-token-after-rollback"),
                "contract-mismatch:token-rollback"),
            (lambda d: d["root_token_schema"]["restart_hint"].__setitem__(
                "never_accepted_as_live_token", False), "contract-mismatch:hint-live-token"),
            (lambda d: d["root_token_schema"]["restart_hint"]["required_revalidation"].pop(),
                "contract-mismatch:hint-revalidation"),
        ]
        for mutation, needle in cases:
            with self.subTest(needle=needle):
                self.semantic_failure(mutation, needle)

    def test_success_replay_residual_and_diagnostics_envelope_is_closed(self) -> None:
        cases = [
            (lambda d: d["phase_caloric_schema"]["required_success_result_fields"].remove(
                "parent_q_fingerprint"), "contract-mismatch:success-result-fields"),
            (lambda d: d["phase_caloric_schema"]["required_success_result_fields"].remove(
                "parent_phase_payload_fingerprint"), "contract-mismatch:success-result-fields"),
            (lambda d: d["first_recovery_and_residuals"]["replay_required_fields"].remove(
                "compressibility_z"), "contract-mismatch:replay-fields"),
            (lambda d: d["first_recovery_and_residuals"]["replay_required_fields"].remove(
                "specific_enthalpy_j_kg"), "contract-mismatch:replay-fields"),
            (lambda d: d["first_recovery_and_residuals"]["residual_schema"].__setitem__(
                "units", "unspecified"), "contract-mismatch:residual-units"),
            (lambda d: d["first_recovery_and_residuals"]["residual_schema"][
                "required_residuals"].pop(), "contract-mismatch:required-residuals"),
            (lambda d: d["first_recovery_and_residuals"]["residual_schema"][
                "acceptance_max"].__setitem__("eos_root", 1.0),
                "contract-mismatch:residual-acceptance-max"),
            (lambda d: d["phase_caloric_schema"]["failure_result_policy"].__setitem__(
                "result_kind", "partial-state"), "contract-mismatch:failure-kind"),
            (lambda d: d["phase_caloric_schema"]["failure_result_policy"][
                "forbidden_on_failure"].remove("accepted_root_token"),
                "contract-mismatch:failure-forbidden-fields"),
        ]
        for mutation, needle in cases:
            with self.subTest(needle=needle):
                self.semantic_failure(mutation, needle)

    def test_phase_and_writeback_policies_are_semantically_checked(self) -> None:
        cases = [
            (lambda d: d["phase_local_validity"].__setitem__(
                "metastable_vapor_statement", "global-equilibrium"),
                "contract-mismatch:metastable-statement"),
            (lambda d: d["chemistry_substep"].__setitem__(
                "constant_pressure_chemistry", "allowed"), "contract-mismatch:constant-pressure"),
            (lambda d: d["chemistry_substep"].__setitem__(
                "specific_internal_energy_relative_displacement_max", 1.0),
                "contract-mismatch:energy-displacement"),
            (lambda d: d["chemistry_substep"].__setitem__(
                "mass_fraction_l1_displacement_max", 1.0),
                "contract-mismatch:composition-displacement"),
            (lambda d: d["phase_local_validity"].__setitem__(
                "composition_sum_absolute_tolerance", 1.0),
                "contract-mismatch:composition-tolerance"),
            (lambda d: d["conservative_writeback_and_publication"].__setitem__(
                "parent_vapor_volume_fraction_authority", "recompute-every-rhs"),
                "contract-mismatch:alpha-authority"),
            (lambda d: d["conservative_writeback_and_publication"][
                "final_reflash_requirements"].pop(), "contract-mismatch:final-reflash"),
            (lambda d: d["generic_solver_exclusion"].__setitem__(
                "phase_caloric_result_generic_solver_admission", "allowed"),
                "contract-mismatch:solver-exclusion:phase_caloric_result_generic_solver_admission"),
            (lambda d: d.__setitem__("stop_conditions", ["allow-everything"]),
                "contract-mismatch:stop-conditions"),
        ]
        for mutation, needle in cases:
            with self.subTest(needle=needle):
                self.semantic_failure(mutation, needle)

    def test_scope_and_negative_matrix_are_semantically_checked(self) -> None:
        self.semantic_failure(
            lambda d: d["negative_cases_required"].pop(),
            "contract-mismatch:negative-cases",
        )
        self.semantic_failure(
            lambda d: d["evidence_boundary"]["source_files_allowed"].append(
                "src/thermo/EquilibriumThermoModel.H"),
            "contract-mismatch:allowed-files",
        )
        self.semantic_failure(
            lambda d: d["phase_caloric_schema"].__setitem__(
                "conversion_to_equilibrium_state", "allowed"),
            "contract-mismatch:result-conversion",
        )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        PhaseCaloricFeasibilityContractV3Tests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print("PHASE_CALORIC_FEASIBILITY_CONTRACT_V3_TEST_PASS checks="
              + str(PhaseCaloricFeasibilityContractV3Tests.checks))
    raise SystemExit(0 if result.wasSuccessful() else 1)
