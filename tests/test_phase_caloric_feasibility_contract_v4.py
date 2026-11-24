#!/usr/bin/env python3
"""Focused adversarial tests for the blocked phase-caloric v4 design."""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "tests/reference/phase_caloric_feasibility_contract_v4.json"
ORACLE = ROOT / "tools/phase_caloric_lifecycle_oracle_v4.py"
VALIDATOR = ROOT / "tools/validate_phase_caloric_feasibility_contract_v4.py"
DOC = ROOT / "docs/phase_caloric_feasibility_contract_v4.md"
ARTIFACT_PATHS = (
    "tests/reference/phase_caloric_feasibility_contract_v4.json",
    "tools/phase_caloric_lifecycle_oracle_v4.py",
    "tools/validate_phase_caloric_feasibility_contract_v4.py",
    "tests/test_phase_caloric_feasibility_contract_v4.py",
    "docs/phase_caloric_feasibility_contract_v4.md",
)
EXPECTED_VALIDATOR_SHA256 = (
    "fbf0e07afe0ce9c3afad64d83d9c0f361efa9ea195c3dcccee018dd70e7af30c"
)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


if sha256_bytes(VALIDATOR.read_bytes()) != EXPECTED_VALIDATOR_SHA256:
    raise RuntimeError("phase-caloric-v4-validator-identity-mismatch")

spec = importlib.util.spec_from_file_location(
    "phase_caloric_v4_validator_verified", VALIDATOR
)
assert spec is not None and spec.loader is not None
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)

DOCUMENT, _, oracle = validator.validate(ROOT)
assert oracle is not None


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("ascii")).hexdigest()


def f64be(value: float) -> str:
    return struct.pack(">d", value).hex()


def parent_binding(
    owner: str,
    *,
    epoch: int = 7,
    generation: int = 11,
    q: str | None = None,
    phase_payload: str | None = None,
    phase: str = "VAPOR",
    species_count: int = 2,
) -> object:
    return oracle.ParentBinding(
        parent_owner_id=owner,
        committed_parent_q_fingerprint=q or digest(owner + "-q"),
        parent_phase_payload_fingerprint=(
            phase_payload or digest(owner + "-phase")
        ),
        parent_phase_identity=phase,
        root_identity=(
            "PR_VAPOR_ROOT"
            if phase == "VAPOR"
            else "PR_UNIQUE_SUPERCRITICAL_ROOT"
        ),
        model_id="test-pr-owner-v4",
        phase_caloric_schema_id="rfam.phase-caloric/v4",
        phase_caloric_options_fingerprint=digest("options"),
        ordered_species_count=species_count,
        ordered_species_fingerprint=digest("species"),
        energy_convention=(
            "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT"
        ),
        process_runtime_epoch=epoch,
        parent_generation=generation,
    )


def new_lifecycle(*bindings: object) -> object:
    epoch = bindings[0].process_runtime_epoch
    lifecycle = oracle.PhaseCaloricLifecycleOracle(epoch)
    for binding in bindings:
        lifecycle.register_parent(binding)
    return lifecycle


def token_dict(token: object) -> dict[str, object]:
    return dataclasses.asdict(token)


def zero_residuals() -> dict[str, float]:
    return {name: 0.0 for name in DOCUMENT["residual_formulas"]["required_residuals"]}


def canonical_request(token: object | None = None) -> dict[str, object]:
    binding = parent_binding("cell-A")
    lifecycle = new_lifecycle(binding)
    transaction_id = lifecycle.begin_transaction("cell-A")
    return {
        "kind": "PHASE_CALORIC_REQUEST_V4",
        "parent_owner_id": "cell-A",
        "committed_parent_q_fingerprint": binding.committed_parent_q_fingerprint,
        "parent_phase_payload_fingerprint": binding.parent_phase_payload_fingerprint,
        "parent_phase_identity": "VAPOR",
        "process_runtime_epoch": 7,
        "parent_generation": 11,
        "transaction_id": transaction_id,
        "phase_caloric_schema_id": "rfam.phase-caloric/v4",
        "phase_caloric_options_fingerprint": digest("options"),
        "model_id": "test-pr-owner-v4",
        "ordered_species_count": 2,
        "ordered_species_fingerprint": digest("species"),
        "energy_convention": "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT",
        "phase_density_kg_m3": 12.0,
        "phase_specific_internal_energy_j_kg": 250000.0,
        "phase_mass_fractions": [0.8, 0.2],
        "pressure_seed_pa": 2000000.0,
        "temperature_seed_k": 300.0,
        "previous_accepted_root_token": (
            None if token is None else token_dict(token)
        ),
    }


def canonical_success() -> dict[str, object]:
    binding = parent_binding("cell-A")
    lifecycle = new_lifecycle(binding)
    transaction_id = lifecycle.begin_transaction("cell-A")
    token = lifecycle.issue_initial_token("cell-A")
    return {
        "kind": "PHASE_CALORIC_SUCCESS_V4",
        "status": "SUCCESS",
        "parent_owner_id": "cell-A",
        "committed_parent_q_fingerprint": binding.committed_parent_q_fingerprint,
        "parent_phase_payload_fingerprint": binding.parent_phase_payload_fingerprint,
        "parent_phase_identity": "VAPOR",
        "process_runtime_epoch": 7,
        "parent_generation": 11,
        "transaction_id": transaction_id,
        "transaction_sequence": 1,
        "phase_caloric_schema_id": "rfam.phase-caloric/v4",
        "phase_caloric_options_fingerprint": digest("options"),
        "model_id": "test-pr-owner-v4",
        "ordered_species_count": 2,
        "ordered_species_fingerprint": digest("species"),
        "energy_convention": "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT",
        "pressure_pa": 2000000.0,
        "temperature_k": 300.0,
        "compressibility_z": 0.85,
        "root_identity": "PR_VAPOR_ROOT",
        "density_kg_m3": 12.0,
        "specific_internal_energy_j_kg": 250000.0,
        "specific_enthalpy_j_kg": 416666.6666666667,
        "mass_fractions": [0.8, 0.2],
        "cv_j_kg_k": 900.0,
        "dedt_rho_j_kg_k": 900.0,
        "dpdt_rho_pa_k": 7000.0,
        "dpdrho_t_pa_m3_kg": 150000.0,
        "residuals": zero_residuals(),
        "accepted_root_token": token_dict(token),
    }


def canonical_failure() -> dict[str, object]:
    return {
        "kind": "PHASE_CALORIC_FAILURE_V4",
        "status": "ROOT_LOSS",
        "failure_stage": "MID_RHS_RECOVERY",
        "parent_binding_digest": digest("parent-binding"),
        "process_runtime_epoch": 7,
        "parent_generation": 11,
        "transaction_id": digest("transaction"),
        "diagnostic_code": "explicit-root-disappeared",
        "residual_summary": {
            "all_finite": False,
            "maximum_normalized_residual": None,
            "failed_residual": "eos_root",
        },
    }


class Snapshot:
    def __init__(self) -> None:
        self._temporary = tempfile.TemporaryDirectory(
            prefix="rfam-phase-caloric-v4-"
        )
        self.root = Path(self._temporary.name)
        for relative in ARTIFACT_PATHS:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / relative, target)

    @property
    def contract(self) -> Path:
        return self.root / ARTIFACT_PATHS[0]

    def document(self) -> dict[str, object]:
        return json.loads(self.contract.read_text(encoding="ascii"))

    def write_document(self, document: dict[str, object]) -> str:
        raw = (json.dumps(document, indent=2) + "\n").encode("ascii")
        self.contract.write_bytes(raw)
        return sha256_bytes(raw)

    def close(self) -> None:
        self._temporary.cleanup()


class PhaseCaloricFeasibilityContractV4Tests(unittest.TestCase):
    checks = 0

    @classmethod
    def bump(cls, amount: int = 1) -> None:
        cls.checks += amount

    def assert_contract_error(self, callable_object, needle: str) -> None:
        with self.assertRaises(
            (validator.ContractError, oracle.LifecycleError)
        ) as caught:
            callable_object()
        self.assertIn(needle, str(caught.exception))
        self.bump()

    def test_canonical_gate_and_typed_objects_are_blocked_and_closed(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                "-B",
                str(VALIDATOR),
                "--root",
                str(ROOT),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PHASE_CALORIC_FEASIBILITY_CONTRACT_V4_BLOCKED", result.stdout)
        for marker in (
            "runtime_owner=0",
            "thermo_recovery=0",
            "two_phase_chemistry=0",
            "numerical_limits_validated=0",
            "production_ready=0",
        ):
            self.assertIn(marker, result.stdout)
        validator.validate_typed_instance(DOCUMENT, "request", canonical_request())
        validator.validate_typed_instance(DOCUMENT, "success", canonical_success())
        validator.validate_typed_instance(DOCUMENT, "failure", canonical_failure())
        self.bump(4)

    def test_json_path_and_section_rebinding_fail_closed(self) -> None:
        alternate = ROOT / "tests/reference/phase_caloric_feasibility_contract_v3.json"
        self.assert_contract_error(
            lambda: validator.validate(ROOT, alternate),
            "contract-canonical-path-mismatch",
        )

        for raw, needle in (
            (
                CONTRACT.read_bytes().replace(
                    b'  "scope": ', b'  "scope": "runtime",\n  "scope": ', 1
                ),
                "duplicate-json-key:scope",
            ),
            (
                CONTRACT.read_bytes().replace(
                    b'"phase_density_relative_displacement_max": 0.0',
                    b'"phase_density_relative_displacement_max": NaN',
                    1,
                ),
                "nonfinite-json-number:NaN",
            ),
        ):
            snapshot = Snapshot()
            try:
                snapshot.contract.write_bytes(raw)
                self.assert_contract_error(
                    lambda snapshot=snapshot, raw=raw: validator.validate(
                        snapshot.root,
                        snapshot.contract,
                        expected_contract_sha256=sha256_bytes(raw),
                        enforce_runtime_path=False,
                        import_oracle=False,
                    ),
                    needle,
                )
            finally:
                snapshot.close()

        mutations = (
            (
                lambda d: d["approvals"].__setitem__("production_ready", True),
                "approval-not-false:production_ready",
            ),
            (
                lambda d: d["lifecycle_ownership"][
                    "parent_generation"
                ].__setitem__("owner", "global-owner"),
                "section-sha256-mismatch:lifecycle_ownership",
            ),
            (
                lambda d: d["event_order"]["rollback"].reverse(),
                "section-sha256-mismatch:event_order",
            ),
            (
                lambda d: d["residual_formulas"].__setitem__(
                    "relative_formula", "undefined"
                ),
                "section-sha256-mismatch:residual_formulas",
            ),
            (
                lambda d: d["applicability_envelope"].__setitem__(
                    "numerical_limit_validation_status", "VALIDATED"
                ),
                "section-sha256-mismatch:applicability_envelope",
            ),
            (
                lambda d: d["publication_boundary"].__setitem__(
                    "phase_local_result_generic_solver_admission", "allowed"
                ),
                "section-sha256-mismatch:publication_boundary",
            ),
        )
        for mutation, needle in mutations:
            snapshot = Snapshot()
            try:
                document = snapshot.document()
                mutation(document)
                expected = snapshot.write_document(document)
                self.assert_contract_error(
                    lambda snapshot=snapshot, expected=expected: validator.validate(
                        snapshot.root,
                        snapshot.contract,
                        expected_contract_sha256=expected,
                        enforce_runtime_path=False,
                        import_oracle=False,
                    ),
                    needle,
                )
            finally:
                snapshot.close()

    def test_request_success_failure_and_nested_types_reject_aliases(self) -> None:
        cases: list[tuple[str, dict[str, object], str]] = []
        request = canonical_request()
        request.pop("parent_phase_identity")
        cases.append(("request", request, "typed-instance-missing"))
        request = canonical_request()
        request["process_runtime_epoch"] = True
        cases.append(("request", request, "invalid-uint64"))
        request = canonical_request()
        request["phase_density_kg_m3"] = 12
        cases.append(("request", request, "invalid-binary64"))
        request = canonical_request()
        request["ordered_species_count"] = True
        cases.append(("request", request, "invalid-uint64"))
        request = canonical_request()
        request["species_count"] = request["ordered_species_count"]
        cases.append(("request", request, "typed-instance-unknown"))

        success = canonical_success()
        success["parent_q_fingerprint"] = success[
            "committed_parent_q_fingerprint"
        ]
        cases.append(("success", success, "typed-instance-unknown"))
        success = canonical_success()
        success["root_identity"] = "PR_UNIQUE_SUPERCRITICAL_ROOT"
        cases.append(("success", success, "typed-phase-root-mismatch"))
        success = canonical_success()
        success["accepted_root_token"]["state"] = {}
        cases.append(("success", success, "typed-instance-unknown:live_root_token"))
        success = canonical_success()
        success["accepted_root_token"]["parent_owner_id"] = "cell-B"
        cases.append(("success", success, "typed-live-token-binding:success"))
        success = canonical_success()
        success["accepted_root_token"]["transaction_sequence"] = 2
        cases.append(("success", success, "typed-live-token-binding:success"))
        success = canonical_success()
        success["accepted_root_token"]["token_digest"] = digest("wrong-token")
        cases.append(("success", success, "typed-live-token-digest:success"))
        success = canonical_success()
        success["ordered_species_count"] = "2"
        cases.append(("success", success, "invalid-uint64"))

        other_binding = parent_binding("cell-B")
        other_lifecycle = new_lifecycle(other_binding)
        other_lifecycle.begin_transaction("cell-B")
        other_token = other_lifecycle.issue_initial_token("cell-B")
        request = canonical_request(other_token)
        cases.append(("request", request, "typed-live-token-binding:request"))

        hint_binding = parent_binding("cell-A")
        hint_lifecycle = new_lifecycle(hint_binding)
        hint_lifecycle.begin_transaction("cell-A")
        hint_token = hint_lifecycle.issue_initial_token("cell-A")
        hint = dataclasses.asdict(hint_lifecycle.make_restart_hint(hint_token))
        hint["hint_digest"] = digest("wrong-hint")
        cases.append(("restart_hint", hint, "typed-restart-hint-digest"))

        for alias in (
            "state",
            "accepted_root_token",
            "token",
            "quality",
            "sound_speed",
            "publication_candidate",
        ):
            failure = canonical_failure()
            failure[alias] = None
            cases.append(("failure", failure, "typed-instance-unknown:failure"))

        for schema_name, instance, needle in cases:
            with self.subTest(schema=schema_name, needle=needle):
                self.assert_contract_error(
                    lambda schema_name=schema_name, instance=instance: (
                        validator.validate_typed_instance(
                            DOCUMENT, schema_name, instance
                        )
                    ),
                    needle,
                )

    def test_canonical_fingerprint_coverage_and_binary64_encoding(self) -> None:
        q_fields = DOCUMENT["canonical_serialization"][
            "committed_parent_q_fingerprint"
        ]["required_fields"]
        q_values: dict[str, object] = {
            "schema": "rfam.committed-parent-q/v4",
            "parent_owner_id": "cell-A",
            "rho_kg_m3_f64be": f64be(10.0),
            "rho_u_kg_m2_s_f64be_xyz": [f64be(1.0), f64be(2.0), f64be(3.0)],
            "rho_e_j_m3_f64be": f64be(2500000.0),
            "rho_y_kg_m3_f64be_ordered_species": [f64be(8.0), f64be(2.0)],
            "rho_z_kg_m3_f64be_ordered_scalars": [f64be(0.0)],
            "rho_c_kg_m3_f64be_ordered_scalars": [f64be(0.0)],
            "ordered_species_count": 2,
            "ordered_species_fingerprint": digest("species"),
            "energy_convention": "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT",
        }
        q_payload = {name: q_values[name] for name in q_fields}
        q_fingerprint = validator.fingerprint_payload(
            DOCUMENT, oracle, "committed_parent_q_fingerprint", q_payload
        )
        self.assertRegex(q_fingerprint, r"^[0-9a-f]{64}$")

        phase_fields = DOCUMENT["canonical_serialization"][
            "parent_phase_payload_fingerprint"
        ]["required_fields"]
        phase_values: dict[str, object] = {
            "schema": "rfam.parent-phase-payload/v4",
            "parent_owner_id": "cell-A",
            "parent_phase_identity": "VAPOR",
            "pressure_pa_f64be": f64be(2000000.0),
            "temperature_k_f64be": f64be(300.0),
            "compressibility_z_f64be": f64be(0.85),
            "density_kg_m3_f64be": f64be(12.0),
            "specific_internal_energy_j_kg_f64be": f64be(250000.0),
            "specific_enthalpy_j_kg_f64be": f64be(416666.6666666667),
            "mass_fractions_f64be_ordered_species": [f64be(0.8), f64be(0.2)],
            "root_identity": "PR_VAPOR_ROOT",
            "model_id": "test-pr-owner-v4",
            "phase_caloric_schema_id": "rfam.phase-caloric/v4",
            "phase_caloric_options_fingerprint": digest("options"),
            "ordered_species_count": 2,
            "ordered_species_fingerprint": digest("species"),
            "energy_convention": "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT",
        }
        phase_payload = {name: phase_values[name] for name in phase_fields}
        phase_fingerprint = validator.fingerprint_payload(
            DOCUMENT, oracle, "parent_phase_payload_fingerprint", phase_payload
        )
        self.assertNotEqual(q_fingerprint, phase_fingerprint)

        direct_a = oracle.canonical_fingerprint(
            "rfam.test.canonical/v4", {"b": 2, "a": 1}
        )
        direct_b = oracle.canonical_fingerprint(
            "rfam.test.canonical/v4", {"a": 1, "b": 2}
        )
        self.assertEqual(direct_a, direct_b)
        reversed_q_fingerprint = validator.fingerprint_payload(
            DOCUMENT,
            oracle,
            "committed_parent_q_fingerprint",
            dict(reversed(list(q_payload.items()))),
        )
        self.assertEqual(q_fingerprint, reversed_q_fingerprint)
        missing_q = dict(q_payload)
        missing_q.pop("rho_e_j_m3_f64be")
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "committed_parent_q_fingerprint",
                missing_q,
            ),
            "fingerprint-field-coverage",
        )
        extra_q = dict(q_payload)
        extra_q["rho_quality_f64be"] = f64be(0.0)
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "committed_parent_q_fingerprint",
                extra_q,
            ),
            "fingerprint-field-coverage",
        )
        invalid_phase = dict(phase_payload)
        invalid_phase["root_identity"] = "PR_UNIQUE_SUPERCRITICAL_ROOT"
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "parent_phase_payload_fingerprint",
                invalid_phase,
            ),
            "fingerprint-phase-root-mismatch",
        )
        invalid_q = dict(q_payload)
        invalid_q["rho_kg_m3_f64be"] = "10.0"
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "committed_parent_q_fingerprint",
                invalid_q,
            ),
            "fingerprint-f64be",
        )
        invalid_q = dict(q_payload)
        invalid_q["rho_e_j_m3_f64be"] = "7ff8000000000000"
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "committed_parent_q_fingerprint",
                invalid_q,
            ),
            "fingerprint-f64be-nonfinite",
        )
        invalid_phase = dict(phase_payload)
        invalid_phase["mass_fractions_f64be_ordered_species"] = [
            f64be(0.8),
            "7ff0000000000000",
        ]
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "parent_phase_payload_fingerprint",
                invalid_phase,
            ),
            "fingerprint-f64be-nonfinite",
        )
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT, oracle, "unknown_fingerprint", {}
            ),
            "unknown-fingerprint-name",
        )
        malformed_q_count = dict(q_payload)
        malformed_q_count["ordered_species_count"] = False
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "committed_parent_q_fingerprint",
                malformed_q_count,
            ),
            "invalid-uint64",
        )
        malformed_phase_count = dict(phase_payload)
        malformed_phase_count["ordered_species_count"] = "2"
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "parent_phase_payload_fingerprint",
                malformed_phase_count,
            ),
            "invalid-uint64",
        )
        self.bump(10)

    def test_canonical_json_serialize_parse_and_fingerprint_round_trip(self) -> None:
        payload = {
            "zeta": [3, 2, 1],
            "alpha": {"z": "last", "a": 1},
        }
        serialized = oracle.canonical_json_bytes(payload)
        self.assertEqual(
            serialized,
            b'{"alpha":{"a":1,"z":"last"},"zeta":[3,2,1]}',
        )
        parsed = oracle.parse_canonical_json_bytes(serialized)
        self.assertEqual(oracle.canonical_json_bytes(parsed), serialized)
        fingerprint = oracle.canonical_fingerprint(
            "rfam.test.canonical-round-trip/v4", payload
        )
        self.assertEqual(
            fingerprint,
            oracle.canonical_fingerprint(
                "rfam.test.canonical-round-trip/v4", parsed
            ),
        )
        self.assert_contract_error(
            lambda: oracle.parse_canonical_json_bytes(
                b'{"zeta":[3,2,1],"alpha":{"a":1,"z":"last"}}'
            ),
            "canonical-json-bytes-not-normalized",
        )
        self.assert_contract_error(
            lambda: oracle.parse_canonical_json_bytes(b'{"a":1,"a":1}'),
            "duplicate-canonical-json-key:a",
        )
        self.bump(3)

    def test_ordered_species_count_binds_all_four_species_vectors(self) -> None:
        q_fields = DOCUMENT["canonical_serialization"][
            "committed_parent_q_fingerprint"
        ]["required_fields"]
        q_values: dict[str, object] = {
            "schema": "rfam.committed-parent-q/v4",
            "parent_owner_id": "cell-A",
            "rho_kg_m3_f64be": f64be(10.0),
            "rho_u_kg_m2_s_f64be_xyz": [f64be(1.0), f64be(2.0), f64be(3.0)],
            "rho_e_j_m3_f64be": f64be(2500000.0),
            "rho_y_kg_m3_f64be_ordered_species": [f64be(10.0)],
            "rho_z_kg_m3_f64be_ordered_scalars": [f64be(0.0)],
            "rho_c_kg_m3_f64be_ordered_scalars": [f64be(0.0)],
            "ordered_species_count": 2,
            "ordered_species_fingerprint": digest("species"),
            "energy_convention": "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT",
        }
        q_payload = {name: q_values[name] for name in q_fields}

        phase_fields = DOCUMENT["canonical_serialization"][
            "parent_phase_payload_fingerprint"
        ]["required_fields"]
        phase_values: dict[str, object] = {
            "schema": "rfam.parent-phase-payload/v4",
            "parent_owner_id": "cell-A",
            "parent_phase_identity": "VAPOR",
            "pressure_pa_f64be": f64be(2000000.0),
            "temperature_k_f64be": f64be(300.0),
            "compressibility_z_f64be": f64be(0.85),
            "density_kg_m3_f64be": f64be(12.0),
            "specific_internal_energy_j_kg_f64be": f64be(250000.0),
            "specific_enthalpy_j_kg_f64be": f64be(416666.6666666667),
            "mass_fractions_f64be_ordered_species": [f64be(1.0)],
            "root_identity": "PR_VAPOR_ROOT",
            "model_id": "test-pr-owner-v4",
            "phase_caloric_schema_id": "rfam.phase-caloric/v4",
            "phase_caloric_options_fingerprint": digest("options"),
            "ordered_species_count": 2,
            "ordered_species_fingerprint": digest("species"),
            "energy_convention": "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT",
        }
        phase_payload = {name: phase_values[name] for name in phase_fields}

        request = canonical_request()
        request["phase_mass_fractions"] = [1.0]
        success = canonical_success()
        success["mass_fractions"] = [1.0]
        cases = (
            lambda: validator.fingerprint_payload(
                DOCUMENT, oracle, "committed_parent_q_fingerprint", q_payload
            ),
            lambda: validator.fingerprint_payload(
                DOCUMENT, oracle, "parent_phase_payload_fingerprint", phase_payload
            ),
            lambda: validator.validate_typed_instance(DOCUMENT, "request", request),
            lambda: validator.validate_typed_instance(DOCUMENT, "success", success),
        )
        for case in cases:
            self.assert_contract_error(case, "ordered-species-count")
        aliased_q = dict(q_payload)
        aliased_q["species_count"] = aliased_q.pop("ordered_species_count")
        self.assert_contract_error(
            lambda: validator.fingerprint_payload(
                DOCUMENT,
                oracle,
                "committed_parent_q_fingerprint",
                aliased_q,
            ),
            "fingerprint-field-coverage",
        )
        self.assert_contract_error(
            lambda: parent_binding("cell-count-alias", species_count=True),
            "invalid-uint64:ordered-species-count",
        )
        self.bump(6)

    def test_rollback_and_final_reflash_invalidate_stale_tokens(self) -> None:
        binding = parent_binding("cell-A")
        lifecycle = new_lifecycle(binding)
        lifecycle.begin_transaction("cell-A")
        token1 = lifecycle.issue_initial_token("cell-A")
        token2 = lifecycle.continue_from_token(token1)
        self.assertFalse(lifecycle.accepts_live_token(token1))
        self.assertTrue(lifecycle.accepts_live_token(token2))
        rolled_back = lifecycle.rollback(
            "cell-A",
            restored_committed_parent_q_fingerprint=(
                binding.committed_parent_q_fingerprint
            ),
            restored_parent_phase_payload_fingerprint=(
                binding.parent_phase_payload_fingerprint
            ),
        )
        self.assertEqual(rolled_back.parent_generation, 12)
        self.assertEqual(
            rolled_back.committed_parent_q_fingerprint,
            binding.committed_parent_q_fingerprint,
        )
        self.assertFalse(lifecycle.accepts_live_token(token2))

        lifecycle.begin_transaction("cell-A")
        token3 = lifecycle.issue_initial_token("cell-A")
        published = lifecycle.final_reflash_publish(
            "cell-A",
            committed_parent_q_fingerprint=digest("cell-A-new-q"),
            parent_phase_payload_fingerprint=digest("cell-A-new-phase"),
            parent_phase_identity="SUPERCRITICAL",
            root_identity="PR_UNIQUE_SUPERCRITICAL_ROOT",
        )
        self.assertEqual(published.parent_generation, 13)
        self.assertFalse(lifecycle.accepts_live_token(token3))
        self.assert_contract_error(
            lambda: lifecycle.issue_initial_token("cell-A"),
            "no-active-parent-transaction",
        )
        self.bump(7)

    def test_restart_hint_requires_fresh_epoch_generation_replay_and_token(self) -> None:
        binding = parent_binding("cell-A")
        lifecycle = new_lifecycle(binding)
        lifecycle.begin_transaction("cell-A")
        token = lifecycle.issue_initial_token("cell-A")
        hint = lifecycle.make_restart_hint(token)
        restarted = lifecycle.restart()
        self.assertFalse(restarted.accepts_live_token(token))
        self.assertFalse(restarted.accepts_live_token(hint))
        fresh = restarted.revalidate_restart_hint(hint)
        self.assertTrue(restarted.accepts_live_token(fresh))
        self.assertEqual(fresh.process_runtime_epoch, 8)
        self.assertEqual(fresh.parent_generation, 12)
        self.assertEqual(fresh.transaction_sequence, 1)
        self.assertNotEqual(fresh.token_digest, token.token_digest)

        second_restart = restarted.restart()
        self.assertFalse(second_restart.accepts_live_token(token))
        self.assertFalse(second_restart.accepts_live_token(fresh))
        reused = second_restart.revalidate_restart_hint(hint)
        self.assertTrue(second_restart.accepts_live_token(reused))
        self.assertEqual(reused.process_runtime_epoch, 9)
        self.assertEqual(reused.parent_generation, 13)
        self.assertEqual(reused.transaction_sequence, 1)
        self.assertNotEqual(reused.token_digest, token.token_digest)
        self.assertNotEqual(reused.token_digest, fresh.token_digest)

        mismatched = dataclasses.replace(
            hint,
            committed_parent_q_fingerprint=digest("other-q"),
        )
        mismatch_restart = lifecycle.restart()
        self.assert_contract_error(
            lambda: mismatch_restart.revalidate_restart_hint(mismatched),
            "restart-hint-parent-revalidation-failed",
        )
        self.bump(14)

    def test_two_parent_independence_order_determinism_and_overflow(self) -> None:
        a = parent_binding("cell-A")
        b = parent_binding("cell-B")

        first = new_lifecycle(a, b)
        first.begin_transaction("cell-A")
        token_a_first = first.issue_initial_token("cell-A")
        first.begin_transaction("cell-B")
        token_b_first = first.issue_initial_token("cell-B")

        second = new_lifecycle(a, b)
        second.begin_transaction("cell-B")
        token_b_second = second.issue_initial_token("cell-B")
        second.begin_transaction("cell-A")
        token_a_second = second.issue_initial_token("cell-A")
        self.assertEqual(token_a_first.token_digest, token_a_second.token_digest)
        self.assertEqual(token_b_first.token_digest, token_b_second.token_digest)

        first.rollback(
            "cell-A",
            restored_committed_parent_q_fingerprint=a.committed_parent_q_fingerprint,
            restored_parent_phase_payload_fingerprint=(
                a.parent_phase_payload_fingerprint
            ),
        )
        self.assertFalse(first.accepts_live_token(token_a_first))
        self.assertTrue(first.accepts_live_token(token_b_first))
        self.assertEqual(first.parent("cell-B").parent_generation, 11)
        self.assert_contract_error(
            lambda: first.begin_transaction("cell-B"),
            "parent-transaction-already-active",
        )

        epoch_max = oracle.PhaseCaloricLifecycleOracle(oracle.UINT64_MAX)
        epoch_max.register_parent(
            parent_binding("cell-max", epoch=oracle.UINT64_MAX)
        )
        self.assert_contract_error(
            epoch_max.restart, "uint64-overflow:process-runtime-epoch"
        )

        generation_max = parent_binding(
            "cell-generation-max", generation=oracle.UINT64_MAX
        )
        generation_lifecycle = new_lifecycle(generation_max)
        generation_lifecycle.begin_transaction("cell-generation-max")
        generation_token = generation_lifecycle.issue_initial_token(
            "cell-generation-max"
        )
        self.assert_contract_error(
            lambda: generation_lifecycle.rollback(
                "cell-generation-max",
                restored_committed_parent_q_fingerprint=(
                    generation_max.committed_parent_q_fingerprint
                ),
                restored_parent_phase_payload_fingerprint=(
                    generation_max.parent_phase_payload_fingerprint
                ),
            ),
            "uint64-overflow:parent-generation",
        )
        self.assertTrue(generation_lifecycle.accepts_live_token(generation_token))
        continued = generation_lifecycle.continue_from_token(generation_token)
        self.assertEqual(continued.transaction_sequence, 2)
        self.assertTrue(generation_lifecycle.accepts_live_token(continued))

        sequence_binding = parent_binding("cell-sequence-max")
        sequence_lifecycle = new_lifecycle(sequence_binding)
        sequence_id = sequence_lifecycle.begin_transaction("cell-sequence-max")
        sequence_transaction = sequence_lifecycle._transactions[
            "cell-sequence-max"
        ]
        sequence_transaction.current_sequence = oracle.UINT64_MAX
        sequence_digest = sequence_lifecycle._token_digest(
            sequence_binding, sequence_transaction, oracle.UINT64_MAX
        )
        sequence_transaction.current_token_digest = sequence_digest
        sequence_token = oracle.LiveToken(
            parent_owner_id="cell-sequence-max",
            transaction_id=sequence_id,
            process_runtime_epoch=sequence_binding.process_runtime_epoch,
            parent_generation=sequence_binding.parent_generation,
            transaction_sequence=oracle.UINT64_MAX,
            parent_binding_digest=sequence_binding.identity_digest(),
            root_identity=sequence_binding.root_identity,
            token_digest=sequence_digest,
        )
        self.assertTrue(sequence_lifecycle.accepts_live_token(sequence_token))
        self.assert_contract_error(
            lambda: sequence_lifecycle.continue_from_token(sequence_token),
            "uint64-overflow:transaction-sequence:cell-sequence-max",
        )
        self.assertTrue(sequence_lifecycle.accepts_live_token(sequence_token))
        self.bump(11)

    def test_artifact_path_hash_rebinding_and_substitute_oracle_fail_preload(self) -> None:
        path_snapshot = Snapshot()
        try:
            document = path_snapshot.document()
            paths = document["artifact_binding"]["canonical_relative_paths"]
            paths["lifecycle_oracle"], paths["documentation"] = (
                paths["documentation"],
                paths["lifecycle_oracle"],
            )
            expected = path_snapshot.write_document(document)
            self.assert_contract_error(
                lambda: validator.validate(
                    path_snapshot.root,
                    path_snapshot.contract,
                    expected_contract_sha256=expected,
                    enforce_runtime_path=False,
                    import_oracle=False,
                ),
                "section-sha256-mismatch:artifact_binding",
            )
        finally:
            path_snapshot.close()

        oracle_snapshot = Snapshot()
        marker = oracle_snapshot.root / "substitute-oracle-executed"
        try:
            substitute = (
                "from pathlib import Path\n"
                f"Path({str(marker)!r}).write_text('executed')\n"
            ).encode("ascii")
            oracle_path = oracle_snapshot.root / ARTIFACT_PATHS[1]
            oracle_path.write_bytes(substitute)
            document = oracle_snapshot.document()
            document["artifact_binding"]["lifecycle_oracle_sha256"] = (
                sha256_bytes(substitute)
            )
            expected = oracle_snapshot.write_document(document)
            self.assert_contract_error(
                lambda: validator.validate(
                    oracle_snapshot.root,
                    oracle_snapshot.contract,
                    expected_contract_sha256=expected,
                    enforce_runtime_path=False,
                    import_oracle=True,
                ),
                "section-sha256-mismatch:artifact_binding",
            )
            self.assertFalse(marker.exists())
            self.bump()
        finally:
            oracle_snapshot.close()

        docs_snapshot = Snapshot()
        try:
            docs_path = docs_snapshot.root / ARTIFACT_PATHS[4]
            docs_path.write_bytes(docs_path.read_bytes() + b"tampered\n")
            document = docs_snapshot.document()
            document["artifact_binding"]["documentation_sha256"] = sha256_bytes(
                docs_path.read_bytes()
            )
            expected = docs_snapshot.write_document(document)
            self.assert_contract_error(
                lambda: validator.validate(
                    docs_snapshot.root,
                    docs_snapshot.contract,
                    expected_contract_sha256=expected,
                    enforce_runtime_path=False,
                    import_oracle=False,
                ),
                "section-sha256-mismatch:artifact_binding",
            )
        finally:
            docs_snapshot.close()


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        PhaseCaloricFeasibilityContractV4Tests
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(
            "PHASE_CALORIC_FEASIBILITY_CONTRACT_V4_TEST_PASS checks="
            + str(PhaseCaloricFeasibilityContractV4Tests.checks)
        )
    raise SystemExit(0 if result.wasSuccessful() else 1)
