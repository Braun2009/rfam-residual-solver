#!/usr/bin/env python3
"""Focused fail-close tests for the O2/C3H8 PR kij pre-fit contract."""

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


REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
TOOL = REPO / "tools/validate_kij_fit_contract.py"
CONTRACT = REPO / "tests/reference/kij_fit_contract_v1.json"
BASE_TEXT = CONTRACT.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)

sys.path.insert(0, str(REPO / "tools"))
import validate_kij_fit_contract as validator  # noqa: E402


class KijFitContractTests(unittest.TestCase):
    checks = 0

    def _fixture(self, root: Path) -> Path:
        for item in validator.EXPECTED_UPSTREAM:
            source = WORKSPACE / item["path"]
            target = root / item["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        target = root / "realFluidAllMach/tests/reference/kij_fit_contract_v1.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        return target

    def _run(
        self,
        mutate: Callable[[dict[str, Any]], None] | None = None,
        raw_text: str | None = None,
        mutate_workspace: Callable[[Path], None] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-fit-contract-") as temporary:
            root = Path(temporary)
            contract = self._fixture(root)
            document = copy.deepcopy(BASE_DOCUMENT)
            if mutate is not None:
                mutate(document)
            with contract.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(
                    raw_text
                    if raw_text is not None
                    else json.dumps(document, indent=2) + "\n"
                )
            if mutate_workspace is not None:
                mutate_workspace(root)
            return subprocess.run(
                [
                    sys.executable,
                    str(TOOL),
                    "--root",
                    str(root),
                    "--contract",
                    str(contract),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

    def _require_failure(
        self, result: subprocess.CompletedProcess[str], reason: str
    ) -> None:
        self.__class__.checks += 1
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, output)
        self.assertIn(validator.FAIL_MARKER, result.stdout)
        self.assertIn(reason, result.stdout)
        self.assertNotIn(validator.BLOCKED_MARKER, output)
        self.assertNotIn("Traceback", output)

    def test_frozen_contract_is_consistent_but_blocked(self) -> None:
        result = self._run()
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(validator.BLOCKED_MARKER, result.stdout)
        self.assertIn("raw_dataset_selected=0", result.stdout)
        self.assertIn("fit_executed=0", result.stdout)
        self.assertIn("validated_kij=0", result.stdout)
        self.assertIn("production_ready=0", result.stdout)

    def test_duplicate_approval_key_is_rejected(self) -> None:
        needle = '    "fit_executed": false,'
        replacement = '    "fit_executed": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:fit_executed",
        )

    def test_duplicate_nested_fit_key_is_rejected(self) -> None:
        needle = '    "status": "not-selected-not-executed",'
        replacement = '    "status": "executed",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:status",
        )

    def test_noncanonical_json_is_rejected(self) -> None:
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace("{\n", "{ \n", 1)),
            "contract-not-canonical-json",
        )

    def test_contract_version_bool_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("contract_version", True)),
            "contract-type-mismatch:root.contract_version",
        )

    def test_unknown_root_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("unreviewed", True)),
            "contract-unknown-keys:root:unreviewed",
        )

    def test_fit_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["approvals"].__setitem__(
                    "fit_executed", True
                )
            ),
            "contract-value-mismatch:root.approvals.fit_executed",
        )

    def test_runtime_value_cannot_be_claimed(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["pair_identity"].__setitem__(
                    "parameter_state", "runtime-value-0.055"
                )
            ),
            "contract-value-mismatch:root.pair_identity.parameter_state",
        )

    def test_pair_order_cannot_be_reclassified(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["pair_identity"].__setitem__(
                    "canonical_pair_key", "O2|C3H8"
                )
            ),
            "contract-value-mismatch:root.pair_identity.canonical_pair_key",
        )

    def test_symmetry_requirement_cannot_be_weakened(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["pair_identity"].__setitem__(
                    "symmetry_requirement", "not-required"
                )
            ),
            "contract-value-mismatch:root.pair_identity.symmetry_requirement",
        )

    def test_source_candidate_cannot_be_promoted_to_fit_input(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["data_provenance"].__setitem__(
                    "source_candidate_role", "training"
                )
            ),
            "contract-value-mismatch:root.data_provenance.source_candidate_role",
        )

    def test_unbound_reference_backend_values_cannot_be_allowed(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["data_provenance"].__setitem__(
                    "no_unbound_reference_backend_values", False
                )
            ),
            "contract-value-mismatch:root.data_provenance.no_unbound_reference_backend_values",
        )

    def test_rights_cannot_be_marked_resolved(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["rights"].__setitem__(
                    "fit_and_redistribution_permitted", True
                )
            ),
            "contract-value-mismatch:root.rights.fit_and_redistribution_permitted",
        )

    def test_phase_label_requirement_cannot_be_removed(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["measurement_schema"]["phase_label_enum"].remove("solid-risk")

        self._require_failure(
            self._run(mutate),
            "contract-list-length-mismatch:root.measurement_schema.phase_label_enum",
        )

    def test_vapor_composition_requirement_cannot_be_weakened(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["measurement_schema"].__setitem__(
                    "vle_fit_requires_both-liquid-and-vapor-composition", False
                )
            ),
            "contract-value-mismatch:root.measurement_schema.vle_fit_requires_both-liquid-and-vapor-composition",
        )

    def test_uncertainty_zero_substitution_cannot_be_allowed(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["measurement_schema"].__setitem__(
                    "missing_uncertainty_policy", "treat-as-zero"
                )
            ),
            "contract-value-mismatch:root.measurement_schema.missing_uncertainty_policy",
        )

    def test_pressure_unit_cannot_be_relabelled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["measurement_schema"].__setitem__(
                    "pressure_unit", "kPa"
                )
            ),
            "contract-value-mismatch:root.measurement_schema.pressure_unit",
        )

    def test_training_holdout_split_cannot_be_inserted_without_new_freeze(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["partition"].__setitem__(
                    "training_dataset_ids", ["dataset-a"]
                )
            ),
            "contract-list-length-mismatch:root.partition.training_dataset_ids",
        )

    def test_holdout_requirement_cannot_be_disabled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["partition"].__setitem__(
                    "holdout_required_before_validation_claim", False
                )
            ),
            "contract-value-mismatch:root.partition.holdout_required_before_validation_claim",
        )

    def test_fit_form_cannot_be_selected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["fit_definition"].__setitem__(
                    "selected_form_id", "constant"
                )
            ),
            "contract-type-mismatch:root.fit_definition.selected_form_id",
        )

    def test_legacy_value_cannot_be_reintroduced(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["fit_definition"].__setitem__(
                    "parameter_values", [0.055]
                )
            ),
            "contract-list-length-mismatch:root.fit_definition.parameter_values",
        )

    def test_envelope_cannot_be_claimed_without_new_freeze(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["applicability_envelope"].__setitem__(
                    "validated_envelope", {"temperature_k": [110.0, 140.0]}
                )
            ),
            "contract-type-mismatch:root.applicability_envelope.validated_envelope",
        )

    def test_missing_pair_policy_cannot_fall_back_to_single_root(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["applicability_envelope"].__setitem__(
                    "outside_envelope_policy", "single-root-fallback"
                )
            ),
            "contract-value-mismatch:root.applicability_envelope.outside_envelope_policy",
        )

    def test_uncertainty_cannot_be_claimed_available(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["uncertainty"].__setitem__(
                    "parameter_uncertainty", {"c0": 0.001}
                )
            ),
            "contract-type-mismatch:root.uncertainty.parameter_uncertainty",
        )

    def test_nonlinear_sensitivity_requirement_cannot_be_removed(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["nonlinear_sensitivity"]["required_observables"].remove(
                "mass_quality"
            )

        self._require_failure(
            self._run(mutate),
            "contract-list-length-mismatch:root.nonlinear_sensitivity.required_observables",
        )

    def test_sensitivity_cannot_be_claimed_complete(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["nonlinear_sensitivity"].__setitem__(
                    "sensitivity_complete", True
                )
            ),
            "contract-value-mismatch:root.nonlinear_sensitivity.sensitivity_complete",
        )

    def test_acceptance_thresholds_cannot_be_erased(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["acceptance"]["required_before_fit"].pop()

        self._require_failure(
            self._run(mutate),
            "contract-list-length-mismatch:root.acceptance.required_before_fit",
        )

    def test_cantera_cannot_be_reclassified_as_kij_authority(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["limitations"][2] = "Cantera validates kij"

        self._require_failure(
            self._run(mutate),
            "contract-value-mismatch:root.limitations[2]",
        )

    def test_upstream_path_and_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first, second = document["upstream"]
            first["path"], second["path"] = second["path"], first["path"]
            first["sha256"], second["sha256"] = second["sha256"], first["sha256"]

        self._require_failure(
            self._run(mutate), "contract-value-mismatch:root.upstream[0].path"
        )

    def test_live_upstream_tamper_is_rejected(self) -> None:
        relative = validator.EXPECTED_UPSTREAM[0]["path"]

        def mutate_workspace(root: Path) -> None:
            path = root / relative
            path.write_bytes(path.read_bytes() + b"\n")

        self._require_failure(
            self._run(mutate_workspace=mutate_workspace),
            "upstream-sha256-mismatch:wp3_sensitivity_contract",
        )

    def test_self_consistent_upstream_rehash_is_rejected(self) -> None:
        relative = validator.EXPECTED_UPSTREAM[0]["path"]
        replacement = b"self-consistent-but-unreviewed-upstream\n"

        def mutate(document: dict[str, Any]) -> None:
            import hashlib

            document["upstream"][0]["sha256"] = hashlib.sha256(replacement).hexdigest()

        def mutate_workspace(root: Path) -> None:
            (root / relative).write_bytes(replacement)

        self._require_failure(
            self._run(mutate=mutate, mutate_workspace=mutate_workspace),
            "contract-value-mismatch:root.upstream[0].sha256",
        )

    def test_path_escape_is_rejected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["upstream"][0].__setitem__(
                    "path", "../outside"
                )
            ),
            "contract-value-mismatch:root.upstream[0].path",
        )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(KijFitContractTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"RFAM_KIJ_FIT_CONTRACT_TEST_PASS checks={KijFitContractTests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
