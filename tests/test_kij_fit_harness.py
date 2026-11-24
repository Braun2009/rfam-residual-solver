#!/usr/bin/env python3
"""Focused fail-close tests for the O2/C3H8 kij dry-run harness."""

from __future__ import annotations

import copy
import csv
import hashlib
import io
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
TOOL = REPO / "tools/kij_fit_harness.py"
MANIFEST = REPO / "tests/reference/kij_fit_harness_synthetic_manifest_v1.json"
BASE_TEXT = MANIFEST.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)
DB_SOURCE = REPO / "src/thermo/MixingRules/BinaryInteractionDB.C"

sys.path.insert(0, str(REPO / "tools"))
import kij_fit_harness as harness  # noqa: E402


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class KijFitHarnessTests(unittest.TestCase):
    checks = 0

    def _fixture(self, root: Path) -> Path:
        paths = [
            BASE_DOCUMENT["rights"]["terms_artifact"]["path"],
            *(item["path"] for item in BASE_DOCUMENT["artifacts"]),
        ]
        for relative in paths:
            source = WORKSPACE / relative
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        manifest = root / "realFluidAllMach/tests/reference/kij_fit_harness_synthetic_manifest_v1.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        return manifest

    @staticmethod
    def _rewrite_csv(
        root: Path,
        document: dict[str, Any],
        artifact_index: int,
        mutate: Callable[[list[list[str]]], None],
    ) -> None:
        authority = document["artifacts"][artifact_index]
        path = root / authority["path"]
        rows = list(csv.reader(io.StringIO(path.read_text(encoding="ascii")), strict=True))
        mutate(rows)
        stream = io.StringIO(newline="")
        csv.writer(stream, lineterminator="\n").writerows(rows)
        payload = stream.getvalue().encode("ascii")
        path.write_bytes(payload)
        authority["sha256"] = sha256(payload)
        authority["byte_count"] = len(payload)

    @classmethod
    def _set_csv_cell(
        cls,
        root: Path,
        document: dict[str, Any],
        artifact_index: int,
        row_index: int,
        column: str,
        value: str,
    ) -> None:
        def mutate(rows: list[list[str]]) -> None:
            column_index = rows[0].index(column)
            rows[row_index][column_index] = value

        cls._rewrite_csv(root, document, artifact_index, mutate)

    def _run(
        self,
        mutate_document: Callable[[dict[str, Any]], None] | None = None,
        mutate_workspace: Callable[[Path, dict[str, Any]], None] | None = None,
        raw_text: str | None = None,
        allow_synthetic: bool = True,
        bad_expected_digest: bool = False,
        audit: bool = False,
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, Any] | None]:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-fit-harness-") as temporary:
            root = Path(temporary)
            manifest = self._fixture(root)
            document = copy.deepcopy(BASE_DOCUMENT)
            if mutate_document is not None:
                mutate_document(document)
            if mutate_workspace is not None:
                mutate_workspace(root, document)
            payload = (
                raw_text.encode("ascii")
                if raw_text is not None
                else (json.dumps(document, indent=2) + "\n").encode("ascii")
            )
            manifest.write_bytes(payload)
            expected_digest = "0" * 64 if bad_expected_digest else sha256(payload)
            command = [
                sys.executable,
                str(TOOL),
                "--root",
                str(root),
                "--manifest",
                str(manifest),
                "--expected-manifest-sha256",
                expected_digest,
            ]
            if allow_synthetic:
                command.append("--allow-synthetic-test-fixture")
            audit_path = root / "audit.json"
            if audit:
                command.extend(("--audit-output", str(audit_path)))
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            audit_document = json.loads(audit_path.read_text(encoding="ascii")) if audit_path.exists() else None
            return result, audit_document

    def _require_failure(
        self, result: subprocess.CompletedProcess[str], reason: str
    ) -> None:
        self.__class__.checks += 1
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, output)
        self.assertIn(harness.FAIL_MARKER, result.stdout)
        self.assertIn(reason, result.stdout)
        self.assertNotIn(harness.PASS_MARKER, output)
        self.assertNotIn(harness.TEST_PASS_MARKER, output)
        self.assertNotIn("Traceback", output)

    def test_synthetic_fixture_passes_only_in_explicit_test_mode(self) -> None:
        result, audit = self._run(audit=True)
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(harness.TEST_PASS_MARKER, result.stdout)
        self.assertIn("fit_executed=0", result.stdout)
        self.assertIn("runtime_kij_emitted=0", result.stdout)
        self.assertIsNotNone(audit)
        assert audit is not None
        self.assertEqual(audit["status"], "synthetic-test-only-inputs-validated")
        self.assertFalse(audit["fit_executed"])
        self.assertFalse(audit["runtime_kij_emitted"])
        self.assertFalse(audit["binary_interaction_db_modified"])
        self.assertNotIn("parameter_values", audit)
        self.assertNotIn("kij", audit)

    def test_synthetic_fixture_is_rejected_by_default(self) -> None:
        result, _ = self._run(allow_synthetic=False)
        self._require_failure(result, "synthetic-fixture-requires-explicit-test-only-flag")

    def test_manifest_digest_is_an_external_required_binding(self) -> None:
        result, _ = self._run(bad_expected_digest=True)
        self._require_failure(result, "manifest-sha256-mismatch")

    def test_duplicate_root_key_is_rejected(self) -> None:
        needle = '  "data_class": "synthetic-nonphysical-test-only",'
        replacement = '  "data_class": "rights-cleared-experimental-vle-pvt",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        result, _ = self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1))
        self._require_failure(result, "duplicate-json-key:data_class")

    def test_duplicate_nested_rights_key_is_rejected(self) -> None:
        needle = '    "numerical_fit_use_permitted": true,'
        replacement = '    "numerical_fit_use_permitted": false,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        result, _ = self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1))
        self._require_failure(result, "duplicate-json-key:numerical_fit_use_permitted")

    def test_missing_rights_object_is_rejected(self) -> None:
        result, _ = self._run(lambda document: document.pop("rights"))
        self._require_failure(result, "missing-keys:root:rights")

    def test_rights_permission_is_required(self) -> None:
        result, _ = self._run(
            lambda document: document["rights"].__setitem__("numerical_fit_use_permitted", False)
        )
        self._require_failure(result, "value-mismatch:rights.numerical_fit_use_permitted")

    def test_rights_review_date_must_be_a_real_iso_date(self) -> None:
        result, _ = self._run(
            lambda document: document["rights"].__setitem__(
                "confirmed_on", "2026-99-99"
            )
        )
        self._require_failure(result, "invalid-date:rights.confirmed_on")

    def test_rights_workflow_must_select_exactly_one_mode(self) -> None:
        result, _ = self._run(
            lambda document: document["rights"].__setitem__("external_no_copy_workflow", True)
        )
        self._require_failure(result, "rights-workflow-must-select-exactly-one-mode")

    def test_reviewer_alias_only_physical_counterexample_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            document["data_class"] = harness.PHYSICAL_CLASS
            document["watermark"] = "RIGHTS-CLEARED EXPERIMENTAL INPUT CANDIDATE"
            document["rights"]["status"] = "resolved-for-declared-numerical-fit-use"
            document["rights"]["attribution"] = (
                "Unverified source with no provider, DOI, or URI in the schema."
            )
            document["rights"]["modification_notice"] = (
                "Training and holdout aliases originate from one campaign."
            )
            document["rights"]["confirmed_by"] = "self-asserted-unreviewed-identity"
            terms = root / document["rights"]["terms_artifact"]["path"]
            terms_payload = (
                b"USE PROHIBITED: no numerical fitting or redistribution permission.\n"
                b"Training and holdout are one source campaign under different aliases.\n"
            )
            terms.write_bytes(terms_payload)
            document["rights"]["terms_artifact"]["sha256"] = sha256(terms_payload)
            document["rights"]["terms_artifact"]["byte_count"] = len(terms_payload)

            aliases = (
                ("UNBOUND_TRAIN_ARTIFACT", "UNBOUND_TRAIN_DATASET", "SAME_CAMPAIGN_ALIAS_A"),
                ("UNBOUND_HOLDOUT_ARTIFACT", "UNBOUND_HOLDOUT_DATASET", "SAME_CAMPAIGN_ALIAS_B"),
            )
            for index, (artifact_id, dataset_id, campaign_id) in enumerate(aliases):
                authority = document["artifacts"][index]
                authority["artifact_id"] = artifact_id
                authority["dataset_id"] = dataset_id
                authority["campaign_id"] = campaign_id

                def rewrite(rows: list[list[str]]) -> None:
                    header = rows[0]
                    updates = {
                        "watermark": "RIGHTS_CLEARED_EXPERIMENTAL_DATA",
                        "artifact_id": artifact_id,
                        "dataset_id": dataset_id,
                        "campaign_id": campaign_id,
                    }
                    for row in rows[1:]:
                        for key, value in updates.items():
                            row[header.index(key)] = value

                self._rewrite_csv(root, document, index, rewrite)
            document["partition"]["training_dataset_ids"] = [aliases[0][1]]
            document["partition"]["holdout_dataset_ids"] = [aliases[1][1]]

        result, _ = self._run(mutate_workspace=mutate, allow_synthetic=False)
        self._require_failure(
            result,
            "missing-keys:artifacts[0]:provider,citation_or_dataset_doi,"
            "immutable_source_uri,retrieved_on,license_or_terms_url,"
            "license_or_terms_artifact,canonical_source_authority_sha256,"
            "canonical_campaign_id,canonical_campaign_authority_sha256,"
            "rights_disposition",
        )

    def test_canonical_authority_digest_ignores_free_form_aliases(self) -> None:
        source_fields = {
            "provider": "provider",
            "citation_or_dataset_doi": "doi:example",
            "immutable_source_uri": "https://example.invalid/archive.bin",
        }
        source_a = harness._canonical_authority_sha256(source_fields)
        source_b = harness._canonical_authority_sha256(dict(source_fields))
        campaign_fields = {
            "canonical_source_authority_sha256": source_a,
            "canonical_campaign_id": "campaign-from-provider",
        }
        self.__class__.checks += 1
        self.assertEqual(source_a, source_b)
        self.assertEqual(
            harness._canonical_authority_sha256(campaign_fields),
            harness._canonical_authority_sha256(dict(campaign_fields)),
        )

    def test_overlapping_canonical_authority_requires_reviewed_justification(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-overlap-review-") as temporary:
            with self.assertRaisesRegex(
                harness.KijFitHarnessError,
                "type-mismatch:partition.same_authority_nonoverlap_justification",
            ):
                harness._validate_nonoverlap_justification(
                    Path(temporary),
                    None,
                    ["training-dataset"],
                    ["holdout-dataset"],
                    ["0" * 64],
                    [],
                )
        self.__class__.checks += 1

    def test_artifact_rights_disposition_cannot_record_prohibited_use(self) -> None:
        artifact = {
            "artifact_id": "artifact",
            "sha256": "1" * 64,
            "provider": "provider",
            "citation_or_dataset_doi": "doi:example",
            "immutable_source_uri": "https://example.invalid/data",
            "license_or_terms_url": "https://example.invalid/terms",
        }
        disposition = {
            "schema": "rfam.kij-artifact-rights-disposition/v1",
            "contract_version": 1,
            "status": "use-prohibited",
            "artifact_id": "artifact",
            "artifact_sha256": "1" * 64,
            "provider": "provider",
            "citation_or_dataset_doi": "doi:example",
            "immutable_source_uri": "https://example.invalid/data",
            "license_or_terms_url": "https://example.invalid/terms",
            "license_or_terms_sha256": "2" * 64,
            "decision": {
                "numerical_fit_use_permitted": False,
                "repository_redistribution_permitted": False,
                "external_no_copy_workflow": True,
            },
            "reviewed_by": "reviewer",
            "reviewed_on": "2026-07-22",
        }
        payload = (json.dumps(disposition, indent=2) + "\n").encode("ascii")
        with self.assertRaisesRegex(
            harness.KijFitHarnessError,
            "value-mismatch:rights-disposition:artifact.status",
        ):
            harness._validate_rights_disposition(
                payload,
                artifact,
                "2" * 64,
                False,
                True,
            )
        self.__class__.checks += 1

    def test_terms_artifact_tamper_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            path = root / document["rights"]["terms_artifact"]["path"]
            path.write_bytes(path.read_bytes() + b"tamper\n")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "byte-count-mismatch:terms")

    def test_canonical_pair_is_required(self) -> None:
        result, _ = self._run(
            lambda document: document["pair"].__setitem__("canonical_pair_key", "O2|C3H8")
        )
        self._require_failure(result, "value-mismatch:pair.canonical_pair_key")

    def test_canonical_pressure_unit_is_required(self) -> None:
        result, _ = self._run(
            lambda document: document["units"].__setitem__("pressure", "kPa")
        )
        self._require_failure(result, "value-mismatch:units.pressure")

    def test_approval_integer_cannot_alias_false(self) -> None:
        result, _ = self._run(
            lambda document: document["approvals"].__setitem__("production_ready", 0)
        )
        self._require_failure(result, "value-mismatch:approvals.production_ready")

    def test_runtime_emission_cannot_be_enabled(self) -> None:
        result, _ = self._run(
            lambda document: document["fit_plan"].__setitem__("emit_runtime_kij", True)
        )
        self._require_failure(result, "value-mismatch:fit_plan.emit_runtime_kij")

    def test_binary_interaction_db_write_cannot_be_enabled(self) -> None:
        result, _ = self._run(
            lambda document: document["fit_plan"].__setitem__("write_binary_interaction_db", True)
        )
        self._require_failure(result, "value-mismatch:fit_plan.write_binary_interaction_db")

    def test_missing_phase_label_is_rejected_after_self_consistent_rehash(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(root, document, 0, 1, "measured_phase_label", "")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "missing-phase-label")

    def test_missing_liquid_tie_line_is_rejected_after_self_consistent_rehash(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(root, document, 0, 1, "liquid_x_o2", "")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "missing-or-malformed-number")

    def test_missing_vapor_tie_line_is_rejected_after_self_consistent_rehash(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(root, document, 0, 1, "vapor_y_o2", "")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "missing-or-malformed-number")

    def test_missing_uncertainty_is_rejected_after_self_consistent_rehash(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(
                root, document, 0, 1, "vapor_y_o2_standard_uncertainty", ""
            )

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "missing-or-malformed-number")

    def test_zero_uncertainty_is_rejected_after_self_consistent_rehash(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(
                root, document, 0, 1, "pressure_standard_uncertainty_pa", "0"
            )

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "nonpositive-uncertainty")

    def test_training_holdout_measurement_overlap_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            training_path = root / document["artifacts"][0]["path"]
            training_rows = list(csv.reader(io.StringIO(training_path.read_text(encoding="ascii"))))
            header = training_rows[0]
            source = dict(zip(header, training_rows[1]))

            def replace_measurement(rows: list[list[str]]) -> None:
                for field in (
                    "measurement_class",
                    "temperature_k",
                    "pressure_pa",
                    "measured_phase_label",
                    "liquid_x_o2",
                    "vapor_y_o2",
                ):
                    rows[1][header.index(field)] = source[field]

            self._rewrite_csv(root, document, 1, replace_measurement)

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "training-holdout-measurement-overlap")

    def test_training_holdout_campaign_overlap_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            campaign = document["artifacts"][0]["campaign_id"]
            document["artifacts"][1]["campaign_id"] = campaign
            self._set_csv_cell(root, document, 1, 1, "campaign_id", campaign)

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "training-holdout-campaign-overlap")

    def test_duplicate_record_identity_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            def duplicate(rows: list[list[str]]) -> None:
                rows.append(rows[1].copy())

            self._rewrite_csv(root, document, 0, duplicate)

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "duplicate-record-identity")

    def test_record_outside_temperature_envelope_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(root, document, 0, 1, "temperature_k", "150")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "temperature-outside-envelope")

    def test_invalid_pressure_envelope_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["envelope"]["pressure_pa"]["minimum"] = -1.0

        result, _ = self._run(mutate)
        self._require_failure(result, "nonpositive-range:envelope.pressure_pa")

    def test_artifact_payload_tamper_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            path = root / document["artifacts"][0]["path"]
            path.write_bytes(path.read_bytes() + b"tamper\n")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "byte-count-mismatch:artifact")

    def test_artifact_path_escape_is_rejected(self) -> None:
        result, _ = self._run(
            lambda document: document["artifacts"][0].__setitem__("path", "../outside.csv")
        )
        self._require_failure(result, "path-escapes-root")

    def test_unknown_nested_key_is_rejected(self) -> None:
        result, _ = self._run(
            lambda document: document["rights"].__setitem__("unreviewed", True)
        )
        self._require_failure(result, "unknown-keys:rights:unreviewed")

    def test_synthetic_row_watermark_cannot_be_erased(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            self._set_csv_cell(root, document, 0, 1, "watermark", "")

        result, _ = self._run(mutate_workspace=mutate)
        self._require_failure(result, "value-mismatch:artifact")

    def test_harness_does_not_modify_binary_interaction_database(self) -> None:
        before = sha256(DB_SOURCE.read_bytes())
        result, _ = self._run(audit=True)
        after = sha256(DB_SOURCE.read_bytes())
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(before, after)

    def test_audit_output_cannot_target_binary_interaction_database(self) -> None:
        before = sha256(DB_SOURCE.read_bytes())
        digest = sha256(MANIFEST.read_bytes())
        result = subprocess.run(
            [
                sys.executable,
                str(TOOL),
                "--root",
                str(WORKSPACE),
                "--manifest",
                str(MANIFEST),
                "--expected-manifest-sha256",
                digest,
                "--allow-synthetic-test-fixture",
                "--audit-output",
                str(DB_SOURCE),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        after = sha256(DB_SOURCE.read_bytes())
        self._require_failure(result, "audit-output-targets-binary-interaction-db")
        self.assertEqual(before, after)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(KijFitHarnessTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"RFAM_KIJ_FIT_HARNESS_TEST_PASS checks={KijFitHarnessTests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
