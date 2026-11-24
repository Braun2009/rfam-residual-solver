#!/usr/bin/env python3
"""Focused tests for the internal-only Lam 2010 propane IDT intake."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable, Dict, Optional


REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
VALIDATOR = REPO / "tools" / "validate_r2_lam2010_intake.py"
GENERATOR = REPO / "tools" / "generate_r2_lam2010_intake.py"
MANIFEST = REPO / "tests" / "reference" / "r2_lam2010_proci_intake_v1.json"
CSV = REPO / "tests" / "reference" / "r2_lam2010_proci_intake_v1.csv"
SOURCE = WORKSPACE / "tmp" / "lam-2010-proci-supplement-mmc1.doc"
RIGHTS_EVIDENCE = WORKSPACE / "tmp" / "crossref-lam-2010.json"
BASE_TEXT = MANIFEST.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)
BASE_CSV = CSV.read_bytes()
BASE_RIGHTS_EVIDENCE = RIGHTS_EVIDENCE.read_bytes()

SPEC = importlib.util.spec_from_file_location("r2_lam2010_intake", VALIDATOR)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import Lam intake validator")
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class Lam2010IntakeTests(unittest.TestCase):
    checks = 0

    def _run(
        self,
        mutate: Optional[Callable[[Dict[str, Any]], None]] = None,
        raw_text: Optional[str] = None,
        mutate_root: Optional[Callable[[Path], None]] = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="r2-lam2010-intake-") as temporary:
            temp = Path(temporary)
            root = temp / "workspace"
            source_target = root / Path(validator.SOURCE_PATH)
            rights_evidence_target = root / Path(validator.RIGHTS_EVIDENCE_PATH)
            generator_target = root / Path(validator.GENERATOR_PATH)
            csv_target = root / Path(validator.CSV_PATH)
            source_target.parent.mkdir(parents=True)
            rights_evidence_target.parent.mkdir(parents=True, exist_ok=True)
            generator_target.parent.mkdir(parents=True)
            csv_target.parent.mkdir(parents=True)
            shutil.copyfile(SOURCE, source_target)
            shutil.copyfile(RIGHTS_EVIDENCE, rights_evidence_target)
            shutil.copyfile(GENERATOR, generator_target)
            shutil.copyfile(CSV, csv_target)

            if raw_text is None:
                document = json.loads(BASE_TEXT)
                if mutate is not None:
                    mutate(document)
                raw_text = json.dumps(document, indent=2, ensure_ascii=True) + "\n"
            manifest_path = temp / "manifest.json"
            manifest_path.write_bytes(raw_text.encode("ascii"))
            if mutate_root is not None:
                mutate_root(root)
            return subprocess.run(
                [
                    sys.executable,
                    str(VALIDATOR),
                    "--root",
                    str(root),
                    "--manifest",
                    str(manifest_path),
                ],
                text=True,
                capture_output=True,
                check=False,
            )

    def _require_failure(self, result: subprocess.CompletedProcess[str], reason: str) -> None:
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertIn('"status":"FAIL"', result.stderr)
        self.assertIn(reason, result.stderr)
        self.assertNotIn(validator.BLOCKED_MARKER, result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_frozen_intake_is_consistent_but_strictly_blocked(self) -> None:
        result = self._run(raw_text=BASE_TEXT)
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stderr, "")
        self.assertIn(validator.BLOCKED_MARKER, result.stdout)
        self.assertIn("tables=6", result.stdout)
        self.assertIn("rows=95", result.stdout)
        self.assertIn("internal_screening_only=1", result.stdout)
        self.assertIn("fit_eligible=0", result.stdout)
        self.assertIn("redistribution=0", result.stdout)

    def test_csv_identity_and_table_semantics(self) -> None:
        rows = validator._parse_csv(BASE_CSV, validator.EXPECTED_TABLES)
        self.__class__.checks += 1
        self.assertEqual(len(BASE_CSV), validator.CSV_BYTE_COUNT)
        self.assertEqual(hashlib.sha256(BASE_CSV).hexdigest(), validator.CSV_SHA256)
        self.assertEqual(len(rows), 95)
        self.assertEqual([sum(row["table_id"] == "S%d" % index for row in rows) for index in range(1, 7)], [17, 18, 19, 18, 10, 13])
        self.assertEqual(sum(row["facility_mode"] == "conventional-pressure-rise" for row in rows), 46)
        self.assertEqual(sum(row["facility_mode"] == "near-constant-volume" for row in rows), 49)
        self.assertEqual(sum(row["dp5_dt_status"] == "reported" for row in rows), 23)
        self.assertEqual(sum(row["dp5_dt_status"] == "not-available-cell" for row in rows), 23)
        self.assertEqual(sum(row["dp5_dt_status"] == "not-tabulated-caption-approximately-zero" for row in rows), 49)

    def test_rights_boundary_is_exact(self) -> None:
        rights = BASE_DOCUMENT["rights"]
        self.__class__.checks += 1
        self.assertEqual(rights["status"], "unresolved")
        self.assertIs(rights["internal_screening_only"], True)
        self.assertIs(rights["fit_eligible"], False)
        self.assertIs(rights["redistribution"], False)
        self.assertIs(rights["commercial_or_production_use_authorized"], False)
        self.assertIs(rights["crossref_metadata_is_authorization"], False)
        self.assertFalse(any(BASE_DOCUMENT["approvals"].values()))

    def test_reported_composition_is_not_normalized_or_reclassified(self) -> None:
        experiment = BASE_DOCUMENT["experiment"]
        self.__class__.checks += 1
        self.assertEqual(experiment["reported_initial_mixture"], "0.8% C3H8/ 8% O2/ Ar")
        self.assertEqual(experiment["reported_composition_basis"], "not-stated-in-supplement")
        self.assertEqual(experiment["reported_species_values"]["Ar"], "balance-unspecified")
        self.assertIs(experiment["normalized_composition_available"], False)
        self.assertIs(experiment["equivalence_ratio_reported"], False)
        self.assertNotIn("initial_composition", experiment)

    def test_duplicate_root_key_is_rejected(self) -> None:
        needle = '  "status": "blocked-internal-screening-only",'
        replacement = '  "status": "fit-ready",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:status",
        )

    def test_duplicate_source_key_is_rejected(self) -> None:
        needle = '    "sha256": "154b690336f9bd2aa6a56c417386d620d5dc360e5e996a8b41900238bd17f977",'
        replacement = '    "sha256": "' + "0" * 64 + '",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:sha256",
        )

    def test_escaped_duplicate_approval_key_is_rejected(self) -> None:
        needle = '    "production_ready": false'
        replacement = '    "production_\\u0072eady": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:production_ready",
        )

    def test_unknown_root_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("unreviewed", True)),
            "contract-unknown-keys:root:unreviewed",
        )

    def test_unknown_source_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document["source"].__setitem__("alternate", True)),
            "contract-unknown-keys:root.source:alternate",
        )

    def test_unknown_table_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document["experiment"]["tables"][0].__setitem__("override", True)),
            "contract-unknown-keys:root.experiment.tables[0]:override",
        )

    def test_bool_contract_version_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("contract_version", True)),
            "contract-type-mismatch:root.contract_version",
        )

    def test_manifest_byte_drift_is_rejected(self) -> None:
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace("{\n", "{ \n", 1)),
            "manifest-sha256-mismatch",
        )

    def test_source_authority_cannot_be_self_consistently_rebound(self) -> None:
        substitute = b"x" * validator.SOURCE_BYTE_COUNT
        alternate = "tmp/substitute.doc"

        def mutate(document: Dict[str, Any]) -> None:
            source = document["source"]
            source["url"] = "https://example.invalid/substitute.doc"
            source["local_path"] = alternate
            source["sha256"] = hashlib.sha256(substitute).hexdigest()

        def mutate_root(root: Path) -> None:
            path = root / Path(alternate)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(substitute)

        self._require_failure(
            self._run(mutate=mutate, mutate_root=mutate_root),
            "contract-value-mismatch:root.source.url",
        )

    def test_live_source_tamper_is_rejected(self) -> None:
        def mutate_root(root: Path) -> None:
            path = root / Path(validator.SOURCE_PATH)
            payload = bytearray(path.read_bytes())
            payload[100] ^= 1
            path.write_bytes(payload)

        self._require_failure(
            self._run(raw_text=BASE_TEXT, mutate_root=mutate_root),
            "file-sha256-mismatch:source",
        )

    def test_rights_evidence_authority_cannot_be_self_consistently_rebound(self) -> None:
        substitute = b'{"status":"ok","message":{}}\n'
        alternate = "tmp/substitute-crossref.json"

        def mutate(document: Dict[str, Any]) -> None:
            evidence = document["rights_evidence"]
            evidence["url"] = "https://example.invalid/rights.json"
            evidence["local_path"] = alternate
            evidence["byte_count"] = len(substitute)
            evidence["sha256"] = hashlib.sha256(substitute).hexdigest()

        def mutate_root(root: Path) -> None:
            path = root / Path(alternate)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(substitute)

        self._require_failure(
            self._run(mutate=mutate, mutate_root=mutate_root),
            "contract-value-mismatch:root.rights_evidence.url",
        )

    def test_live_rights_evidence_tamper_is_rejected(self) -> None:
        def mutate_root(root: Path) -> None:
            path = root / Path(validator.RIGHTS_EVIDENCE_PATH)
            path.write_bytes(path.read_bytes() + b"\n")

        self._require_failure(
            self._run(raw_text=BASE_TEXT, mutate_root=mutate_root),
            "file-byte-count-mismatch:rights-evidence",
        )

    def test_rights_evidence_semantics_reject_open_license_substitution(self) -> None:
        document = json.loads(BASE_RIGHTS_EVIDENCE.decode("utf-8"))
        document["message"]["license"][0]["content-version"] = "vor"
        document["message"]["license"][0]["URL"] = "https://creativecommons.org/licenses/by/4.0/"
        payload = json.dumps(document, ensure_ascii=True, separators=(",", ":")).encode("ascii")
        self.__class__.checks += 1
        with self.assertRaisesRegex(
            validator.IntakeError,
            r"contract-value-mismatch:rights_evidence.message.license\[0\].content-version",
        ):
            validator._validate_rights_evidence(payload)

    def test_generator_authority_cannot_be_self_consistently_rebound(self) -> None:
        substitute = b"#!/usr/bin/env python3\nprint('fake')\n"
        alternate = "realFluidAllMach/tools/fake_lam_generator.py"

        def mutate(document: Dict[str, Any]) -> None:
            extraction = document["extraction"]
            extraction["generator_path"] = alternate
            extraction["generator_sha256"] = hashlib.sha256(substitute).hexdigest()

        def mutate_root(root: Path) -> None:
            path = root / Path(alternate)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(substitute)

        self._require_failure(
            self._run(mutate=mutate, mutate_root=mutate_root),
            "contract-value-mismatch:root.extraction.generator_path",
        )

    def test_live_generator_tamper_is_rejected(self) -> None:
        def mutate_root(root: Path) -> None:
            path = root / Path(validator.GENERATOR_PATH)
            path.write_bytes(path.read_bytes() + b"\n")

        self._require_failure(
            self._run(raw_text=BASE_TEXT, mutate_root=mutate_root),
            "file-sha256-mismatch:generator",
        )

    def test_csv_authority_cannot_be_self_consistently_rebound(self) -> None:
        substitute = BASE_CSV.replace(b"1405", b"1406", 1)
        alternate = "realFluidAllMach/tests/reference/substitute_lam.csv"

        def mutate(document: Dict[str, Any]) -> None:
            extraction = document["extraction"]
            extraction["csv_path"] = alternate
            extraction["csv_byte_count"] = len(substitute)
            extraction["csv_sha256"] = hashlib.sha256(substitute).hexdigest()

        def mutate_root(root: Path) -> None:
            path = root / Path(alternate)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(substitute)

        self._require_failure(
            self._run(mutate=mutate, mutate_root=mutate_root),
            "contract-value-mismatch:root.extraction.csv_path",
        )

    def test_live_csv_tamper_is_rejected(self) -> None:
        def mutate_root(root: Path) -> None:
            path = root / Path(validator.CSV_PATH)
            path.write_bytes(path.read_bytes().replace(b"1405", b"1406", 1))

        self._require_failure(
            self._run(raw_text=BASE_TEXT, mutate_root=mutate_root),
            "file-sha256-mismatch:csv",
        )

    def test_internal_screening_only_cannot_be_disabled(self) -> None:
        self._require_failure(
            self._run(lambda document: document["rights"].__setitem__("internal_screening_only", False)),
            "contract-value-mismatch:root.rights.internal_screening_only",
        )

    def test_fit_or_redistribution_cannot_be_enabled(self) -> None:
        for key in ("fit_eligible", "redistribution"):
            with self.subTest(key=key):
                self._require_failure(
                    self._run(lambda document, selected=key: document["rights"].__setitem__(selected, True)),
                    "contract-value-mismatch:root.rights.%s" % key,
                )

    def test_approval_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(lambda document: document["approvals"].__setitem__("r2_physical_fit_approved", True)),
            "contract-value-mismatch:root.approvals.r2_physical_fit_approved",
        )

    def test_table_count_and_facility_mode_cannot_drift(self) -> None:
        self._require_failure(
            self._run(lambda document: document["experiment"].__setitem__("table_count", 5)),
            "contract-value-mismatch:root.experiment.table_count",
        )
        self._require_failure(
            self._run(lambda document: document["experiment"]["tables"][0].__setitem__("facility_mode", "near-constant-volume")),
            "contract-value-mismatch:root.experiment.tables[0].facility_mode",
        )

    def test_parser_rejects_near_volume_zero_imputation(self) -> None:
        mutated = BASE_CSV.replace(
            b"S2,1,near-constant-volume,1381,5.29,,not-tabulated-caption-approximately-zero,0.15",
            b"S2,1,near-constant-volume,1381,5.29,0,not-tabulated-caption-approximately-zero,0.15",
            1,
        )
        self.__class__.checks += 1
        with self.assertRaisesRegex(validator.IntakeError, "csv-near-volume-dpdt-mismatch:S2-R1"):
            validator._parse_csv(mutated, validator.EXPECTED_TABLES)

    def test_parser_rejects_na_imputation(self) -> None:
        mutated = BASE_CSV.replace(
            b"S1,1,conventional-pressure-rise,1405,5.46,,not-available-cell,0.083",
            b"S1,1,conventional-pressure-rise,1405,5.46,0,not-available-cell,0.083",
            1,
        )
        self.__class__.checks += 1
        with self.assertRaisesRegex(validator.IntakeError, "csv-missing-dpdt-has-value:S1-R1"):
            validator._parse_csv(mutated, validator.EXPECTED_TABLES)

    def test_parser_rejects_row_and_field_drift(self) -> None:
        mutations = [
            (BASE_CSV.replace(b"S1,2,", b"S1,3,", 1), "csv-row-index-mismatch:S1-R2"),
            (BASE_CSV.replace(b"conventional-pressure-rise", b"near-constant-volume", 1), "csv-facility-mode-mismatch:S1-R1"),
            (BASE_CSV.replace(b"t5_k", b"temperature_k", 1), "csv-header-mismatch"),
        ]
        for payload, reason in mutations:
            with self.subTest(reason=reason):
                self.__class__.checks += 1
                with self.assertRaisesRegex(validator.IntakeError, reason):
                    validator._parse_csv(payload, validator.EXPECTED_TABLES)

    def test_generator_rejects_wrong_source_before_word(self) -> None:
        with tempfile.TemporaryDirectory(prefix="r2-lam-generator-negative-") as temporary:
            temp = Path(temporary)
            source = temp / "bad.doc"
            source.write_bytes(b"not-the-frozen-source")
            result = subprocess.run(
                [
                    sys.executable,
                    str(GENERATOR),
                    "--source",
                    str(source),
                    "--output",
                    str(temp / "out.csv"),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertIn("source-byte-count-mismatch", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Lam2010IntakeTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if outcome.wasSuccessful():
        print("R2_LAM2010_INTAKE_TEST_PASS checks=%d" % Lam2010IntakeTests.checks)
    raise SystemExit(0 if outcome.wasSuccessful() else 1)
