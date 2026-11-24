#!/usr/bin/env python3
"""Focused fail-closed tests for the bounded SciExpeM R2 intake."""

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
TOOL = REPO / "tools" / "validate_r2_experimental_intake.py"
MANIFEST = REPO / "tests" / "reference" / "r2_experimental_sciexpem_intake_v1.json"
SOURCE_DIR = REPO / "tests" / "reference" / "r2_experimental_sciexpem_v1"
SOURCE_RELATIVE_DIR = Path("realFluidAllMach/tests/reference/r2_experimental_sciexpem_v1")
BASE_TEXT = MANIFEST.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)

SPEC = importlib.util.spec_from_file_location("r2_experimental_intake", TOOL)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import R2 experimental intake validator")
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class R2ExperimentalIntakeTests(unittest.TestCase):
    checks = 0

    def _run(
        self,
        mutate: Optional[Callable[[Dict[str, Any]], None]] = None,
        raw_text: Optional[str] = None,
        mutate_root: Optional[Callable[[Path], None]] = None,
        output: bool = False,
        overwrite: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="r2-experimental-intake-") as temporary:
            temp = Path(temporary)
            root = temp / "workspace"
            target_sources = root / SOURCE_RELATIVE_DIR
            target_sources.parent.mkdir(parents=True)
            shutil.copytree(SOURCE_DIR, target_sources)
            manifest_path = temp / "manifest.json"
            if raw_text is None:
                document = json.loads(BASE_TEXT)
                if mutate is not None:
                    mutate(document)
                raw_text = json.dumps(document, indent=2, ensure_ascii=True) + "\n"
            manifest_path.write_bytes(raw_text.encode("ascii"))
            if mutate_root is not None:
                mutate_root(root)
            command = [
                sys.executable,
                str(TOOL),
                "--root",
                str(root),
                "--manifest",
                str(manifest_path),
            ]
            output_path = temp / "result.json"
            if output:
                command.extend(["--output", str(output_path)])
            if overwrite:
                command.append("--overwrite")
            result = subprocess.run(command, text=True, capture_output=True, check=False)
            if output_path.exists():
                result.output_payload = output_path.read_bytes()  # type: ignore[attr-defined]
            return result

    def _require_failure(self, result: subprocess.CompletedProcess[str], reason: str) -> None:
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertIn('"status":"FAIL"', result.stderr)
        self.assertIn(reason, result.stderr)
        self.assertNotIn(validator.BLOCKED_MARKER, result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_frozen_intake_is_consistent_but_not_fit_eligible(self) -> None:
        result = self._run()
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stderr, "")
        self.assertIn(validator.BLOCKED_MARKER, result.stdout)
        self.assertIn("records=6", result.stdout)
        self.assertIn("points=48", result.stdout)
        self.assertIn("c3h8_curves=0", result.stdout)
        self.assertIn("training_eligible=0", result.stdout)
        self.assertIn("holdout_eligible=0", result.stdout)

    def test_canonical_output_has_frozen_identity_and_scientific_disposition(self) -> None:
        result = self._run(output=True)
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        payload = result.output_payload  # type: ignore[attr-defined]
        self.assertEqual(len(payload), validator.EXPECTED_EXTRACTION["expected_byte_count"])
        self.assertEqual(
            hashlib.sha256(payload).hexdigest(),
            validator.EXPECTED_EXTRACTION["expected_sha256"],
        )
        document = json.loads(payload)
        self.assertEqual(document["summary"]["record_count"], 6)
        self.assertEqual(document["summary"]["point_count"], 48)
        self.assertEqual(document["summary"]["species_set"], ["AR", "H2", "O2"])
        self.assertEqual(document["summary"]["c3h8_curve_count"], 0)
        self.assertFalse(document["eligibility"]["training_eligible"])
        self.assertFalse(document["eligibility"]["holdout_eligible"])
        self.assertTrue(all(not item["contains_c3h8"] for item in document["datasets"]))

    def test_output_is_byte_deterministic(self) -> None:
        first = self._run(output=True)
        second = self._run(output=True)
        self.__class__.checks += 1
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertEqual(
            first.output_payload,  # type: ignore[attr-defined]
            second.output_payload,  # type: ignore[attr-defined]
        )

    def test_last_record_property_id_order_is_not_misread(self) -> None:
        result = self._run(output=True)
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        document = json.loads(result.output_payload)  # type: ignore[attr-defined]
        last = document["datasets"][-1]
        self.assertEqual(last["record_id"], 4639759)
        self.assertEqual(last["points"][0]["temperature_k"], "1176.0000000000")
        self.assertEqual(last["points"][0]["ignition_delay_us"], "36.0000000000")

    def test_duplicate_root_key_is_rejected(self) -> None:
        needle = '  "status": "blocked-no-propane-curves",'
        replacement = '  "status": "fit-ready",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:status",
        )

    def test_duplicate_nested_approval_key_is_rejected(self) -> None:
        needle = '    "production_ready": false'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:production_ready",
        )

    def test_duplicate_nested_source_key_is_rejected(self) -> None:
        needle = '      "record_id": 4639716,'
        replacement = '      "record_id": 1,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:record_id",
        )

    def test_unknown_root_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("unreviewed", True)),
            "contract-unknown-keys:root:unreviewed",
        )

    def test_unknown_nested_source_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["sources"][0].__setitem__(
                    "alternate_authority", True
                )
            ),
            "contract-unknown-keys:root.sources[0]:alternate_authority",
        )

    def test_contract_version_bool_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("contract_version", True)),
            "contract-type-mismatch:root.contract_version",
        )

    def test_noncanonical_manifest_is_rejected(self) -> None:
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace("{\n", "{ \n", 1)),
            "manifest-not-canonical-json",
        )

    def test_source_url_hash_and_bytes_cannot_be_self_consistently_rebound(self) -> None:
        payload = b"x" * 2328

        def mutate(document: Dict[str, Any]) -> None:
            source = document["sources"][0]
            source["source_url"] = "https://example.invalid/substitute.xml"
            source["byte_count"] = len(payload)
            source["sha256"] = hashlib.sha256(payload).hexdigest()
            source["zenodo_md5"] = hashlib.md5(payload).hexdigest()

        def mutate_root(root: Path) -> None:
            path = root / Path(BASE_DOCUMENT["sources"][0]["repository_path"])
            path.write_bytes(payload)

        self._require_failure(
            self._run(mutate=mutate, mutate_root=mutate_root),
            "contract-value-mismatch:root.sources[0].source_url",
        )

    def test_source_path_hash_and_bytes_cannot_be_self_consistently_rebound(self) -> None:
        payload = b"substitute-source\n"
        alternate = "realFluidAllMach/tests/reference/r2_experimental_sciexpem_v1/substitute.xml"

        def mutate(document: Dict[str, Any]) -> None:
            source = document["sources"][0]
            source["repository_path"] = alternate
            source["byte_count"] = len(payload)
            source["sha256"] = hashlib.sha256(payload).hexdigest()
            source["zenodo_md5"] = hashlib.md5(payload).hexdigest()

        def mutate_root(root: Path) -> None:
            path = root / Path(alternate)
            path.write_bytes(payload)

        self._require_failure(
            self._run(mutate=mutate, mutate_root=mutate_root),
            "contract-value-mismatch:root.sources[0].repository_path",
        )

    def test_source_entries_cannot_be_semantically_swapped(self) -> None:
        def mutate(document: Dict[str, Any]) -> None:
            document["sources"][0], document["sources"][1] = (
                document["sources"][1],
                document["sources"][0],
            )

        self._require_failure(
            self._run(mutate),
            "contract-value-mismatch:root.sources[0].record_id",
        )

    def test_live_source_tamper_is_rejected_before_xml_parse(self) -> None:
        def mutate_root(root: Path) -> None:
            path = root / Path(BASE_DOCUMENT["sources"][0]["repository_path"])
            payload = bytearray(path.read_bytes())
            payload[0] = ord("X")
            path.write_bytes(payload)

        self._require_failure(
            self._run(mutate_root=mutate_root),
            "source-sha256-mismatch:4639716",
        )

    def test_missing_source_is_rejected(self) -> None:
        def mutate_root(root: Path) -> None:
            path = root / Path(BASE_DOCUMENT["sources"][3]["repository_path"])
            path.unlink()

        self._require_failure(
            self._run(mutate_root=mutate_root),
            "source-stat-failed:4639741",
        )

    def test_license_cannot_be_reclassified(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["license"].__setitem__(
                    "record_license_id", "cc0-1.0"
                )
            ),
            "contract-value-mismatch:root.license.record_license_id",
        )

    def test_digitized_provenance_cannot_be_erased(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["observed_structure"].__setitem__(
                    "ignition_delay_source_type", "reported"
                )
            ),
            "contract-value-mismatch:root.observed_structure.ignition_delay_source_type",
        )

    def test_uncertainty_cannot_be_claimed_present(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["uncertainty_availability"].__setitem__(
                    "temperature", True
                )
            ),
            "contract-value-mismatch:root.uncertainty_availability.temperature",
        )

    def test_training_and_production_approvals_cannot_be_enabled(self) -> None:
        for key in ("experimental_training_data_selected", "production_ready"):
            with self.subTest(key=key):
                self._require_failure(
                    self._run(
                        lambda document, selected=key: document["approvals"].__setitem__(
                            selected, True
                        )
                    ),
                    "contract-value-mismatch:root.approvals.%s" % key,
                )

    def test_blocking_reason_cannot_be_removed(self) -> None:
        self._require_failure(
            self._run(lambda document: document["eligibility"]["blocking_reasons"].pop()),
            "contract-list-length-mismatch:root.eligibility.blocking_reasons",
        )

    def test_extraction_identity_cannot_be_self_asserted(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["extraction"].__setitem__(
                    "expected_sha256", "0" * 64
                )
            ),
            "contract-value-mismatch:root.extraction.expected_sha256",
        )

    def test_overwrite_requires_output(self) -> None:
        self._require_failure(self._run(overwrite=True), "overwrite-requires-output")


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(R2ExperimentalIntakeTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if outcome.wasSuccessful():
        print(
            "R2_EXPERIMENTAL_INTAKE_TEST_PASS checks=%d"
            % R2ExperimentalIntakeTests.checks
        )
    raise SystemExit(0 if outcome.wasSuccessful() else 1)
