#!/usr/bin/env python3
"""Focused fail-close tests for the WP3 O2/C3H8 source candidate gate."""

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
TOOL = REPO / "tools/validate_kij_source_candidate.py"
MANIFEST = REPO / "tests/reference/kij_o2_c3h8_source_candidate_v1.json"
BASE_TEXT = MANIFEST.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)

sys.path.insert(0, str(REPO / "tools"))
import validate_kij_source_candidate as validator  # noqa: E402


class KijSourceCandidateTests(unittest.TestCase):
    checks = 0

    def _fixture(self, root: Path) -> Path:
        for item in validator.EXPECTED_UPSTREAM:
            source = WORKSPACE / item["path"]
            target = root / item["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        target = (
            root
            / "realFluidAllMach/tests/reference/kij_o2_c3h8_source_candidate_v1.json"
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        return target

    def _run(
        self,
        mutate: Callable[[dict[str, Any]], None] | None = None,
        raw_text: str | None = None,
        source_payload: bytes | None = None,
        mutate_workspace: Callable[[Path], None] | None = None,
        extra_args: list[str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-source-") as temporary:
            root = Path(temporary)
            manifest = self._fixture(root)
            document = copy.deepcopy(BASE_DOCUMENT)
            if mutate is not None:
                mutate(document)
            with manifest.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(
                    raw_text
                    if raw_text is not None
                    else json.dumps(document, indent=2) + "\n"
                )
            if mutate_workspace is not None:
                mutate_workspace(root)
            command = [
                sys.executable,
                str(TOOL),
                "--root",
                str(root),
                "--manifest",
                str(manifest),
            ]
            if source_payload is not None:
                source = root / "source.xml"
                source.write_bytes(source_payload)
                command.extend(["--source-xml", str(source)])
            if extra_args:
                command.extend(extra_args)
            return subprocess.run(command, capture_output=True, text=True, check=False)

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

    def test_frozen_candidate_is_consistent_but_blocked(self) -> None:
        result = self._run()
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(validator.BLOCKED_MARKER, result.stdout)
        self.assertIn("live_source_verified=0", result.stdout)
        self.assertIn("source_selected=0", result.stdout)
        self.assertIn("fit_approved=0", result.stdout)
        self.assertIn("production_ready=0", result.stdout)

    def test_duplicate_approval_key_is_rejected(self) -> None:
        needle = '    "fit_approved": false,'
        replacement = '    "fit_approved": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:fit_approved",
        )

    def test_duplicate_nested_rights_key_is_rejected(self) -> None:
        needle = (
            '    "classification": '
            '"unresolved-artifact-specific-redistribution-and-fit-use",'
        )
        replacement = '    "classification": "resolved",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:classification",
        )

    def test_noncanonical_json_is_rejected(self) -> None:
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace("{\n", "{ \n", 1)),
            "manifest-not-canonical-json",
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

    def test_fit_approval_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["approvals"].__setitem__(
                    "fit_approved", True
                )
            ),
            "contract-value-mismatch:root.approvals.fit_approved",
        )

    def test_source_selection_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["approvals"].__setitem__(
                    "source_selected", True
                )
            ),
            "contract-value-mismatch:root.approvals.source_selected",
        )

    def test_functional_form_cannot_be_silently_selected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["pair"].__setitem__(
                    "fit_functional_form_selected", True
                )
            ),
            "contract-value-mismatch:root.pair.fit_functional_form_selected",
        )

    def test_upstream_path_and_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first, second = document["upstream"][0], document["upstream"][1]
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

    def test_source_url_change_is_rejected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["source_artifact"].__setitem__(
                    "url", "https://example.invalid/source.xml"
                )
            ),
            "contract-value-mismatch:root.source_artifact.url",
        )

    def test_source_hash_and_matching_fake_bytes_are_rejected_preparse(self) -> None:
        payload = b"x" * validator.SOURCE_BYTE_COUNT

        def mutate(document: dict[str, Any]) -> None:
            import hashlib

            document["source_artifact"]["sha256"] = hashlib.sha256(payload).hexdigest()

        self._require_failure(
            self._run(mutate=mutate, source_payload=payload),
            "contract-value-mismatch:root.source_artifact.sha256",
        )

    def test_same_length_source_tamper_is_rejected_before_xml_parse(self) -> None:
        self._require_failure(
            self._run(source_payload=b"x" * validator.SOURCE_BYTE_COUNT),
            "source-sha256-mismatch",
        )

    def test_source_byte_count_change_is_rejected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["source_artifact"].__setitem__(
                    "byte_count", validator.SOURCE_BYTE_COUNT + 1
                )
            ),
            "contract-value-mismatch:root.source_artifact.byte_count",
        )

    def test_repository_copy_cannot_be_claimed(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["source_artifact"].__setitem__(
                    "repository_copy", True
                )
            ),
            "contract-value-mismatch:root.source_artifact.repository_copy",
        )

    def test_rights_cannot_be_marked_resolved(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["rights"].__setitem__(
                    "classification", "resolved"
                )
            ),
            "contract-value-mismatch:root.rights.classification",
        )

    def test_rights_evidence_url_cannot_be_redirected(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["rights"]["evidence"][1].__setitem__(
                    "url", "https://example.invalid/terms"
                )
            ),
            "contract-value-mismatch:root.rights.evidence[1].url",
        )

    def test_missing_temperature_uncertainty_cannot_be_erased(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["experiment"]["missing_measurement_fields"].remove(
                "temperature_uncertainty_k"
            )

        self._require_failure(
            self._run(mutate),
            "contract-list-length-mismatch:root.experiment.missing_measurement_fields",
        )

    def test_vapor_composition_cannot_be_claimed_present(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["experiment"]["missing_measurement_fields"].remove("vapor_y_o2")

        self._require_failure(
            self._run(mutate),
            "contract-list-length-mismatch:root.experiment.missing_measurement_fields",
        )

    def test_candidate_csv_cannot_be_reclassified_as_fit_input(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["candidate_extraction"].__setitem__(
                    "fit_input", True
                )
            ),
            "contract-value-mismatch:root.candidate_extraction.fit_input",
        )

    def test_candidate_hash_cannot_be_self_asserted(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["candidate_extraction"].__setitem__(
                    "sha256", "a" * 64
                )
            ),
            "contract-value-mismatch:root.candidate_extraction.sha256",
        )

    def test_training_eligibility_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["eligibility"].__setitem__(
                    "training_eligible", True
                )
            ),
            "contract-value-mismatch:root.eligibility.training_eligible",
        )

    def test_holdout_eligibility_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(
                lambda document: document["eligibility"].__setitem__(
                    "holdout_eligible", True
                )
            ),
            "contract-value-mismatch:root.eligibility.holdout_eligible",
        )

    def test_blocking_reason_cannot_be_removed(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["eligibility"]["blocking_reasons"].pop()

        self._require_failure(
            self._run(mutate),
            "contract-list-length-mismatch:root.eligibility.blocking_reasons",
        )

    def test_candidate_output_requires_verified_source(self) -> None:
        self._require_failure(
            self._run(extra_args=["--candidate-output", "candidate.csv"]),
            "candidate-output-requires-source-verification",
        )

    def test_overwrite_requires_output_path(self) -> None:
        self._require_failure(
            self._run(extra_args=["--overwrite-candidate-output"]),
            "overwrite-requires-candidate-output",
        )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(KijSourceCandidateTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"RFAM_KIJ_SOURCE_CANDIDATE_TEST_PASS checks={KijSourceCandidateTests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
