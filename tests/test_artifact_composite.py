#!/usr/bin/env python3
"""Focused fail-close tests for the artifact composite toolkit."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List


REPO = Path(__file__).resolve().parents[1]
TOOL = REPO / "tools/artifact_composite.py"
sys.path.insert(0, str(REPO / "tools"))
import artifact_composite as composite  # noqa: E402


def canonical_bytes(document: Dict[str, Any]) -> bytes:
    text = json.dumps(document, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    return text.encode("ascii")


def manual_composite(root: Path, paths: List[str]) -> str:
    digest = hashlib.sha256()
    for relative in sorted(paths):
        payload = root.joinpath(*relative.split("/")).read_bytes()
        digest.update(relative.encode("ascii"))
        digest.update(b"\0")
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


class ArtifactCompositeTests(unittest.TestCase):
    checks = 0

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rfam-composite-")
        self.root = Path(self.temporary.name)
        (self.root / "nested").mkdir()
        (self.root / "alpha.txt").write_bytes(b"alpha\n")
        (self.root / "nested/beta.bin").write_bytes(bytes(range(32)))
        self.paths = ["nested/beta.bin", "alpha.txt"]
        self.manifest = self.root / "freeze/manifest.json"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_tool(
        self, *arguments: str
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(TOOL), *arguments],
            capture_output=True,
            check=False,
            text=True,
        )

    def generate(self, print_composite: bool = False) -> subprocess.CompletedProcess[str]:
        arguments = [
            "generate",
            "--root",
            str(self.root),
            "--manifest",
            str(self.manifest),
        ]
        for path in self.paths:
            arguments.extend(["--artifact", path])
        if print_composite:
            arguments.append("--print-composite")
        return self.run_tool(*arguments)

    def verify(
        self, expected: str, print_composite: bool = False
    ) -> subprocess.CompletedProcess[str]:
        arguments = [
            "verify",
            "--root",
            str(self.root),
            "--manifest",
            str(self.manifest),
            "--expected-composite",
            expected,
        ]
        if print_composite:
            arguments.append("--print-composite")
        return self.run_tool(*arguments)

    def document(self) -> Dict[str, Any]:
        return json.loads(self.manifest.read_text(encoding="ascii"))

    def rewrite(self, document: Dict[str, Any]) -> None:
        self.manifest.write_bytes(canonical_bytes(document))

    def require_failure(self, result: subprocess.CompletedProcess[str], reason: str) -> None:
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertNotIn("Traceback", result.stderr)
        line = result.stderr.strip()
        payload = json.loads(line)
        self.assertEqual(payload["schema"], composite.RESULT_SCHEMA)
        self.assertEqual(payload["status"], "FAIL")
        self.assertIn(reason, payload["reason"])

    def test_generate_is_deterministic_canonical_and_silent(self) -> None:
        result = self.generate()
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr, "")
        document = self.document()
        expected = manual_composite(self.root, self.paths)
        self.assertEqual(document["composite_sha256"], expected)
        self.assertEqual(document["artifact_count"], 2)
        self.assertEqual(
            [item["path"] for item in document["artifacts"]], sorted(self.paths)
        )
        for item in document["artifacts"]:
            self.assertEqual(item["id"], "artifact:" + item["path"])
        self.assertEqual(self.manifest.read_bytes(), canonical_bytes(document))

        first = self.manifest.read_bytes()
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.manifest.read_bytes(), first)

    def test_print_composite_is_explicit_opt_in(self) -> None:
        expected = manual_composite(self.root, self.paths)
        result = self.generate(print_composite=True)
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, expected + "\n")
        self.assertEqual(result.stderr, "")
        result = self.verify(expected, print_composite=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, expected + "\n")

    def test_verify_success_is_silent(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        result = self.verify(expected)
        self.__class__.checks += 1
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr, "")

    def test_duplicate_json_key_is_rejected(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        text = self.manifest.read_text(encoding="ascii")
        needle = '  "schema": "rfam.artifact-composite/v1"\n'
        self.assertEqual(text.count(needle), 1)
        with self.manifest.open("w", encoding="ascii", newline="\n") as stream:
            stream.write(
                text.replace(needle, '  "schema": "unreviewed",\n' + needle, 1)
            )
        self.require_failure(
            self.verify(manual_composite(self.root, self.paths)),
            "duplicate-json-key:schema",
        )

    def test_schema_is_closed_and_contract_version_is_strict_integer(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        document = self.document()
        document["unreviewed"] = True
        self.rewrite(document)
        self.require_failure(self.verify(expected), "unknown-keys:root:unreviewed")

        self.assertEqual(self.generate().returncode, 0)
        document = self.document()
        document["contract_version"] = True
        self.rewrite(document)
        self.require_failure(self.verify(expected), "contract-version-mismatch")

    def test_artifact_schema_and_id_path_binding_are_closed(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        document = self.document()
        document["artifacts"][0]["unreviewed"] = True
        self.rewrite(document)
        self.require_failure(self.verify(expected), "unknown-keys:artifacts[0]")

        self.assertEqual(self.generate().returncode, 0)
        document = self.document()
        document["artifacts"][0]["id"] = document["artifacts"][1]["id"]
        self.rewrite(document)
        self.require_failure(
            self.verify(expected), "artifact-id-path-binding-mismatch"
        )

    def test_duplicate_case_alias_and_unsorted_paths_are_rejected(self) -> None:
        result = self.run_tool(
            "generate",
            "--root",
            str(self.root),
            "--manifest",
            str(self.manifest),
            "--artifact",
            "alpha.txt",
            "--artifact",
            "ALPHA.TXT",
        )
        self.require_failure(result, "duplicate-or-case-aliased-path")

        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        document = self.document()
        document["artifacts"].reverse()
        self.rewrite(document)
        self.require_failure(self.verify(expected), "artifacts-not-strictly-path-sorted")

    def test_path_spelling_matches_the_filesystem_on_case_insensitive_roots(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix=".rfam-composite-case-", dir=str(REPO)
        ) as directory:
            root = Path(directory)
            (root / "CaseDir").mkdir()
            (root / "CaseDir/Artifact.bin").write_bytes(b"identity\n")
            alias = root / "casedir/artifact.BIN"
            if not alias.exists():
                self.skipTest("repository filesystem is case-sensitive")
            manifest = root / "manifest.json"
            result = self.run_tool(
                "generate",
                "--root",
                str(root),
                "--manifest",
                str(manifest),
                "--artifact",
                "casedir/artifact.BIN",
            )
            self.require_failure(result, "path-case-mismatch")
            self.assertFalse(manifest.exists())

    def test_duplicate_physical_file_is_rejected_when_hardlinks_are_supported(self) -> None:
        alias = self.root / "alpha-hardlink.txt"
        try:
            os.link(str(self.root / "alpha.txt"), str(alias))
        except (OSError, NotImplementedError) as exc:
            self.skipTest("hardlink creation unavailable: %s" % exc)
        result = self.run_tool(
            "generate",
            "--root",
            str(self.root),
            "--manifest",
            str(self.manifest),
            "--artifact",
            "alpha.txt",
            "--artifact",
            "alpha-hardlink.txt",
        )
        self.require_failure(result, "duplicate-physical-file")

    def test_algorithm_golden_vector(self) -> None:
        golden_root = self.root / "golden"
        (golden_root / "dir").mkdir(parents=True)
        (golden_root / "a.txt").write_bytes(b"A")
        (golden_root / "dir/b.bin").write_bytes(b"\x00\xff")
        golden_manifest = self.root / "golden-manifest.json"
        actual = composite.generate_manifest(
            golden_root, golden_manifest, ["a.txt", "dir/b.bin"]
        )
        self.__class__.checks += 1
        self.assertEqual(
            actual,
            "b1e5f0766f68a8b8f5614d4a3cb302e417268c1fec4e2eb704db895b63ba59cd",
        )
        self.assertEqual(
            composite.verify_manifest(golden_root, golden_manifest, actual), actual
        )

    def test_root_escape_and_nonportable_paths_are_rejected(self) -> None:
        for path, reason in (
            ("../outside.txt", "non-canonical-path-segment"),
            ("nested\\beta.bin", "non-relative-posix-path"),
            ("nested//beta.bin", "empty-path-segment"),
            ("NUL.txt", "windows-reserved-path-segment"),
        ):
            with self.subTest(path=path):
                result = self.run_tool(
                    "generate",
                    "--root",
                    str(self.root),
                    "--manifest",
                    str(self.manifest),
                    "--artifact",
                    path,
                )
                self.require_failure(result, reason)

    def test_file_size_and_hash_drift_fail_closed(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        (self.root / "alpha.txt").write_bytes(b"longer payload\n")
        self.require_failure(self.verify(expected), "artifact-size-mismatch:alpha.txt")

        (self.root / "alpha.txt").write_bytes(b"alpha\n")
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        (self.root / "alpha.txt").write_bytes(b"ALPHA\n")
        self.require_failure(self.verify(expected), "artifact-sha256-mismatch:alpha.txt")

    def test_manifest_size_and_hash_tamper_fail_closed(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        document = self.document()
        document["artifacts"][0]["size"] += 1
        self.rewrite(document)
        self.require_failure(self.verify(expected), "artifact-size-mismatch")

        self.assertEqual(self.generate().returncode, 0)
        document = self.document()
        document["artifacts"][0]["sha256"] = "1" * 64
        self.rewrite(document)
        self.require_failure(self.verify(expected), "artifact-sha256-mismatch")

    def test_external_expected_composite_rejects_self_consistent_rehash(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        frozen = self.document()["composite_sha256"]
        (self.root / "alpha.txt").write_bytes(b"replacement\n")
        self.assertEqual(self.generate().returncode, 0)
        replacement = self.document()["composite_sha256"]
        self.assertNotEqual(replacement, frozen)
        self.require_failure(self.verify(frozen), "expected-composite-mismatch")

    def test_external_pin_rejects_self_consistent_path_id_and_hash_rebind(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        frozen = self.document()["composite_sha256"]

        # Keep the bytes identical while moving their semantic path. A freshly
        # generated manifest updates path, derived ID, per-file hash, and the
        # embedded composite together; the external pin must still reject it.
        (self.root / "replacement.txt").write_bytes(b"alpha\n")
        self.paths = ["nested/beta.bin", "replacement.txt"]
        self.assertEqual(self.generate().returncode, 0)
        rebound = self.document()
        self.assertEqual(
            [item["id"] for item in rebound["artifacts"]],
            ["artifact:" + item["path"] for item in rebound["artifacts"]],
        )
        self.assertNotEqual(rebound["composite_sha256"], frozen)
        self.require_failure(self.verify(frozen), "expected-composite-mismatch")

    def test_nested_unicode_escaped_duplicate_key_is_rejected(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        text = self.manifest.read_text(encoding="ascii")
        needle = '      "path": "alpha.txt",\n'
        self.assertEqual(text.count(needle), 1)
        replacement = (
            '      "pa\\u0074h": "replacement.txt",\n' + needle
        )
        with self.manifest.open("w", encoding="ascii", newline="\n") as stream:
            stream.write(text.replace(needle, replacement, 1))
        self.require_failure(self.verify(expected), "duplicate-json-key:path")

    def test_embedded_and_cli_composites_are_validated(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        self.require_failure(self.verify("x" * 64), "invalid-sha256:expected-composite")

        document = self.document()
        document["composite_sha256"] = "0" * 64
        self.rewrite(document)
        self.require_failure(self.verify(expected), "invalid-sha256:composite_sha256")

    def test_noncanonical_manifest_json_is_rejected(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        document = self.document()
        with self.manifest.open("w", encoding="ascii", newline="\n") as stream:
            stream.write(json.dumps(document, separators=(",", ":")) + "\n")
        self.require_failure(self.verify(expected), "manifest-not-canonical-json")

    def test_manifest_cannot_include_itself(self) -> None:
        self.manifest.parent.mkdir(parents=True)
        self.manifest.write_bytes(b"old manifest\n")
        relative = self.manifest.relative_to(self.root).as_posix()
        result = self.run_tool(
            "generate",
            "--root",
            str(self.root),
            "--manifest",
            str(self.manifest),
            "--artifact",
            relative,
        )
        self.require_failure(result, "manifest-cannot-be-an-artifact")

    def test_symlink_artifact_fails_closed_when_supported(self) -> None:
        target = self.root / "alpha.txt"
        link = self.root / "alpha-link.txt"
        try:
            os.symlink(str(target), str(link))
        except (OSError, NotImplementedError) as exc:
            self.skipTest("symlink creation unavailable: %s" % exc)
        result = self.run_tool(
            "generate",
            "--root",
            str(self.root),
            "--manifest",
            str(self.manifest),
            "--artifact",
            "alpha-link.txt",
        )
        self.require_failure(result, "symlink-or-reparse-point")

    def test_symlink_manifest_fails_closed_when_supported(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        expected = self.document()["composite_sha256"]
        link = self.root / "manifest-link.json"
        try:
            os.symlink(str(self.manifest), str(link))
        except (OSError, NotImplementedError) as exc:
            self.skipTest("symlink creation unavailable: %s" % exc)
        original = self.manifest
        self.manifest = link
        try:
            self.require_failure(self.verify(expected), "symlink-or-reparse-point:manifest")
        finally:
            self.manifest = original

    def test_argument_errors_are_structured(self) -> None:
        self.require_failure(self.run_tool("verify"), "argument-error")


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ArtifactCompositeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print("ARTIFACT_COMPOSITE_TEST_PASS checks=%d" % ArtifactCompositeTests.checks)
    raise SystemExit(0 if result.wasSuccessful() else 1)
