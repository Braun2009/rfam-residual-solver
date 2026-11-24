#!/usr/bin/env python3
"""Focused fail-close tests for recovered UCSD internal-screen provenance."""

from __future__ import annotations

import copy
import hashlib
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
TOOL = REPO / "tools" / "validate_r2_ucsd_internal_screen_provenance.py"
CONTRACT = REPO / "tests" / "reference" / "r2_ucsd_internal_screen_provenance_v1.json"
BASE_TEXT = CONTRACT.read_text(encoding="utf-8")
BASE = json.loads(BASE_TEXT)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class R2UcsdInternalScreenProvenanceTests(unittest.TestCase):
    def run_contract(self, contract: Path, root: Path = WORKSPACE, *extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(TOOL), "--root", str(root), "--contract", str(contract), *extra],
            capture_output=True,
            text=True,
            check=False,
        )

    def mutated(self, mutation: Callable[[dict[str, Any]], None]) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(BASE)
        mutation(document)
        with tempfile.TemporaryDirectory(prefix="rfam-ucsd-contract-") as directory:
            contract = Path(directory) / "contract.json"
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            return self.run_contract(contract)

    def raw(self, text: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-ucsd-contract-") as directory:
            contract = Path(directory) / "contract.json"
            contract.write_text(text, encoding="utf-8")
            return self.run_contract(contract)

    def fixture(self, directory: Path) -> Path:
        for artifact in BASE["artifacts"]:
            source = WORKSPACE / artifact["path"]
            target = directory / artifact["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        contract = directory / "realFluidAllMach/tests/reference/r2_ucsd_internal_screen_provenance_v1.json"
        contract.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(CONTRACT, contract)
        return contract

    def require_failure(self, result: subprocess.CompletedProcess[str], needle: str) -> None:
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("R2_UCSD_INTERNAL_SCREEN_PROVENANCE_FAIL", result.stderr)
        self.assertIn(needle, result.stderr)
        self.assertNotIn("R2_UCSD_INTERNAL_SCREEN_PROVENANCE_BLOCKED", result.stdout + result.stderr)
        self.assertNotIn("Traceback", result.stdout + result.stderr)

    def test_canonical_contract_is_blocked(self) -> None:
        result = self.run_contract(CONTRACT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_UCSD_INTERNAL_SCREEN_PROVENANCE_BLOCKED checks=15", result.stdout)

    def test_every_approval_remains_false(self) -> None:
        for key in BASE["approvals"]:
            with self.subTest(key=key):
                self.require_failure(
                    self.mutated(lambda document, key=key: document["approvals"].__setitem__(key, True)),
                    f"contract.approvals.{key}",
                )

    def test_contract_version_bool_is_not_integer(self) -> None:
        self.require_failure(
            self.mutated(lambda document: document.__setitem__("contract_version", True)),
            "type-mismatch:contract.contract_version",
        )

    def test_duplicate_root_scope_is_rejected(self) -> None:
        needle = '  "scope": "provenance-recovery-only-no-fit-no-physical-validation",'
        replacement = '  "scope": "production-approved",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self.require_failure(self.raw(BASE_TEXT.replace(needle, replacement, 1)), "duplicate-json-key:scope")

    def test_duplicate_nested_approval_is_rejected(self) -> None:
        needle = '    "production_ready": false,'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self.require_failure(self.raw(BASE_TEXT.replace(needle, replacement, 1)), "duplicate-json-key:production_ready")

    def test_unknown_authorization_is_rejected(self) -> None:
        self.require_failure(
            self.mutated(lambda document: document.__setitem__("unreviewed_authorization", True)),
            "key-mismatch:contract",
        )

    def test_artifact_path_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            first = document["artifacts"][3]
            second = document["artifacts"][5]
            for key in ("path", "sha256", "size", "kind"):
                first[key], second[key] = second[key], first[key]

        self.require_failure(self.mutated(mutate), "contract.artifacts[3].path")

    def test_self_consistent_yaml_rehash_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-ucsd-fixture-") as directory:
            root = Path(directory)
            contract = self.fixture(root)
            yaml_path = root / BASE["artifacts"][3]["path"]
            yaml_path.write_bytes(yaml_path.read_bytes() + b"\n")
            document = json.loads(contract.read_text(encoding="utf-8"))
            document["artifacts"][3]["sha256"] = sha256(yaml_path)
            document["artifacts"][3]["size"] = yaml_path.stat().st_size
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            self.require_failure(self.run_contract(contract, root), "contract.artifacts[3].sha256")

    def test_self_consistent_generator_substitution_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-ucsd-fixture-") as directory:
            root = Path(directory)
            contract = self.fixture(root)
            generator = root / BASE["artifacts"][5]["path"]
            generator.write_text("print('substitute')\n", encoding="utf-8")
            document = json.loads(contract.read_text(encoding="utf-8"))
            document["artifacts"][5]["sha256"] = sha256(generator)
            document["artifacts"][5]["size"] = generator.stat().st_size
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            self.require_failure(self.run_contract(contract, root), "contract.artifacts[5].sha256")

    def test_self_consistent_protocol_module_substitution_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-ucsd-fixture-") as directory:
            root = Path(directory)
            contract = self.fixture(root)
            module = root / BASE["artifacts"][8]["path"]
            module.write_text("raise SystemExit('substitute')\n", encoding="utf-8")
            document = json.loads(contract.read_text(encoding="utf-8"))
            document["artifacts"][8]["sha256"] = sha256(module)
            document["artifacts"][8]["size"] = module.stat().st_size
            contract.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            self.require_failure(self.run_contract(contract, root), "contract.artifacts[8].sha256")

    def test_missing_reference_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-ucsd-fixture-") as directory:
            root = Path(directory)
            contract = self.fixture(root)
            (root / BASE["artifacts"][6]["path"]).unlink()
            self.require_failure(self.run_contract(contract, root), "artifact-missing:detailed_reference")

    def test_runtime_identity_rebinding_is_rejected(self) -> None:
        self.require_failure(
            self.mutated(lambda document: document["runtime_identity"].__setitem__("cantera_version", "3.2.0")),
            "contract.runtime_identity.cantera_version",
        )

    def test_rights_or_fit_status_cannot_advance(self) -> None:
        for key in ("fit_use_allowed", "redistribution_rights_clear"):
            with self.subTest(key=key):
                self.require_failure(
                    self.mutated(lambda document, key=key: document["conversion_rights"].__setitem__(key, True)),
                    f"contract.conversion_rights.{key}",
                )


if __name__ == "__main__":
    unittest.main()
