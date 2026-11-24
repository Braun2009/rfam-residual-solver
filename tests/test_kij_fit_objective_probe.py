#!/usr/bin/env python3
"""Focused fail-close tests for the synthetic RFAM kij objective probe."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable


REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
TOOL = REPO / "tools/validate_kij_fit_objective_probe.py"
CONTRACT = REPO / "tests/reference/kij_fit_objective_contract_v1.json"
REFERENCE = REPO / "tests/reference/kij_fit_objective_synthetic_v1.csv"
BASE_TEXT = CONTRACT.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)

sys.path.insert(0, str(REPO / "tools"))
import validate_kij_fit_objective_probe as validator  # noqa: E402


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class KijFitObjectiveProbeTests(unittest.TestCase):
    checks = 0

    def _fixture(self, root: Path) -> Path:
        paths = set(validator.SOURCE_PATHS)
        paths.update(item["path"] for item in BASE_DOCUMENT["upstream"])
        paths.add(validator.REFERENCE_PATH)
        paths.add(validator.HARNESS_MANIFEST_PATH)
        manifest = json.loads((WORKSPACE / validator.HARNESS_MANIFEST_PATH).read_text(encoding="ascii"))
        paths.add(manifest["rights"]["terms_artifact"]["path"])
        paths.update(item["path"] for item in manifest["artifacts"])
        for relative in paths:
            source = WORKSPACE / relative
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        contract = root / validator.CONTRACT_PATH
        contract.parent.mkdir(parents=True, exist_ok=True)
        return contract

    def _run(
        self,
        mutate: Callable[[dict[str, Any]], None] | None = None,
        raw_text: str | None = None,
        mutate_workspace: Callable[[Path, dict[str, Any]], None] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-objective-") as temporary:
            root = Path(temporary)
            contract = self._fixture(root)
            document = copy.deepcopy(BASE_DOCUMENT)
            if mutate is not None:
                mutate(document)
            if mutate_workspace is not None:
                mutate_workspace(root, document)
            with contract.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(
                    raw_text
                    if raw_text is not None
                    else json.dumps(document, indent=2) + "\n"
                )
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

    def _require_failure(self, result: subprocess.CompletedProcess[str], reason: str) -> None:
        self.__class__.checks += 1
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, output)
        self.assertIn(validator.FAIL_MARKER, result.stdout)
        self.assertIn(reason, result.stdout)
        self.assertNotIn(validator.REFERENCE_ONLY_MARKER, output)
        self.assertNotIn(validator.PASS_MARKER, output)
        self.assertNotIn("Traceback", output)

    def _identities(self) -> list[tuple[str, str]]:
        _, _, identities = validator._manifest_and_input(WORKSPACE)
        return identities

    def _mutated_reference(self, mutate: Callable[[list[str]], None]) -> bytes:
        lines = REFERENCE.read_text(encoding="ascii").splitlines()
        mutate(lines)
        return ("\n".join(lines) + "\n").encode("ascii")

    def test_superseded_reference_fails_closed_after_source_drift(self) -> None:
        self._require_failure(self._run(), "source-bundle-sha256-mismatch")

    def test_duplicate_approval_key_is_rejected(self) -> None:
        needle = '    "fit_execution_approved": false,'
        replacement = '    "fit_execution_approved": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:fit_execution_approved",
        )

    def test_unknown_root_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("unreviewed", True)),
            "unknown-keys:root:unreviewed",
        )

    def test_fit_approval_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(lambda document: document["approvals"].__setitem__("fit_execution_approved", True)),
            "value-mismatch:root.approvals.fit_execution_approved",
        )

    def test_parameter_publication_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(lambda document: document["objective"].__setitem__("parameter_values_publishable", True)),
            "value-mismatch:root.objective.parameter_values_publishable",
        )

    def test_cantera_cannot_be_inserted_as_objective_authority(self) -> None:
        self._require_failure(
            self._run(lambda document: document["objective"].__setitem__("model_backend", "Cantera")),
            "unknown-keys:root.objective:model_backend",
        )

    def test_test_kij_values_are_frozen_as_nonphysical_endpoints(self) -> None:
        self._require_failure(
            self._run(lambda document: document["objective"].__setitem__("trial_test_kij", 0.055)),
            "value-mismatch:root.objective.trial_test_kij",
        )

    def test_split_policy_cannot_be_weakened(self) -> None:
        self._require_failure(
            self._run(lambda document: document["objective"].__setitem__("split_reporting", "combined")),
            "value-mismatch:root.objective.split_reporting",
        )

    def test_source_path_and_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            files = document["source_bundle"]["files"]
            files[0], files[1] = files[1], files[0]

        self._require_failure(
            self._run(mutate),
            "value-mismatch:root.source_bundle.files[0]",
        )

    def test_self_consistent_source_tamper_is_rejected(self) -> None:
        def mutate_workspace(root: Path, document: dict[str, Any]) -> None:
            source = root / validator.SOURCE_PATHS[0]
            source.write_bytes(source.read_bytes() + b"\n")
            document["source_bundle"]["sha256"] = validator.source_bundle_sha(root)

        self._require_failure(
            self._run(mutate_workspace=mutate_workspace),
            "value-mismatch:root.source_bundle.sha256",
        )

    def test_self_consistent_manifest_rehash_is_rejected(self) -> None:
        replacement = b"unreviewed synthetic manifest replacement\n"

        def mutate_workspace(root: Path, document: dict[str, Any]) -> None:
            (root / validator.HARNESS_MANIFEST_PATH).write_bytes(replacement)
            document["upstream"][1]["sha256"] = sha256(replacement)

        self._require_failure(
            self._run(mutate_workspace=mutate_workspace),
            "value-mismatch:root.upstream[1].sha256",
        )

    def test_self_consistent_reference_rehash_is_rejected(self) -> None:
        replacement = b"unreviewed reference replacement\n"

        def mutate_workspace(root: Path, document: dict[str, Any]) -> None:
            (root / validator.REFERENCE_PATH).write_bytes(replacement)
            document["reference"]["sha256"] = sha256(replacement)
            document["reference"]["byte_count"] = len(replacement)

        self._require_failure(
            self._run(mutate_workspace=mutate_workspace),
            "value-mismatch:root.reference.sha256",
        )

    def test_training_holdout_identity_relabel_is_rejected(self) -> None:
        payload = self._mutated_reference(
            lambda lines: lines.__setitem__(3, lines[3].replace(",training,", ",holdout,", 1))
        )
        with self.assertRaisesRegex(validator.ObjectiveContractError, "value-mismatch:output.identity"):
            validator.validate_output(payload, self._identities())
        self.__class__.checks += 1

    def test_negative_objective_component_is_rejected(self) -> None:
        def mutate(lines: list[str]) -> None:
            fields = lines[3].split(",")
            fields[17] = "-1"
            lines[3] = ",".join(fields)

        with self.assertRaisesRegex(validator.ObjectiveContractError, "negative-objective-component"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_total_objective_must_equal_weighted_components(self) -> None:
        def mutate(lines: list[str]) -> None:
            fields = lines[3].split(",")
            fields[22] = "1"
            lines[3] = ",".join(fields)

        with self.assertRaisesRegex(validator.ObjectiveContractError, "objective-sum-mismatch"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_rho_component_must_be_exercised(self) -> None:
        def mutate(lines: list[str]) -> None:
            for row in range(3, 6):
                fields = lines[row].split(",")
                old = float(fields[20])
                fields[20] = "0"
                fields[22] = repr(float(fields[22]) - old)
                lines[row] = ",".join(fields)
            summary = lines[-1]
            training = sum(float(lines[row].split(",")[22]) for row in (3, 4))
            holdout = float(lines[5].split(",")[22])
            lines[-1] = re.sub(r"training_objective=[^ ]+", f"training_objective={training!r}", summary)
            lines[-1] = re.sub(r"holdout_objective=[^ ]+", f"holdout_objective={holdout!r}", lines[-1])

        with self.assertRaisesRegex(validator.ObjectiveContractError, "objective-component-not-exercised"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_energy_component_must_be_exercised(self) -> None:
        def mutate(lines: list[str]) -> None:
            for row in range(3, 6):
                fields = lines[row].split(",")
                old = float(fields[21])
                fields[21] = "0"
                fields[22] = repr(float(fields[22]) - old)
                lines[row] = ",".join(fields)
            training = sum(float(lines[row].split(",")[22]) for row in (3, 4))
            holdout = float(lines[5].split(",")[22])
            lines[-1] = re.sub(r"training_objective=[^ ]+", f"training_objective={training!r}", lines[-1])
            lines[-1] = re.sub(r"holdout_objective=[^ ]+", f"holdout_objective={holdout!r}", lines[-1])

        with self.assertRaisesRegex(validator.ObjectiveContractError, "objective-component-not-exercised"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_nonfinite_output_is_rejected(self) -> None:
        def mutate(lines: list[str]) -> None:
            fields = lines[3].split(",")
            fields[14] = "nan"
            lines[3] = ",".join(fields)

        with self.assertRaisesRegex(validator.ObjectiveContractError, "nonfinite-number"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_material_closure_is_enforced(self) -> None:
        def mutate(lines: list[str]) -> None:
            fields = lines[3].split(",")
            fields[23] = "0.1"
            lines[3] = ",".join(fields)

        with self.assertRaisesRegex(validator.ObjectiveContractError, "material-residual"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_energy_closure_is_enforced(self) -> None:
        def mutate(lines: list[str]) -> None:
            fields = lines[3].split(",")
            fields[25] = "0.1"
            lines[3] = ",".join(fields)

        with self.assertRaisesRegex(validator.ObjectiveContractError, "energy-residual"):
            validator.validate_output(self._mutated_reference(mutate), self._identities())
        self.__class__.checks += 1

    def test_success_marker_cannot_claim_fit_execution(self) -> None:
        payload = self._mutated_reference(
            lambda lines: lines.__setitem__(-1, lines[-1].replace("fit_executed=0", "fit_executed=1"))
        )
        with self.assertRaisesRegex(validator.ObjectiveContractError, "summary-marker-mismatch"):
            validator.validate_output(payload, self._identities())
        self.__class__.checks += 1

    def test_substitute_executable_is_rejected_before_execution(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-probe-substitute-") as temporary:
            substitute = Path(temporary) / "probe"
            marker = Path(temporary) / "executed"
            with substitute.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(f"#!/bin/sh\ntouch '{marker}'\n")
            if os.name != "nt":
                substitute.chmod(0o755)
            with self.assertRaisesRegex(validator.ObjectiveContractError, "probe-executable-sha256-mismatch"):
                validator.run_probe(substitute, b"")
            self.assertFalse(marker.exists())
        self.__class__.checks += 1

    def test_frozen_executable_when_available(self) -> None:
        executable = Path(os.environ.get("RFAM_KIJ_FIT_OBJECTIVE_PROBE", "/tmp/rfam_kij_fit_objective_probe"))
        if not executable.exists():
            self.skipTest("frozen WSL probe executable is not available")
        _, input_payload, identities = validator._manifest_and_input(WORKSPACE)
        actual = validator.run_probe(executable, input_payload)
        reference = REFERENCE.read_bytes()
        validator.validate_output(actual, identities)
        self.__class__.checks += 1
        self.assertEqual(actual, reference)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(KijFitObjectiveProbeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"RFAM_KIJ_FIT_OBJECTIVE_TEST_PASS checks={KijFitObjectiveProbeTests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
