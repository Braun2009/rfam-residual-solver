#!/usr/bin/env python3
"""Focused fail-close tests for the current-source WP3 kij screen."""

from __future__ import annotations

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
TOOL = REPO / "tools/validate_kij_sensitivity_v2.py"
MANIFEST = REPO / "tests/reference/kij_sensitivity_contract_v2.json"
REFERENCE = REPO / "tests/reference/kij_sensitivity_o2_c3h8_v2.csv"
BASE_TEXT = MANIFEST.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)

sys.path.insert(0, str(REPO / "tools"))
import validate_kij_sensitivity_v2 as validator  # noqa: E402


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class KijSensitivityContractV2Tests(unittest.TestCase):
    checks = 0

    def _fixture(self, root: Path) -> Path:
        required = set(validator.SOURCE_PATHS)
        required.update(path for path, _sha in validator.EXPECTED_UPSTREAM.values())
        required.add(BASE_DOCUMENT["reference"]["path"])
        for relative in sorted(required):
            source = WORKSPACE / relative
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        return root / "realFluidAllMach/tests/reference/kij_sensitivity_contract_v2.json"

    def _run(
        self,
        mutate_document: Callable[[dict[str, Any]], None] | None = None,
        mutate_reference: Callable[[list[str], dict[str, int]], None] | None = None,
        mutate_workspace: Callable[[Path, dict[str, Any]], None] | None = None,
        raw_text: str | None = None,
        probe_path: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-contract-") as temporary:
            root = Path(temporary)
            contract = self._fixture(root)
            document = json.loads(BASE_TEXT)
            if mutate_reference is not None:
                reference_path = root / document["reference"]["path"]
                lines = reference_path.read_text(encoding="ascii").splitlines()
                fields = lines[1].split(",")
                mutate_reference(lines, {name: index for index, name in enumerate(fields)})
                payload = ("\n".join(lines) + "\n").encode("ascii")
                reference_path.write_bytes(payload)
                document["reference"]["sha256"] = sha256(payload)
            if mutate_document is not None:
                mutate_document(document)
            if mutate_workspace is not None:
                mutate_workspace(root, document)
            contract.parent.mkdir(parents=True, exist_ok=True)
            with contract.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(
                    raw_text
                    if raw_text is not None
                    else json.dumps(document, indent=2) + "\n"
                )
            command = [
                sys.executable, str(TOOL), "--root", str(root),
                "--contract", str(contract),
            ]
            if probe_path is not None:
                command.extend(["--probe-executable", str(probe_path)])
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
        self.assertNotIn(validator.SUCCESS_MARKER, output)
        self.assertNotIn("Traceback", output)

    def test_frozen_reference_fails_closed_after_source_drift(self) -> None:
        self._require_failure(self._run(), "source-bundle-sha256-mismatch")

    def test_duplicate_approval_key_fails_before_false_override(self) -> None:
        needle = '    "production_ready": false,'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:production_ready",
        )

    def test_duplicate_upstream_key_is_rejected(self) -> None:
        needle = '      "id": "tp_ph_repair_owner_report",'
        replacement = '      "id": "unreviewed",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:id",
        )

    def test_duplicate_core_provenance_key_is_rejected(self) -> None:
        needle = (
            '      "sha256": '
            '"e16a747f84c6d6891c3154bb385cad8ce39c6233952d8b8e1ce4fee0e07130dc",'
        )
        replacement = '      "sha256": "' + "0" * 64 + '",\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:sha256",
        )

    def test_contract_version_bool_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("contract_version", True)),
            "contract-mismatch:contract_version",
        )

    def test_unknown_root_key_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document.__setitem__("unreviewed", True)),
            "unknown-keys:root:unreviewed",
        )

    def test_approval_cannot_be_enabled(self) -> None:
        self._require_failure(
            self._run(lambda document: document["approvals"].__setitem__("production_ready", True)),
            "contract-mismatch:approvals",
        )

    def test_kij_grid_cannot_be_narrowed(self) -> None:
        self._require_failure(
            self._run(lambda document: document["model"]["sweep"].__setitem__(0, 0.0)),
            "kij-sweep-contract-mismatch:0",
        )

    def test_state_grid_cannot_be_reordered(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["state_grid"][0], document["state_grid"][1] = (
                document["state_grid"][1], document["state_grid"][0]
            )

        self._require_failure(self._run(mutate), "state_grid[0].case_id")

    def test_source_paths_cannot_be_swapped_with_matching_bundle_hash(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            files = document["source_bundle"]["files"]
            files[0], files[1] = files[1], files[0]
            document["source_bundle"]["sha256"] = validator.source_bundle_sha(
                WORKSPACE, tuple(files)
            )

        self._require_failure(self._run(mutate), "contract-mismatch:source_bundle.files")

    def test_source_bundle_hash_tamper_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document["source_bundle"].__setitem__("sha256", "1" * 64)),
            "source-bundle-sha256-mismatch",
        )

    def test_changed_source_and_matching_bundle_hash_are_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            source = root / validator.SOURCE_PATHS[0]
            source.write_bytes(source.read_bytes() + b"\n")
            document["source_bundle"]["sha256"] = validator.source_bundle_sha(root)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "contract-mismatch:source_bundle.sha256",
        )

    def test_each_provenance_report_and_matching_hash_is_rejected(self) -> None:
        for index in (3, 4, 5):
            with self.subTest(index=index):
                def mutate(
                    root: Path, document: dict[str, Any], report_index: int = index
                ) -> None:
                    item = document["upstream"][report_index]
                    report = root / item["path"]
                    payload = report.read_bytes() + b"\n"
                    report.write_bytes(payload)
                    item["sha256"] = sha256(payload)

                self._require_failure(
                    self._run(mutate_workspace=mutate),
                    f"contract-mismatch:upstream[{index}].sha256",
                )

    def test_provenance_report_paths_and_hashes_cannot_be_swapped(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            left = document["upstream"][3]
            right = document["upstream"][5]
            left["path"], right["path"] = right["path"], left["path"]
            left["sha256"], right["sha256"] = right["sha256"], left["sha256"]

        self._require_failure(
            self._run(mutate_document=mutate),
            "contract-mismatch:upstream[3].path",
        )

    def test_tp_ph_source_and_matching_provenance_are_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            item = document["core_provenance"]["tp_ph_flash_source"]
            source = root / item["path"]
            payload = source.read_bytes() + b"\n"
            source.write_bytes(payload)
            item["sha256"] = sha256(payload)
            document["source_bundle"]["sha256"] = validator.source_bundle_sha(root)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "contract-mismatch:core_provenance.tp_ph_flash_source.sha256",
        )

    def test_caloric_source_and_matching_provenance_are_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            item = document["core_provenance"]["ideal_caloric_80k_source"]
            source = root / item["path"]
            payload = source.read_bytes() + b"\n"
            source.write_bytes(payload)
            item["sha256"] = sha256(payload)
            document["source_bundle"]["sha256"] = validator.source_bundle_sha(root)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "contract-mismatch:core_provenance.ideal_caloric_80k_source.sha256",
        )

    def test_core_provenance_path_and_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            left = document["core_provenance"]["tp_ph_flash_source"]
            right = document["core_provenance"]["ideal_caloric_80k_source"]
            left["path"], right["path"] = right["path"], left["path"]
            left["sha256"], right["sha256"] = right["sha256"], left["sha256"]

        self._require_failure(
            self._run(mutate_document=mutate),
            "contract-mismatch:core_provenance.tp_ph_flash_source.path",
        )

    def test_core_provenance_role_and_admission_cannot_be_weakened(self) -> None:
        for provenance_id, field in (
            ("tp_ph_flash_source", "admission"),
            ("ideal_caloric_80k_source", "role"),
        ):
            with self.subTest(provenance_id=provenance_id, field=field):
                def mutate(
                    document: dict[str, Any],
                    item_id: str = provenance_id,
                    item_field: str = field,
                ) -> None:
                    document["core_provenance"][item_id][item_field] = "unreviewed"

                self._require_failure(
                    self._run(mutate_document=mutate),
                    f"contract-mismatch:core_provenance.{provenance_id}.{field}",
                )

    def test_reference_hash_tamper_is_rejected(self) -> None:
        self._require_failure(
            self._run(lambda document: document["reference"].__setitem__("sha256", "1" * 64)),
            "reference-sha256-mismatch",
        )

    def test_valid_changed_csv_and_matching_reference_hash_are_rejected(self) -> None:
        def mutate(lines: list[str], fields: dict[str, int]) -> None:
            row = lines[2].split(",")
            row[fields["e_j_kg"]] = str(float(row[fields["e_j_kg"]]) + 1.0)
            lines[2] = ",".join(row)

        self._require_failure(
            self._run(mutate_reference=mutate),
            "contract-mismatch:reference.sha256",
        )

    def test_duplicate_row_with_self_consistent_hash_is_rejected(self) -> None:
        def mutate(lines: list[str], _fields: dict[str, int]) -> None:
            lines[3] = lines[2]

        self._require_failure(self._run(mutate_reference=mutate), "kij-grid-mismatch:row[1]")

    def test_supercritical_quality_cannot_be_invented(self) -> None:
        def mutate(lines: list[str], fields: dict[str, int]) -> None:
            for index in range(2, len(lines) - 1):
                row = lines[index].split(",")
                if row[0] == "o2_supercritical_control_6mpa":
                    row[fields["quality_defined"]] = "1"
                    row[fields["vapor_mass_fraction"]] = "1"
                    lines[index] = ",".join(row)
                    return
            raise AssertionError("control row not found")

        self._require_failure(
            self._run(mutate_reference=mutate), "supercritical-quality-invented"
        )

    def test_two_phase_quality_out_of_range_is_rejected(self) -> None:
        def mutate(lines: list[str], fields: dict[str, int]) -> None:
            row = lines[2].split(",")
            row[fields["vapor_mass_fraction"]] = "1.1"
            lines[2] = ",".join(row)

        self._require_failure(
            self._run(mutate_reference=mutate), "two-phase-quality-out-of-range"
        )

    def test_material_closure_failure_is_rejected(self) -> None:
        def mutate(lines: list[str], fields: dict[str, int]) -> None:
            row = lines[2].split(",")
            row[fields["material_residual"]] = "1e-6"
            lines[2] = ",".join(row)

        self._require_failure(
            self._run(mutate_reference=mutate), "material-closure-failed"
        )

    def test_control_must_be_independent_of_kij(self) -> None:
        def mutate(lines: list[str], fields: dict[str, int]) -> None:
            matching = [
                index for index in range(2, len(lines) - 1)
                if lines[index].startswith("o2_supercritical_control_6mpa,")
            ]
            row = lines[matching[-1]].split(",")
            row[fields["rho_kg_m3"]] = str(float(row[fields["rho_kg_m3"]]) + 0.01)
            row[fields["single_rho_o2_kg_m3"]] = row[fields["rho_kg_m3"]]
            lines[matching[-1]] = ",".join(row)

        self._require_failure(
            self._run(mutate_reference=mutate), "control-kij-dependence"
        )

    def test_probe_identity_is_checked_before_execution(self) -> None:
        self._require_failure(
            self._run(probe_path=MANIFEST), "probe-executable-sha256-mismatch"
        )

    def test_substitute_probe_and_matching_expected_hash_are_rejected_pre_exec(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-kij-substitute-") as temporary:
            substitute = Path(temporary) / "probe.sh"
            marker = Path(temporary) / "executed"
            with substitute.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(
                    "#!/bin/sh\n"
                    "touch '" + marker.as_posix() + "'\n"
                    "cat '" + REFERENCE.as_posix() + "'\n"
                )
            substitute.chmod(0o755)
            substitute_sha = sha256(substitute.read_bytes())

            def mutate(document: dict[str, Any]) -> None:
                document["build"]["expected_executable_sha256"] = substitute_sha

            self._require_failure(
                self._run(mutate_document=mutate, probe_path=substitute),
                "contract-mismatch:build.expected_executable_sha256",
            )
            self.assertFalse(marker.exists())
            self.__class__.checks += 1


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(KijSensitivityContractV2Tests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"RFAM_KIJ_SENSITIVITY_V2_TEST_PASS checks={KijSensitivityContractV2Tests.checks}")
    raise SystemExit(0 if result.wasSuccessful() else 1)
