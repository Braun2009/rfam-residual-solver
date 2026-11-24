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
TOOL = REPO / "tools/validate_o2_c3h8_topology_envelope.py"
CONTRACT = REPO / "tests/reference/o2_c3h8_topology_envelope_contract_v1.json"
RISK_REGISTRY = REPO / "tests/reference/o2_c3h8_topology_model_risk_v1.json"
BASE_TEXT = CONTRACT.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)
RISK_TEXT = RISK_REGISTRY.read_text(encoding="ascii")

sys.path.insert(0, str(REPO / "tools"))
import validate_o2_c3h8_topology_envelope as validator  # noqa: E402


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def zero_offsets() -> dict[str, validator.Interval]:
    zero = validator.Interval(0.0, 0.0)
    return {
        "pressure_pa": zero,
        "temperature_k": zero,
        "kij": zero,
        "overall_z_o2": zero,
    }


def supported_box() -> dict[str, validator.Interval]:
    return {
        "pressure_pa": validator.Interval(0.0, 10.0e6),
        "temperature_k": validator.Interval(0.0, 1000.0),
        "kij": validator.Interval(-1.0, 1.0),
        "overall_z_o2": validator.Interval(-1.0, 2.0),
    }


class O2C3H8TopologyEnvelopeTests(unittest.TestCase):
    checks = 0

    def _fixture(self, root: Path) -> Path:
        required = {
            validator.TRAJECTORY_PATH,
            validator.RISK_PATH,
            validator.TOPOLOGY_REGRESSION_PATH,
            *validator.TRAJECTORY_SOURCE_PATHS,
            *validator.TOPOLOGY_SOURCE_PATHS,
        }
        for relative in sorted(required):
            source = WORKSPACE / relative
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        return root / validator.CONTRACT_PATH

    def _run(
        self,
        mutate_document: Callable[[dict[str, Any]], None] | None = None,
        mutate_workspace: Callable[[Path, dict[str, Any]], None] | None = None,
        raw_text: str | None = None,
        trajectory_probe: Path | None = None,
        topology_probe: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-topology-envelope-") as temporary:
            root = Path(temporary)
            contract = self._fixture(root)
            document = copy.deepcopy(BASE_DOCUMENT)
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
                sys.executable,
                str(TOOL),
                "--root",
                str(root),
                "--contract",
                str(contract),
            ]
            if trajectory_probe is not None:
                command.extend(["--trajectory-probe", str(trajectory_probe)])
            if topology_probe is not None:
                command.extend(["--topology-probe", str(topology_probe)])
            return subprocess.run(command, capture_output=True, text=True, check=False)

    def _require_failure(
        self, result: subprocess.CompletedProcess[str], reason: str
    ) -> None:
        self.__class__.checks += 1
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, output)
        self.assertIn(validator.FAIL_MARKER, result.stdout)
        self.assertIn(reason, result.stdout)
        for status in (
            validator.STATUS_THREE_PHASE,
            validator.STATUS_OUTSIDE,
            validator.STATUS_CLEAR,
        ):
            self.assertNotIn(
                f"RFAM_O2_C3H8_TOPOLOGY_ENVELOPE_{status}", output
            )
        self.assertNotIn("Traceback", output)

    @staticmethod
    def _anchor(
        temperature: float,
        pressure: float = 5.0e6,
        kij: float = 0.055,
        z_o2: float = 0.5,
        committable: bool = True,
    ) -> validator.Anchor:
        return validator.Anchor(
            identifier="synthetic-anchor",
            group="synthetic-group",
            order=0.0,
            pressure=pressure,
            temperature=temperature,
            kij=kij,
            z_o2=z_o2,
            committable=committable,
        )

    @staticmethod
    def _risk(temperature: float = 150.0) -> validator.RiskRecord:
        return validator.RiskRecord(
            identifier="synthetic-risk",
            pressure=5.0e6,
            temperature=temperature,
            kij=0.055,
            z_interval=validator.Interval(0.0, 1.0),
        )

    def test_nominal_safe_synthetic_tube_is_only_screen_clear(self) -> None:
        evidence = validator.evaluate_tube(
            [self._anchor(140.0)],
            [self._risk()],
            zero_offsets(),
            zero_offsets(),
            zero_offsets(),
            supported_box(),
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_CLEAR)
        self.assertEqual(evidence["risk_box_hits"], 0)
        self.assertFalse(evidence["physical_vlle_validated"])
        self.assertFalse(evidence["three_phase_owner_implemented"])
        self.assertFalse(evidence["production_ready"])

    def test_closed_interval_boundary_touch_requires_three_phase_owner(self) -> None:
        uncertainty = zero_offsets()
        uncertainty["temperature_k"] = validator.Interval(0.0, 1.0)
        evidence = validator.evaluate_tube(
            [self._anchor(149.0)],
            [self._risk()],
            uncertainty,
            zero_offsets(),
            zero_offsets(),
            supported_box(),
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_THREE_PHASE)
        self.assertEqual(evidence["risk_box_hits"], 1)

    def test_out_of_envelope_is_distinct_from_model_risk(self) -> None:
        supported = supported_box()
        supported["pressure_pa"] = validator.Interval(4.0e6, 6.0e6)
        evidence = validator.evaluate_tube(
            [self._anchor(200.0, pressure=7.0e6)],
            [self._risk()],
            zero_offsets(),
            zero_offsets(),
            zero_offsets(),
            supported,
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_OUTSIDE)
        self.assertEqual(evidence["risk_box_hits"], 0)
        self.assertEqual(evidence["out_of_envelope_boxes"], 1)

    def test_three_phase_status_precedes_simultaneous_out_of_envelope(self) -> None:
        supported = supported_box()
        supported["temperature_k"] = validator.Interval(0.0, 149.0)
        evidence = validator.evaluate_tube(
            [self._anchor(150.0)],
            [self._risk(150.0)],
            zero_offsets(),
            zero_offsets(),
            zero_offsets(),
            supported,
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_THREE_PHASE)
        self.assertEqual(evidence["risk_box_hits"], 1)
        self.assertEqual(evidence["out_of_envelope_boxes"], 1)

    def test_uncertainty_widening_cannot_remove_a_model_risk_hit(self) -> None:
        narrow = zero_offsets()
        narrow["temperature_k"] = validator.Interval(-0.4, 0.4)
        wide = zero_offsets()
        wide["temperature_k"] = validator.Interval(-0.4, 1.0)
        narrow_evidence = validator.evaluate_tube(
            [self._anchor(149.0)],
            [self._risk()],
            narrow,
            zero_offsets(),
            zero_offsets(),
            supported_box(),
        )
        wide_evidence = validator.evaluate_tube(
            [self._anchor(149.0)],
            [self._risk()],
            wide,
            zero_offsets(),
            zero_offsets(),
            supported_box(),
        )
        self.__class__.checks += 1
        self.assertEqual(narrow_evidence["risk_box_hits"], 0)
        self.assertGreaterEqual(
            wide_evidence["risk_box_hits"], narrow_evidence["risk_box_hits"]
        )
        self.assertEqual(wide_evidence["status"], validator.STATUS_THREE_PHASE)

    def test_each_nondegenerate_box_samples_lower_midpoint_upper(self) -> None:
        uncertainty = {
            name: validator.Interval(-1.0, 1.0) for name in zero_offsets()
        }
        evidence = validator.evaluate_tube(
            [self._anchor(140.0)],
            [self._risk(900.0)],
            uncertainty,
            zero_offsets(),
            zero_offsets(),
            supported_box(),
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["sample_points"], 3**4)

    def test_frozen_trajectory_fails_closed_after_caloric_source_drift(self) -> None:
        self._require_failure(self._run(), "trajectory-source-bundle-mismatch")

    def test_duplicate_approval_key_fails_without_success_marker(self) -> None:
        needle = '    "production_ready": false\n'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:production_ready",
        )

    def test_duplicate_risk_provenance_key_is_rejected(self) -> None:
        needle = '    "probe_source_path": '
        replacement = '    "probe_source_path": "unreviewed.cpp",\n' + needle
        self.assertEqual(RISK_TEXT.count(needle), 1)
        with tempfile.TemporaryDirectory(prefix="rfam-topology-risk-duplicate-") as temporary:
            path = Path(temporary) / "risk.json"
            with path.open("w", encoding="ascii", newline="\n") as stream:
                stream.write(RISK_TEXT.replace(needle, replacement, 1))
            with self.assertRaisesRegex(
                validator.EnvelopeError, "duplicate-json-key:probe_source_path"
            ):
                validator.validate_risk_registry(WORKSPACE, path)
        self.__class__.checks += 1

    def test_authority_path_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            trajectory = document["authorities"]["trajectory_reference"]
            risk = document["authorities"]["model_risk_registry"]
            for key in ("path", "sha256", "byte_count"):
                trajectory[key], risk[key] = risk[key], trajectory[key]

        self._require_failure(
            self._run(mutate_document=mutate),
            "value-mismatch:trajectory-authority.path",
        )

    def test_synchronized_risk_payload_rehash_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            path = root / validator.RISK_PATH
            risk = json.loads(path.read_text(encoding="ascii"))
            risk["records"][0]["temperature_k"] += 0.1
            payload = (json.dumps(risk, indent=2) + "\n").encode("ascii")
            path.write_bytes(payload)
            authority = document["authorities"]["model_risk_registry"]
            authority["sha256"] = sha256(payload)
            authority["byte_count"] = len(payload)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "value-mismatch:risk-authority.sha256",
        )

    def test_topology_regression_append_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            regression = root / validator.TOPOLOGY_REGRESSION_PATH
            regression.write_bytes(regression.read_bytes() + b"\n")

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "topology-regression-byte-count-mismatch",
        )

    def test_missing_topology_regression_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            (root / validator.TOPOLOGY_REGRESSION_PATH).unlink()

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "topology-regression-read-failed",
        )

    def test_synchronized_topology_regression_rehash_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            regression = root / validator.TOPOLOGY_REGRESSION_PATH
            regression_payload = regression.read_bytes() + b"\n"
            regression.write_bytes(regression_payload)

            risk_path = root / validator.RISK_PATH
            risk = json.loads(risk_path.read_text(encoding="ascii"))
            risk["backend"]["regression_sha256"] = sha256(regression_payload)
            risk["backend"]["regression_byte_count"] = len(regression_payload)
            risk_payload = (json.dumps(risk, indent=2) + "\n").encode("ascii")
            risk_path.write_bytes(risk_payload)

            authority = document["authorities"]["model_risk_registry"]
            authority["sha256"] = sha256(risk_payload)
            authority["byte_count"] = len(risk_payload)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "value-mismatch:risk-authority.sha256",
        )

    def test_synchronized_topology_regression_path_rebind_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            regression = root / validator.TOPOLOGY_REGRESSION_PATH
            alternate_relative = (
                "realFluidAllMach/tests/unreviewed_topology_regression.py"
            )
            alternate = root / alternate_relative
            alternate.write_bytes(regression.read_bytes())

            risk_path = root / validator.RISK_PATH
            risk = json.loads(risk_path.read_text(encoding="ascii"))
            risk["backend"]["regression_path"] = alternate_relative
            risk_payload = (json.dumps(risk, indent=2) + "\n").encode("ascii")
            risk_path.write_bytes(risk_payload)

            authority = document["authorities"]["model_risk_registry"]
            authority["sha256"] = sha256(risk_payload)
            authority["byte_count"] = len(risk_payload)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "value-mismatch:risk-authority.sha256",
        )

    def test_source_list_and_matching_bundle_hash_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            build = document["authorities"]["trajectory_probe_build"]
            sources = build["source_files"]
            sources[0], sources[1] = sources[1], sources[0]
            build["source_bundle_sha256"] = validator.source_bundle_sha(
                WORKSPACE, tuple(sources)
            )

        self._require_failure(
            self._run(mutate_document=mutate),
            "value-mismatch:trajectory_probe_build.source_files",
        )

    def test_changed_source_and_matching_bundle_hash_are_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            source = root / validator.TRAJECTORY_SOURCE_PATHS[0]
            source.write_bytes(source.read_bytes() + b"\n")
            build = document["authorities"]["trajectory_probe_build"]
            build["source_bundle_sha256"] = validator.source_bundle_sha(
                root, validator.TRAJECTORY_SOURCE_PATHS
            )

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "value-mismatch:trajectory_probe_build.source_bundle_sha256",
        )

    def test_substitute_executable_and_matching_contract_hash_fail_pre_exec(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-topology-substitute-") as temporary:
            marker = Path(temporary) / "launched.txt"
            substitute = Path(temporary) / "trajectory-probe.cmd"
            with substitute.open("w", encoding="ascii", newline="") as stream:
                stream.write(f"@echo launched>{marker}\r\n")
            substitute_hash = sha256(substitute.read_bytes())

            def mutate(document: dict[str, Any]) -> None:
                document["authorities"]["trajectory_probe_build"][
                    "expected_executable_sha256"
                ] = substitute_hash

            self._require_failure(
                self._run(mutate_document=mutate, trajectory_probe=substitute),
                "value-mismatch:trajectory_probe_build.expected_executable_sha256",
            )
            self.assertFalse(marker.exists())


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        O2C3H8TopologyEnvelopeTests
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(
            "RFAM_O2_C3H8_TOPOLOGY_ENVELOPE_TEST_PASS "
            f"checks={O2C3H8TopologyEnvelopeTests.checks}"
        )
    raise SystemExit(0 if result.wasSuccessful() else 1)
