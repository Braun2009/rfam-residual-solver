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
TOOL = REPO / "tools/validate_o2_c3h8_topology_uncertainty_tube.py"
CONTRACT = (
    REPO
    / "tests/reference/o2_c3h8_topology_uncertainty_tube_contract_v1.json"
)
BASE_TEXT = CONTRACT.read_text(encoding="ascii")
BASE_DOCUMENT = json.loads(BASE_TEXT)

sys.path.insert(0, str(REPO / "tools"))
import validate_o2_c3h8_topology_uncertainty_tube as validator  # noqa: E402


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def wide_supported() -> validator.Box:
    return {
        "pressure_pa": validator.Interval(0.0, 10.0e6),
        "temperature_k": validator.Interval(0.0, 1000.0),
        "kij": validator.Interval(-1.0, 1.0),
        "overall_z_o2": validator.Interval(0.0, 1.0),
    }


class O2C3H8TopologyUncertaintyTubeTests(unittest.TestCase):
    checks = 0

    @classmethod
    def setUpClass(cls) -> None:
        cls.base = validator.load_base_validator(WORKSPACE)

    def _fixture(self, root: Path) -> Path:
        required = {
            validator.BASE_VALIDATOR_PATH,
            validator.BASE_CONTRACT_PATH,
            validator.TRAJECTORY_PATH,
            validator.RISK_PATH,
            self.base.TOPOLOGY_REGRESSION_PATH,
            *self.base.TRAJECTORY_SOURCE_PATHS,
            *self.base.TOPOLOGY_SOURCE_PATHS,
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
        alternate_contract_path: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(
            prefix="rfam-topology-uncertainty-tube-"
        ) as temporary:
            root = Path(temporary)
            contract = self._fixture(root)
            document = copy.deepcopy(BASE_DOCUMENT)
            if mutate_document is not None:
                mutate_document(document)
            if mutate_workspace is not None:
                mutate_workspace(root, document)
            if alternate_contract_path:
                contract = root / "unbound-contract.json"
            contract.parent.mkdir(parents=True, exist_ok=True)
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
                f"RFAM_O2_C3H8_TOPOLOGY_UNCERTAINTY_TUBE_{status}",
                output,
            )
        self.assertNotIn("Traceback", output)

    @staticmethod
    def _anchor(
        *,
        pressure: float = 5.0e6,
        temperature: float = 140.0,
        kij: float = 0.055,
        z_o2: float = 0.5,
        committable: bool = True,
        identifier: str = "synthetic-anchor",
        order: float = 0.0,
    ) -> validator.Anchor:
        return validator.Anchor(
            identifier=identifier,
            group="synthetic-group",
            order=order,
            pressure=pressure,
            temperature=temperature,
            kij=kij,
            z_o2=z_o2,
            committable=committable,
        )

    @staticmethod
    def _risk(
        *,
        pressure: float = 5.0e6,
        temperature: float = 150.0,
        kij: float = 0.055,
        z_low: float = 0.0,
        z_high: float = 1.0,
        identifier: str = "synthetic-risk",
    ) -> validator.RiskRecord:
        return validator.RiskRecord(
            identifier=identifier,
            pressure=pressure,
            temperature=temperature,
            kij=kij,
            z_interval=validator.Interval(z_low, z_high),
        )

    def _evaluate(
        self,
        anchor: validator.Anchor,
        risk: validator.RiskRecord,
        *,
        trial_offset: validator.Interval | None = None,
        inlet: validator.Box | None = None,
        numerical: validator.Box | None = None,
        detector: validator.Box | None = None,
        supported: validator.Box | None = None,
    ) -> dict[str, Any]:
        return validator.evaluate_tube(
            [anchor],
            [risk],
            validator.Interval(-1.0, 1.0),
            trial_offset or validator.Interval(0.0, 0.0),
            inlet or validator.zero_offsets(),
            numerical or validator.zero_offsets(),
            detector or validator.zero_offsets(),
            supported or wide_supported(),
        )

    def test_frozen_tube_fails_closed_with_stale_base_trajectory(self) -> None:
        self._require_failure(
            self._run(),
            "base-envelope-validation-failed:trajectory-source-bundle-mismatch",
        )

    def test_safe_synthetic_tube_is_screen_only_clear(self) -> None:
        evidence = self._evaluate(self._anchor(), self._risk(temperature=200.0))
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_CLEAR)
        self.assertEqual(evidence["risk_box_hits"], 0)
        self.assertFalse(evidence["sampling_is_admission_authority"])

    def test_interval_contact_blocks_even_when_three_point_sampling_misses(self) -> None:
        inlet = validator.zero_offsets()
        inlet["temperature_k"] = validator.Interval(0.0, 2.0)
        evidence = self._evaluate(
            self._anchor(temperature=149.0),
            self._risk(temperature=149.25),
            inlet=inlet,
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_THREE_PHASE)
        self.assertEqual(evidence["risk_box_hits"], 1)
        self.assertEqual(evidence["sampled_risk_hits"], 0)

    def test_trial_kij_closed_boundary_contact_blocks(self) -> None:
        detector = validator.zero_offsets()
        detector["kij"] = validator.Interval(-0.001, 0.0)
        evidence = self._evaluate(
            self._anchor(kij=0.05),
            self._risk(temperature=140.0, kij=0.055),
            trial_offset=validator.Interval(0.0, 0.004),
            detector=detector,
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_THREE_PHASE)

    def test_inlet_pressure_temperature_and_composition_axes_all_block(self) -> None:
        cases = (
            (
                "pressure",
                self._anchor(pressure=4.95e6, temperature=140.0),
                self._risk(pressure=5.0e6, temperature=140.0),
                "pressure_pa",
                validator.Interval(0.0, 50000.0),
            ),
            (
                "temperature",
                self._anchor(temperature=149.0),
                self._risk(temperature=150.0),
                "temperature_k",
                validator.Interval(0.0, 1.0),
            ),
            (
                "composition",
                self._anchor(temperature=140.0, z_o2=0.49),
                self._risk(
                    temperature=140.0, z_low=0.5, z_high=0.5
                ),
                "overall_z_o2",
                validator.Interval(0.0, 0.01),
            ),
        )
        for name, anchor, risk, axis, offset in cases:
            with self.subTest(axis=name):
                inlet = validator.zero_offsets()
                inlet[axis] = offset
                evidence = self._evaluate(anchor, risk, inlet=inlet)
                self.assertEqual(
                    evidence["status"], validator.STATUS_THREE_PHASE
                )
                self.assertEqual(evidence["risk_box_hits"], 1)
                self.__class__.checks += 1

    def test_each_expansion_component_is_monotone_for_contact(self) -> None:
        scenarios: list[tuple[str, dict[str, Any], dict[str, Any]]] = []

        narrow_inlet = validator.zero_offsets()
        narrow_inlet["temperature_k"] = validator.Interval(0.0, 0.5)
        wide_inlet = validator.zero_offsets()
        wide_inlet["temperature_k"] = validator.Interval(0.0, 1.0)
        scenarios.append(("inlet", {"inlet": narrow_inlet}, {"inlet": wide_inlet}))

        narrow_numerical = validator.zero_offsets()
        narrow_numerical["temperature_k"] = validator.Interval(0.0, 0.5)
        wide_numerical = validator.zero_offsets()
        wide_numerical["temperature_k"] = validator.Interval(0.0, 1.0)
        scenarios.append(
            (
                "numerical",
                {"numerical": narrow_numerical},
                {"numerical": wide_numerical},
            )
        )

        narrow_detector = validator.zero_offsets()
        narrow_detector["temperature_k"] = validator.Interval(-0.5, 0.0)
        wide_detector = validator.zero_offsets()
        wide_detector["temperature_k"] = validator.Interval(-1.0, 0.0)
        scenarios.append(
            (
                "detector",
                {"detector": narrow_detector},
                {"detector": wide_detector},
            )
        )

        scenarios.append(
            (
                "trial-kij",
                {"trial_offset": validator.Interval(0.0, 0.003)},
                {"trial_offset": validator.Interval(0.0, 0.005)},
            )
        )

        for name, narrow_kwargs, wide_kwargs in scenarios:
            with self.subTest(component=name):
                if name == "trial-kij":
                    anchor = self._anchor(temperature=140.0, kij=0.05)
                    risk = self._risk(temperature=140.0, kij=0.055)
                else:
                    anchor = self._anchor(temperature=149.0)
                    risk = self._risk(temperature=150.0)
                narrow = self._evaluate(anchor, risk, **narrow_kwargs)
                wide = self._evaluate(anchor, risk, **wide_kwargs)
                self.assertEqual(narrow["risk_box_hits"], 0)
                self.assertEqual(wide["risk_box_hits"], 1)
                self.assertEqual(wide["status"], validator.STATUS_THREE_PHASE)
                self.__class__.checks += 1

    def test_three_phase_precedes_simultaneous_out_of_envelope(self) -> None:
        supported = wide_supported()
        supported["temperature_k"] = validator.Interval(0.0, 149.0)
        evidence = self._evaluate(
            self._anchor(temperature=150.0),
            self._risk(temperature=150.0),
            supported=supported,
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_THREE_PHASE)
        self.assertEqual(evidence["risk_box_hits"], 1)
        self.assertEqual(evidence["out_of_envelope_boxes"], 1)

    def test_uncommittable_safe_anchor_is_out_of_envelope(self) -> None:
        evidence = self._evaluate(
            self._anchor(committable=False), self._risk(temperature=200.0)
        )
        self.__class__.checks += 1
        self.assertEqual(evidence["status"], validator.STATUS_OUTSIDE)
        self.assertEqual(evidence["uncommittable_anchors"], 1)

    def test_empty_risk_registry_fails_closed(self) -> None:
        with self.assertRaisesRegex(validator.TubeError, "empty-risk-registry"):
            validator.evaluate_tube(
                [self._anchor()],
                [],
                validator.Interval(-1.0, 1.0),
                validator.Interval(0.0, 0.0),
                validator.zero_offsets(),
                validator.zero_offsets(),
                validator.zero_offsets(),
                wide_supported(),
            )
        self.__class__.checks += 1

    def test_anchor_outside_declared_trial_range_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            validator.TubeError, "anchor-kij-outside-trial-range"
        ):
            validator.evaluate_tube(
                [self._anchor(kij=0.2)],
                [self._risk()],
                validator.Interval(-0.05, 0.15),
                validator.Interval(0.0, 0.0),
                validator.zero_offsets(),
                validator.zero_offsets(),
                validator.zero_offsets(),
                wide_supported(),
            )
        self.__class__.checks += 1

    def test_duplicate_approval_key_fails_without_status_publication(self) -> None:
        needle = '    "production_ready": false\n'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace(needle, replacement, 1)),
            "duplicate-json-key:production_ready",
        )

    def test_missing_component_is_rejected(self) -> None:
        self._require_failure(
            self._run(
                mutate_document=lambda document: document["tube_components"].pop(
                    "trial_kij_interval"
                )
            ),
            "value-mismatch:contract.tube-components",
        )

    def test_category_swap_is_rejected_even_when_total_width_is_unchanged(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            inlet = document["tube_components"]["inlet_condition_error"]
            detector = document["tube_components"]["topology_detector_margin"]
            inlet["temperature_k_offset"], detector["temperature_k_offset"] = (
                detector["temperature_k_offset"],
                inlet["temperature_k_offset"],
            )

        self._require_failure(
            self._run(mutate_document=mutate),
            "value-mismatch:contract.tube-components",
        )

    def test_trial_kij_cannot_be_relabelled_as_fitted(self) -> None:
        self._require_failure(
            self._run(
                mutate_document=lambda document: document["tube_components"][
                    "trial_kij_interval"
                ].__setitem__("kij_fitted", True)
            ),
            "value-mismatch:contract.tube-components",
        )

    def test_physical_validation_and_production_cannot_be_enabled(self) -> None:
        for field in ("physical_validation", "production_ready"):
            with self.subTest(field=field):
                self._require_failure(
                    self._run(
                        mutate_document=lambda document, field=field: document[
                            "approvals"
                        ].__setitem__(field, True)
                    ),
                    "value-mismatch:contract.approvals",
                )

    def test_base_validator_drift_and_synchronized_rehash_are_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            path = root / validator.BASE_VALIDATOR_PATH
            payload = path.read_bytes() + b"\n"
            path.write_bytes(payload)
            authority = document["authorities"]["base_topology_envelope_validator"]
            authority["sha256"] = sha256(payload)
            authority["byte_count"] = len(payload)

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "value-mismatch:contract.authorities",
        )

    def test_substitute_base_validator_is_rejected_before_import(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="rfam-topology-validator-marker-"
        ) as temporary:
            marker = Path(temporary) / "executed.txt"

            def mutate(root: Path, document: dict[str, Any]) -> None:
                path = root / validator.BASE_VALIDATOR_PATH
                prefix = (
                    "from pathlib import Path\n"
                    f"Path({str(marker)!r}).write_text('executed')\n"
                ).encode("ascii")
                expected_size = validator.EXPECTED_AUTHORITIES[
                    "base_topology_envelope_validator"
                ]["byte_count"]
                self.assertLess(len(prefix), expected_size)
                path.write_bytes(prefix + b"#" * (expected_size - len(prefix)))

            self._require_failure(
                self._run(mutate_workspace=mutate),
                "base_topology_envelope_validator-sha256-mismatch",
            )
            self.assertFalse(marker.exists())

    def test_base_contract_drift_without_rebind_is_rejected(self) -> None:
        def mutate(root: Path, document: dict[str, Any]) -> None:
            path = root / validator.BASE_CONTRACT_PATH
            path.write_bytes(path.read_bytes() + b"\n")

        self._require_failure(
            self._run(mutate_workspace=mutate),
            "base_topology_envelope_contract-byte-count-mismatch",
        )

    def test_trajectory_and_risk_authority_drift_are_rejected(self) -> None:
        for relative, reason in (
            (
                validator.TRAJECTORY_PATH,
                "trajectory_reference-byte-count-mismatch",
            ),
            (validator.RISK_PATH, "model_risk_registry-byte-count-mismatch"),
        ):
            with self.subTest(path=relative):
                def mutate(
                    root: Path,
                    document: dict[str, Any],
                    relative: str = relative,
                ) -> None:
                    path = root / relative
                    path.write_bytes(path.read_bytes() + b"\n")

                self._require_failure(
                    self._run(mutate_workspace=mutate), reason
                )

    def test_authority_path_swap_is_rejected(self) -> None:
        def mutate(document: dict[str, Any]) -> None:
            left = document["authorities"]["trajectory_reference"]
            right = document["authorities"]["model_risk_registry"]
            for field in ("path", "sha256", "byte_count"):
                left[field], right[field] = right[field], left[field]

        self._require_failure(
            self._run(mutate_document=mutate),
            "value-mismatch:contract.authorities",
        )

    def test_noncanonical_crlf_contract_is_rejected(self) -> None:
        self._require_failure(
            self._run(raw_text=BASE_TEXT.replace("\n", "\r\n")),
            "contract-not-canonical-lf",
        )

    def test_alternate_contract_path_is_rejected(self) -> None:
        self._require_failure(
            self._run(alternate_contract_path=True), "contract-path-mismatch"
        )


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        O2C3H8TopologyUncertaintyTubeTests
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(
            "RFAM_O2_C3H8_TOPOLOGY_UNCERTAINTY_TUBE_TEST_PASS "
            f"checks={O2C3H8TopologyUncertaintyTubeTests.checks}"
        )
    raise SystemExit(0 if result.wasSuccessful() else 1)
