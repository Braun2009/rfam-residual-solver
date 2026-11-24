#!/usr/bin/env python3
"""Focused tests for the blocked C3Mech R2 fit screen."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
TOOLS = ROOT / "tools"
CONTRACT = ROOT / "tests/reference/r2_c3mech_fit_screen_contract_v1.json"
TRAINING = WORKSPACE / "tmp/a408-r2-c3mech-44i0-training-v2.json"
HOLDOUT = WORKSPACE / "tmp/a408-r2-c3mech-44i0-holdout-v2.json"
SCREEN_TOOL = TOOLS / "run_r2_c3mech_fit_screen.py"

sys.path.insert(0, str(TOOLS))
import artifact_composite  # noqa: E402
import run_r2_c3mech_fit_screen as screen  # noqa: E402


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def reference_sample(
    *,
    qdot: float = 0.0,
    temperature: float = 1000.0,
    time_s: float = 0.0,
) -> dict:
    return {
        "density_kg_m3": 1.0,
        "heat_release_rate_w_m3": qdot,
        "major_species_mass_fractions": {
            "C3H8": 0.2,
            "O2": 0.8,
            "CO": 0.0,
            "CO2": 0.0,
            "H2O": 0.0,
            "N2": 0.0,
        },
        "pressure_pa": 1.0e6,
        "status": "eligible-within-mechanism-thermo-intersection",
        "temperature_k": temperature,
        "time_s": time_s,
    }


def prediction_sample(
    *,
    qdot: float = 0.0,
    temperature: float = 1000.0,
    time_s: float = 0.0,
) -> dict:
    return {
        "time_s": time_s,
        "temperature_k": temperature,
        "pressure_pa": 1.0e6,
        "density_kg_m3": 1.0,
        "sensible_enthalpy_j_kg": 0.0,
        "total_enthalpy_j_kg": 0.0,
        "qdot_w_m3": qdot,
        "mass_fractions": [0.2, 0.8, 0.0, 0.0, 0.0, 0.0],
    }


def prepared_case(
    case_id: str,
    reference_idt: float | None,
    *,
    qdot: float = 0.0,
) -> screen.PreparedCase:
    run_times = (0.0, 5.0e-4) if reference_idt is not None else (0.0,)
    return screen.PreparedCase(
        case_id=case_id,
        partition="protocol-training",
        temperature_k=1000.0,
        pressure_pa=1.0e6,
        initial_mass_fractions=(0.2, 0.8, 0.0, 0.0, 0.0, 0.0),
        eligible_samples=(reference_sample(qdot=qdot),),
        runner_times_s=run_times,
        reference_ignition_delay_s=reference_idt,
    )


class R2C3MechFitScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        screen.load_bound_helpers(
            WORKSPACE,
            ROOT / "tests/reference/r2_c3mech_fit_training_inputs_v1_composite.json",
        )

    def test_contract_schema_and_authorization_are_exact(self) -> None:
        document, raw = screen.load_json(CONTRACT, "contract")
        self.assertEqual(digest(raw), screen.CONTRACT_AUTHORITY["sha256"])
        self.assertEqual(len(raw), screen.CONTRACT_AUTHORITY["byte_count"])
        specs, config, starts = screen.validate_contract(document)
        self.assertEqual([spec.name for spec in specs], list(screen.PARAMETER_NAMES))
        self.assertEqual(config.maximum_evaluations, 800)
        self.assertEqual([kind for kind, _ in starts], [
            "midpoint", "halton-1", "halton-2", "halton-3", "halton-4",
        ])

        mutations = {
            "unknown-root-field": lambda value: value.__setitem__("unreviewed", True),
            "runner-phase": lambda value: value["runner"].__setitem__("phase", "constant-volume"),
            "history-reduction": lambda value: value["objective"].__setitem__(
                "history_component", "mean-signed-error"
            ),
            "holdout-timing": lambda value: value["input_composites"]["holdout"].__setitem__(
                "verified_only_after_training_selection_frozen", False
            ),
            "self-consistent-composite-rebind": lambda value: value["input_composites"]
            ["training_and_shared"].update({
                "manifest_sha256": "1" * 64,
                "external_composite_sha256": "2" * 64,
            }),
            "authorization": lambda value: value["approvals"].__setitem__(
                "production_ready", True
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                candidate = copy.deepcopy(document)
                mutate(candidate)
                with self.assertRaises(screen.ScreenError):
                    screen.validate_contract(candidate)

    def test_contract_redirection_is_rejected_before_any_composite_is_opened(self) -> None:
        document, _ = screen.load_json(CONTRACT, "contract")
        redirected = copy.deepcopy(document)
        redirected["input_composites"]["training_and_shared"]["manifest_path"] = (
            "tmp/unreviewed-composite.json"
        )
        with mock.patch.object(screen, "verify_runtime"), mock.patch.object(
            screen, "load_json", return_value=(redirected, b"redirected")
        ), mock.patch.object(screen, "verify_composite") as verifier:
            with self.assertRaisesRegex(screen.ScreenError, "contract-content-mismatch"):
                screen.execute(CONTRACT)
        verifier.assert_not_called()

    def test_duplicate_json_key_and_nonfinite_token_fail_closed(self) -> None:
        fixtures = {
            "duplicate": b'{"approvals":false,"approvals":true}\n',
            "nan": b'{"value":NaN}\n',
        }
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-json-") as directory:
            for label, raw in fixtures.items():
                with self.subTest(label=label):
                    path = Path(directory) / (label + ".json")
                    path.write_bytes(raw)
                    with self.assertRaises(screen.ScreenError):
                        screen.load_json(path, label)

    def test_split_references_produce_per_case_eligible_prefix_and_sentinel(self) -> None:
        training_document, _ = screen.load_json(TRAINING, "training")
        training, training_times = screen.validate_reference(
            training_document,
            "protocol-training",
            screen.EXPECTED_PARTITIONS["training_case_ids"],
        )
        self.assertEqual([len(case.eligible_samples) for case in training], [9, 8, 7])
        self.assertEqual(
            [len(case.runner_times_s) - len(case.eligible_samples) for case in training],
            [1, 1, 1],
        )
        self.assertEqual(len(training_times), 14)

        holdout_document, _ = screen.load_json(HOLDOUT, "holdout")
        holdout, _ = screen.validate_reference(
            holdout_document,
            "protocol-holdout",
            screen.EXPECTED_PARTITIONS["holdout_case_ids"],
        )
        self.assertEqual([len(case.eligible_samples) for case in holdout], [9, 7, 7])
        self.assertEqual(
            [len(case.runner_times_s) - len(case.eligible_samples) for case in holdout],
            [1, 1, 0],
        )
        self.assertIsNone(holdout[-1].reference_ignition_delay_s)

    def test_terminal_sentinel_prediction_is_never_scored(self) -> None:
        case = prepared_case("train_phi0p8_diluted", 2.0e-4, qdot=0.0)
        headers = {"ignition_detected": "1", "ignition_delay_s": "0.0002"}
        eligible_prediction = prediction_sample(qdot=1.0e8)
        normal_terminal = prediction_sample(
            qdot=0.0, temperature=1000.0, time_s=5.0e-4
        )
        extreme_terminal = prediction_sample(
            qdot=1.0e300, temperature=1.0e200, time_s=5.0e-4
        )
        first = screen.score_case(
            case,
            headers,
            [eligible_prediction, normal_terminal],
            screen.EXPECTED_OBJECTIVE["component_weights"],
        )
        second = screen.score_case(
            case,
            headers,
            [eligible_prediction, extreme_terminal],
            screen.EXPECTED_OBJECTIVE["component_weights"],
        )
        self.assertEqual(first, second)
        self.assertAlmostEqual(first[2]["qdot"], 1.0)

    def test_missing_idt_is_omitted_and_active_weights_are_renormalized(self) -> None:
        weights = screen.EXPECTED_OBJECTIVE["component_weights"]
        no_idt = prepared_case("holdout_phi1p4", None)
        no_idt_objective, _, no_idt_components = screen.score_case(
            no_idt,
            {"ignition_detected": "0", "ignition_delay_s": "-1"},
            [prediction_sample(qdot=1.0e8)],
            weights,
        )
        self.assertNotIn("ignition_delay", no_idt_components)
        self.assertAlmostEqual(no_idt_objective, 0.15 / 0.70, places=15)

        with_idt = prepared_case("holdout_phi1p1", 1.0e-4)
        with_idt_objective, _, with_idt_components = screen.score_case(
            with_idt,
            {"ignition_detected": "1", "ignition_delay_s": "0.0001"},
            [prediction_sample(qdot=1.0e8), prediction_sample(time_s=5.0e-4)],
            weights,
        )
        self.assertEqual(with_idt_components["ignition_delay"], 0.0)
        self.assertAlmostEqual(with_idt_objective, 0.15, places=15)

    def test_holdout_install_and_access_are_rejected_until_selection_freezes(self) -> None:
        document, _ = screen.load_json(CONTRACT, "contract")
        specs, _, _ = screen.validate_contract(document)
        training = [prepared_case(case_id, None) for case_id in document["partitions"]["training_case_ids"]]
        holdout = [prepared_case(case_id, None) for case_id in document["partitions"]["holdout_case_ids"]]
        objective = screen.R2C3MechObjective(Path("unused"), specs, training, document)
        with self.assertRaisesRegex(screen.ScreenError, "holdout-accessed-during-selection"):
            objective.evaluate([0.5] * len(specs), "holdout")
        with self.assertRaisesRegex(screen.ScreenError, "holdout-installed-during-selection"):
            objective.install_holdout(holdout)
        self.assertEqual(objective.holdout_access_attempts_during_selection, 1)

        objective.freeze_selection()
        objective.install_holdout(holdout)
        self.assertEqual(objective.case_ids["holdout"], tuple(document["partitions"]["holdout_case_ids"]))
        with self.assertRaisesRegex(screen.ScreenError, "holdout-installed-more-than-once"):
            objective.install_holdout(holdout)

    def test_invalid_search_point_never_replaces_valid_incumbent(self) -> None:
        spec = screen.ParameterSpec("x", "linear", 0.0, 1.0, 0.0, 1.0)
        config = screen.SearchConfig(0.25, 0.5, 0.2, 10, 1)

        def evaluate(point):
            if point[0] == 0.75:
                return screen.ObjectiveEvaluation(False, -100.0, (), {}, "synthetic-invalid")
            return screen.ObjectiveEvaluation(True, 1.0, (1.0,), {"case": 1.0})

        result = screen.run_pattern_search(
            0, "test", (0.5,), (spec,), config, evaluate
        )
        self.assertTrue(result.valid)
        self.assertEqual(result.selected_u, (0.5,))
        self.assertEqual(result.evaluation.objective, 1.0)

    def test_composite_drift_is_rejected_even_when_manifest_pin_is_unchanged(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-composite-") as directory:
            root = Path(directory)
            artifact = root / "inputs/reference.json"
            artifact.parent.mkdir(parents=True)
            artifact.write_bytes(b'{"frozen":true}\n')
            manifest = root / "manifest.json"
            composite = artifact_composite.generate_manifest(
                root, manifest, ["inputs/reference.json"]
            )
            authority = {
                "manifest_path": "manifest.json",
                "manifest_sha256": digest(manifest.read_bytes()),
                "external_composite_sha256": composite,
                "artifact_count": 1,
                "verified_before_selection": True,
            }
            expected_module = root / screen.BOUND_PYTHON_AUTHORITIES["artifact_composite"]["path"]
            expected_module.parent.mkdir(parents=True)
            expected_module.write_bytes(Path(artifact_composite.__file__).read_bytes())
            with mock.patch.object(artifact_composite, "__file__", str(expected_module)):
                self.assertEqual(screen.verify_composite(root, authority, "fixture"), manifest)
                artifact.write_bytes(b'{"frozen":false}\n')
                with self.assertRaisesRegex(screen.ScreenError, "composite-verification-failed"):
                    screen.verify_composite(root, authority, "fixture")

    def test_fit_header_requires_exact_parameter_identity(self) -> None:
        header = "fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n"
        header += "".join(
            "fit_" + name + "=" + format(value, ".17g") + "\n"
            for name, value in zip(screen.PARAMETER_NAMES, screen.DEFAULT_PARAMETERS)
        )
        trajectory = "schema=rfam.r2-0d-rfam-trajectory/v1\nR2_0D_RFAM_RUNNER_OK\n"
        self.assertEqual(
            screen.strip_fit_parameter_header(
                header + trajectory, screen.DEFAULT_PARAMETERS, "case"
            ),
            trajectory,
        )
        with self.assertRaisesRegex(screen.ScreenError, "runner-fit-header-mismatch"):
            screen.strip_fit_parameter_header(
                (header + trajectory).replace("fit_R1_A_SI=10000000000", "fit_R1_A_SI=10000000001"),
                screen.DEFAULT_PARAMETERS,
                "case",
            )

    def test_cli_failure_preserves_existing_output_without_success_or_traceback(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-cli-") as directory:
            temporary = Path(directory)
            contract = temporary / "contract.json"
            contract.write_bytes(CONTRACT.read_bytes())
            output = temporary / "result.json"
            stale = b'{"stale":true}\n'
            output.write_bytes(stale)
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCREEN_TOOL),
                    "--contract",
                    str(contract),
                    "--output",
                    str(output),
                ],
                check=False,
                text=True,
                capture_output=True,
                timeout=10,
            )
            combined = result.stdout + result.stderr
            self.assertEqual(result.returncode, 2, combined)
            self.assertIn(screen.FAIL_MARKER, result.stderr)
            self.assertNotIn(screen.BLOCKED_MARKER, combined)
            self.assertNotIn("Traceback", combined)
            self.assertEqual(output.read_bytes(), stale)

    def test_cli_preflight_failure_creates_no_fresh_output(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-cli-") as directory:
            temporary = Path(directory)
            output = temporary / "result.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCREEN_TOOL),
                    "--contract",
                    str(temporary / "not-canonical.json"),
                    "--output",
                    str(output),
                ],
                check=False,
                text=True,
                capture_output=True,
                timeout=10,
            )
            combined = result.stdout + result.stderr
            self.assertEqual(result.returncode, 2, combined)
            self.assertIn(screen.FAIL_MARKER, result.stderr)
            self.assertNotIn(screen.BLOCKED_MARKER, combined)
            self.assertNotIn("Traceback", combined)
            self.assertFalse(output.exists())

    def test_atomic_publication_never_replaces_existing_target(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-output-") as directory:
            output = Path(directory) / "result.json"
            original = b'{"owner":"existing"}\n'
            output.write_bytes(original)
            with self.assertRaisesRegex(screen.ScreenError, "output-path-already-exists"):
                screen.write_atomic(output, b'{"owner":"replacement"}\n')
            self.assertEqual(output.read_bytes(), original)
            self.assertEqual(list(output.parent.glob(output.name + ".tmp.*")), [])

    def test_atomic_publication_loses_race_without_replacing_winner(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-output-") as directory:
            output = Path(directory) / "result.json"
            winner = b'{"owner":"race-winner"}\n'
            real_link = os.link

            def collide(source, target):
                Path(target).write_bytes(winner)
                return real_link(source, target)

            with mock.patch.object(screen.os, "link", side_effect=collide):
                with self.assertRaises(FileExistsError):
                    screen.write_atomic(output, b'{"owner":"candidate"}\n')
            self.assertEqual(output.read_bytes(), winner)
            self.assertEqual(list(output.parent.glob(output.name + ".tmp.*")), [])

    def test_output_symlink_alias_is_rejected_without_target_mutation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-r2-c3-output-") as directory:
            temporary = Path(directory)
            target = temporary / "authority.json"
            target.write_bytes(b'{"frozen":true}\n')
            output = temporary / "result.json"
            try:
                output.symlink_to(target)
            except OSError as error:
                self.skipTest("symlink unavailable: " + str(error))
            with self.assertRaisesRegex(screen.ScreenError, "output-path-already-exists"):
                screen.admit_output_path(output)
            self.assertEqual(target.read_bytes(), b'{"frozen":true}\n')


if __name__ == "__main__":
    unittest.main()
