#!/usr/bin/env python3
"""Focused contract tests for the WSL C3Mech/RFAM route screen."""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPOSITORY_ROOT.parent
TOOL = REPOSITORY_ROOT / "tools/run_r2_c3mech_wsl_screen.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("rfam_r2_c3mech_wsl_screen", TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError("tool-import-spec-failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class C3MechWslScreenTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tool = load_tool()

    def synthetic_reference(self):
        composition = {
            "C3H8": 0.1,
            "O2": 0.9,
            "CO": 0.0,
            "CO2": 0.0,
            "H2O": 0.0,
            "N2": 0.0,
        }
        cases = []
        ids = (
            self.tool.CASE_IDS["protocol-training"]
            + self.tool.CASE_IDS["protocol-holdout"]
        )
        for case_index, case_id in enumerate(ids):
            samples = []
            for index, time_s in enumerate(self.tool.SAMPLE_TIMES):
                if index < 2:
                    samples.append(
                        {
                            "time_s": time_s,
                            "status": "eligible-within-mechanism-thermo-intersection",
                            "temperature_k": 1200.0 + case_index + index,
                            "pressure_pa": 2.0e6,
                            "density_kg_m3": 10.0,
                            "heat_release_rate_w_m3": -1.0e8,
                            "mass_fractions": dict(composition),
                        }
                    )
                elif index == 2:
                    samples.append(
                        {
                            "time_s": time_s,
                            "status": "outside-thermo-envelope-unscored",
                        }
                    )
                else:
                    samples.append(
                        {
                            "time_s": time_s,
                            "status": "not-evaluated-after-thermo-envelope-exit",
                        }
                    )
            cases.append(
                {
                    "id": case_id,
                    "partition": (
                        "protocol-training" if case_index < 3 else "protocol-holdout"
                    ),
                    "temperature_k": 1200.0 + case_index,
                    "pressure_pa": 2.0e6,
                    "initial_mass_fractions": dict(composition),
                    "ignition_delay": {
                        "definition": "first-1K-limited-linear-crossing-of-T-minus-T0-equals-400K",
                        "status": "scored-within-mechanism-thermo-intersection",
                        "time_s": 1.0e-4,
                    },
                    "samples": samples,
                }
            )
        return {
            "schema": self.tool.REFERENCE_SCHEMA,
            "status": self.tool.REFERENCE_STATUS,
            "generator_sha256": self.tool.sha256(TOOL),
            "authority": {
                "manifest_path": self.tool.INPUT_MANIFEST_REL.as_posix(),
                "composite_sha256": self.tool.INPUT_COMPOSITE_SHA256,
            },
            "runtime_identity": copy.deepcopy(self.tool.RUNTIME_IDENTITY),
            "mechanism": {},
            "rights_boundary": {
                "mechanism_code_license": "CC-BY-4.0-candidate",
                "experimental_training_rights_cleared": False,
                "experimental_holdout_rights_cleared": False,
                "experimental_values_ingested": False,
                "physical_fit_or_validation_claim": False,
            },
            "reactor": {},
            "matrices": [],
            "sample_times_s": list(self.tool.SAMPLE_TIMES),
            "species_order": list(self.tool.SPECIES),
            "observables": [],
            "cases": cases,
            "approvals": {
                "r1_r2_fit_approved": False,
                "physical_validation_approved": False,
                "phase_aware_chemistry_enabled": False,
                "fluent_execution_approved": False,
                "production_ready": False,
            },
        }

    def test_authority_composite_and_frozen_matrices(self) -> None:
        authority = self.tool.verify_authorities(WORKSPACE_ROOT)
        self.assertEqual(
            authority["composite_sha256"], self.tool.INPUT_COMPOSITE_SHA256
        )
        matrix_path, matrix = self.tool.validate_matrix(WORKSPACE_ROOT)
        self.assertEqual(
            matrix_path.relative_to(WORKSPACE_ROOT), self.tool.MATRIX_REL
        )
        self.assertEqual(len(matrix["sample_times_s"]), 54)
        self.assertTrue(
            set(case["id"] for case in matrix["cases"][:3]).isdisjoint(
                case["id"] for case in matrix["cases"][3:]
            )
        )

    def test_authority_rehash_cannot_replace_external_composite(self) -> None:
        with mock.patch.object(self.tool, "INPUT_COMPOSITE_SHA256", "0" * 64):
            with self.assertRaisesRegex(self.tool.ScreenError, "input-composite"):
                self.tool.verify_authorities(WORKSPACE_ROOT)

    def test_duplicate_keys_and_nonfinite_values_fail_closed(self) -> None:
        with self.assertRaisesRegex(self.tool.ScreenError, "duplicate-json-key:a"):
            self.tool.strict_pairs([("a", 1), ("a", 2)])
        with self.assertRaisesRegex(self.tool.ScreenError, "nonfinite-json-constant"):
            self.tool.reject_constant("NaN")

    def test_reference_scope_and_exit_order_are_closed(self) -> None:
        reference = self.synthetic_reference()
        self.tool.validate_reference(reference)
        round_tripped = json.loads(self.tool.canonical_bytes(reference).decode("ascii"))
        self.tool.validate_reference(round_tripped)

        approval = copy.deepcopy(reference)
        approval["approvals"]["production_ready"] = True
        with self.assertRaisesRegex(self.tool.ScreenError, "approval-enabled"):
            self.tool.validate_reference(approval)

        rights = copy.deepcopy(reference)
        rights["rights_boundary"]["experimental_training_rights_cleared"] = True
        with self.assertRaisesRegex(self.tool.ScreenError, "rights-boundary-weakened"):
            self.tool.validate_reference(rights)

        unknown = copy.deepcopy(reference)
        unknown["unreviewed"] = True
        with self.assertRaisesRegex(self.tool.ScreenError, "schema-open:reference"):
            self.tool.validate_reference(unknown)

        reordered = copy.deepcopy(reference)
        reordered["cases"][0]["samples"][3] = copy.deepcopy(
            reordered["cases"][0]["samples"][0]
        )
        reordered["cases"][0]["samples"][3]["time_s"] = self.tool.SAMPLE_TIMES[3]
        with self.assertRaisesRegex(self.tool.ScreenError, "eligible-sample-shape"):
            self.tool.validate_reference(reordered)

    def test_case_metrics_are_deterministic_and_include_six_species(self) -> None:
        detailed_y = {
            "C3H8": 0.1,
            "O2": 0.9,
            "CO": 0.0,
            "CO2": 0.0,
            "H2O": 0.0,
            "N2": 0.0,
        }
        case = {
            "id": "case",
            "partition": "protocol-training",
            "ignition_delay": {
                "status": "scored-within-mechanism-thermo-intersection",
                "time_s": 2.0e-4,
            },
            "samples": [
                {
                    "time_s": 0.0,
                    "status": "eligible-within-mechanism-thermo-intersection",
                    "temperature_k": 1000.0,
                    "heat_release_rate_w_m3": -1.0e8,
                    "mass_fractions": detailed_y,
                },
                {
                    "time_s": 1.0e-6,
                    "status": "eligible-within-mechanism-thermo-intersection",
                    "temperature_k": 1010.0,
                    "heat_release_rate_w_m3": -2.0e8,
                    "mass_fractions": detailed_y,
                },
            ],
        }
        rfam = [
            {
                "temperature_k": 1002.0,
                "qdot_w_m3": -2.0e8,
                "mass_fractions": [0.11, 0.89, 0.0, 0.0, 0.0, 0.0],
            },
            {
                "temperature_k": 1006.0,
                "qdot_w_m3": -1.0e8,
                "mass_fractions": [0.09, 0.91, 0.0, 0.0, 0.0, 0.0],
            },
        ]
        metrics = self.tool.compare_case(
            case,
            {"ignition_detected": "1", "ignition_delay_s": "0.0001"},
            rfam,
        )
        self.assertEqual(metrics["rfam_to_detailed_idt_ratio"], 0.5)
        self.assertEqual(metrics["temperature_rmse_k"], 3.16227766017)
        self.assertEqual(set(metrics["species"]), set(self.tool.SPECIES))
        self.assertEqual(metrics["species"]["C3H8"]["max_abs_mass_fraction"], 0.01)

    def test_publication_never_replaces_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.json"
            output.write_bytes(b"existing\n")
            with self.assertRaisesRegex(self.tool.ScreenError, "already-exists"):
                self.tool.publish_no_replace(output, b"replacement\n")
            self.assertEqual(output.read_bytes(), b"existing\n")


if __name__ == "__main__":
    unittest.main()
