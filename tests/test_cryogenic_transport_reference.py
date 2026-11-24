#!/usr/bin/env python3
"""Validate the frozen pure O2/propane transport reference matrix."""

from __future__ import annotations

import csv
import importlib.util
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "tools" / "generate_cryogenic_transport_reference.py"
REFERENCE = (
    ROOT / "tests" / "reference" / "cryogenic_transport_pure_coolprop72_v1.csv"
)

EXPECTED_FIELDS = [
    "schema_version",
    "case_id",
    "species",
    "coolprop_fluid",
    "cas",
    "backend",
    "coolprop_version",
    "coolprop_git_revision",
    "input_mode",
    "temperature_k",
    "pressure_pa",
    "quality_mole",
    "phase",
    "rho_mass_kg_m3",
    "viscosity_pa_s",
    "conductivity_w_m_k",
    "viscosity_reference",
    "conductivity_reference",
    "transport_json_sha256",
    "source_api",
    "production_ready",
    "notes",
]

TRANSPORT_IDENTITIES = {
    "O2": (
        "Lemmon-IJT-2004",
        "Lemmon-IJT-2004",
        "17d25c4a10f8ccd5ec69d14f0fb1d3327a687423b6bdeea77ee1f79bff449b8a",
    ),
    "C3H8": (
        "Vogel-JPCRD-1998",
        "Marsh-JCED-2002",
        "c3a43f7b18a507f67f04e08fc28f353c0ae5a8370503cb21485b2298782afbfb",
    ),
}

HAS_COOLPROP = importlib.util.find_spec("CoolProp") is not None


class CryogenicTransportReferenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with REFERENCE.open(newline="", encoding="ascii") as stream:
            reader = csv.DictReader(stream)
            cls.fields = reader.fieldnames
            cls.rows = list(reader)
        cls.by_id = {row["case_id"]: row for row in cls.rows}

    @unittest.skipUnless(HAS_COOLPROP, "pinned CoolProp is unavailable")
    def test_frozen_reference_is_current(self) -> None:
        result = subprocess.run(
            [sys.executable, str(GENERATOR), "--check"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("CRYOGENIC_TRANSPORT_REFERENCE_CURRENT", result.stdout)
        self.assertIn("production_ready=0", result.stdout)

    @unittest.skipUnless(HAS_COOLPROP, "pinned CoolProp is unavailable")
    def test_regeneration_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-transport-reference-") as tmp:
            output = Path(tmp) / "reference.csv"
            result = subprocess.run(
                [sys.executable, str(GENERATOR), "--output", str(output)],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(output.read_bytes(), REFERENCE.read_bytes())

    def test_schema_identity_and_scope(self) -> None:
        self.assertEqual(self.fields, EXPECTED_FIELDS)
        self.assertEqual(len(self.rows), 39)
        case_ids = [row["case_id"] for row in self.rows]
        self.assertEqual(case_ids, sorted(case_ids))
        self.assertEqual(len(set(case_ids)), len(case_ids))
        self.assertEqual(set(row["species"] for row in self.rows), {"O2", "C3H8"})
        for row in self.rows:
            self.assertEqual(
                row["schema_version"],
                "rfam.cryogenic-transport-pure-reference/v1",
            )
            self.assertEqual(row["backend"], "HEOS")
            self.assertEqual(row["coolprop_version"], "7.2.0")
            self.assertEqual(
                row["coolprop_git_revision"],
                "98b3523d5daa98454618d381d2ae53f7471d216b",
            )
            self.assertEqual(row["production_ready"], "0")
            self.assertIn("not-mixture-or-effective-two-phase-transport", row["notes"])
            self.assertIn("not-runtime-table", row["notes"])
            viscosity_ref, conductivity_ref, transport_sha = TRANSPORT_IDENTITIES[
                row["species"]
            ]
            self.assertEqual(row["viscosity_reference"], viscosity_ref)
            self.assertEqual(row["conductivity_reference"], conductivity_ref)
            self.assertEqual(row["transport_json_sha256"], transport_sha)
            for field in (
                "temperature_k",
                "pressure_pa",
                "rho_mass_kg_m3",
                "viscosity_pa_s",
                "conductivity_w_m_k",
            ):
                value = float(row[field])
                self.assertTrue(math.isfinite(value) and value > 0.0, (row, field))

    def test_anchor_values(self) -> None:
        oxygen = self.by_id["o2_pt_130k_4mpa"]
        self.assertEqual(oxygen["phase"], "liquid")
        self.assertTrue(
            math.isclose(
                float(oxygen["rho_mass_kg_m3"]),
                920.6598729327242,
                rel_tol=2.0e-14,
            )
        )
        self.assertTrue(
            math.isclose(
                float(oxygen["viscosity_pa_s"]),
                8.210537630511058e-5,
                rel_tol=2.0e-14,
            )
        )
        self.assertTrue(
            math.isclose(
                float(oxygen["conductivity_w_m_k"]),
                0.09635429027916971,
                rel_tol=2.0e-14,
            )
        )

        propane = self.by_id["c3h8_pt_450k_4mpa"]
        self.assertEqual(propane["phase"], "supercritical_gas")
        self.assertTrue(
            math.isclose(
                float(propane["rho_mass_kg_m3"]),
                57.45750719596648,
                rel_tol=2.0e-14,
            )
        )
        self.assertTrue(
            math.isclose(
                float(propane["viscosity_pa_s"]),
                1.3652057927061282e-5,
                rel_tol=2.0e-14,
            )
        )
        self.assertTrue(
            math.isclose(
                float(propane["conductivity_w_m_k"]),
                0.04342724491523844,
                rel_tol=2.0e-14,
            )
        )

    def test_saturation_endpoints_and_one_sided_states(self) -> None:
        for prefix, pressure_labels in (("o2", ("2mpa", "4mpa")), ("c3h8", ("1mpa", "4mpa"))):
            for label in pressure_labels:
                liquid = self.by_id[f"{prefix}_sat_liquid_{label}"]
                vapor = self.by_id[f"{prefix}_sat_vapor_{label}"]
                subcooled = self.by_id[f"{prefix}_sat_subcooled_0p25k_{label}"]
                superheated = self.by_id[f"{prefix}_sat_superheated_0p25k_{label}"]
                self.assertEqual(liquid["input_mode"], "PQ")
                self.assertEqual(vapor["input_mode"], "PQ")
                self.assertEqual(liquid["quality_mole"], "0.00000000000000000e+00")
                self.assertEqual(vapor["quality_mole"], "1.00000000000000000e+00")
                self.assertEqual(liquid["temperature_k"], vapor["temperature_k"])
                self.assertEqual(liquid["pressure_pa"], vapor["pressure_pa"])
                self.assertGreater(
                    float(liquid["rho_mass_kg_m3"]),
                    float(vapor["rho_mass_kg_m3"]),
                )
                self.assertIn(subcooled["phase"], {"liquid", "supercritical_liquid"})
                self.assertIn(superheated["phase"], {"gas", "supercritical_gas"})
                self.assertEqual(subcooled["quality_mole"], "")
                self.assertEqual(superheated["quality_mole"], "")

    def test_critical_point_is_bracketed_not_queried_exactly(self) -> None:
        for prefix in ("o2", "c3h8"):
            below = self.by_id[f"{prefix}_near_critical_below"]
            above = self.by_id[f"{prefix}_near_critical_above"]
            self.assertEqual(below["input_mode"], "PT")
            self.assertEqual(above["input_mode"], "PT")
            self.assertLess(float(below["temperature_k"]), float(above["temperature_k"]))
            self.assertEqual(below["pressure_pa"], above["pressure_pa"])
            self.assertNotEqual(below["phase"], "twophase")
            self.assertNotEqual(above["phase"], "twophase")


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        CryogenicTransportReferenceTest
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(
            "CRYOGENIC_TRANSPORT_REFERENCE_TEST_PASS "
            f"tests={result.testsRun} rows={len(CryogenicTransportReferenceTest.rows)} "
            "production_ready=0"
        )
        raise SystemExit(0)
    print("CRYOGENIC_TRANSPORT_REFERENCE_TEST_FAIL")
    raise SystemExit(1)
