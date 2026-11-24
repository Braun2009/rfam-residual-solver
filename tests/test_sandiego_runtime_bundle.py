import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parent.parent
BUNDLE = ROOT / "tests/reference/sandiego_runtime_bundle_v1.json"
SOURCE = ROOT / "tests/reference/sandiego.yaml"


class SandiegoRuntimeBundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(BUNDLE.read_text(encoding="ascii"))

    def test_identity_and_counts(self):
        document = self.document
        self.assertEqual(document["schema"], "rfam.sandiego-runtime-bundle/v1")
        self.assertEqual(document["species_count"], 59)
        self.assertEqual(document["reaction_count"], 271)
        self.assertEqual(len(document["species"]), 59)
        self.assertEqual(len(document["reactions"]), 271)
        digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        self.assertEqual(document["source"]["sha256"], digest)

    def test_species_thermo_is_complete(self):
        names = [row["name"] for row in self.document["species"]]
        self.assertEqual(len(names), len(set(names)))
        self.assertIn("C3H8", names)
        self.assertIn("OH", names)
        for row in self.document["species"]:
            self.assertGreater(row["molecular_weight_kg_per_mol"], 0.0)
            self.assertEqual(row["thermo"]["model"], "NASA7")
            self.assertEqual(len(row["thermo"]["temperature-ranges"]), 3)
            self.assertEqual(len(row["thermo"]["data"]), 2)
            self.assertTrue(all(len(coefficients) == 7 for coefficients in row["thermo"]["data"]))

    def test_reaction_semantics_cover_pressure_dependence(self):
        counts = self.document["reaction_type_counts"]
        self.assertEqual(counts["Arrhenius"], 233)
        self.assertEqual(counts["three-body-Arrhenius"], 15)
        self.assertEqual(counts["falloff-Troe"], 23)
        self.assertEqual(sum(counts.values()), 271)
        self.assertEqual(
            [row["index"] for row in self.document["reactions"]],
            list(range(271)),
        )
        self.assertTrue(any(row["duplicate"] for row in self.document["reactions"]))
        self.assertTrue(all(row["reactants"] and row["products"] for row in self.document["reactions"]))

    def test_rate_references_are_dimensionally_complete(self):
        references = self.document["rate_reference_states"]
        self.assertEqual([row["id"] for row in references], ["lean-hot", "stoich-hot", "rich-hot"])
        for row in references:
            self.assertEqual(len(row["forward_rates_of_progress_kmol_m3_s"]), 271)
            self.assertEqual(len(row["reverse_rates_of_progress_kmol_m3_s"]), 271)
            self.assertEqual(len(row["net_rates_of_progress_kmol_m3_s"]), 271)
            self.assertEqual(len(row["net_production_rates_kmol_m3_s"]), 59)
            self.assertGreater(row["pressure_Pa"], 0.0)
            self.assertGreater(row["temperature_K"], 0.0)


if __name__ == "__main__":
    unittest.main()
