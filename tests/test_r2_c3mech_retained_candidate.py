import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = ROOT / "tools" / "build_r2_c3mech_retained_candidate.py"
SPEC = importlib.util.spec_from_file_location("retained_candidate", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FakeThirdBody:
    def __init__(self, efficiencies):
        self.efficiencies = efficiencies


class FakeReaction:
    def __init__(
        self,
        equation,
        reactants,
        products,
        *,
        duplicate=False,
        orders=None,
        efficiencies=None,
    ):
        self.equation = equation
        self.reactants = reactants
        self.products = products
        self.duplicate = duplicate
        self.orders = orders or {}
        self.third_body = (
            FakeThirdBody(efficiencies) if efficiencies is not None else None
        )


class RetainedCandidateContractTests(unittest.TestCase):
    def test_duplicate_json_key_rejected(self):
        with self.assertRaises(MODULE.CandidateError):
            json.loads('{"a":1,"a":2}', object_pairs_hook=MODULE.reject_duplicate_keys)

    def test_duplicate_reaction_group_and_dependencies_are_closed(self):
        reactions = [
            FakeReaction("A + B <=> C", {"A": 1, "B": 1}, {"C": 1}, duplicate=True),
            FakeReaction("A + B <=> C", {"A": 1, "B": 1}, {"C": 1}, duplicate=True),
            FakeReaction(
                "C + M <=> D + M",
                {"C": 1},
                {"D": 1},
                orders={"Q": 0.5},
                efficiencies={"N2": 2.0, "AR": 0.7},
            ),
        ]
        indices, species = MODULE.close_selected_reactions(reactions, [0, 2])
        self.assertEqual(indices, [0, 1, 2])
        self.assertTrue({"A", "B", "C", "D", "Q", "N2", "AR"}.issubset(species))

    def test_threshold_contract_requires_descending_zero_fallback(self):
        self.assertEqual(MODULE.validate_thresholds((0.1, 0.01, 0.0)), (0.1, 0.01, 0.0))
        for values in ((0.1, 0.01), (0.01, 0.1, 0.0), (0.1, -1.0, 0.0)):
            with self.subTest(values=values), self.assertRaises(MODULE.CandidateError):
                MODULE.validate_thresholds(values)

    def test_case_metrics_exposes_all_admission_components(self):
        samples = []
        for index, rise in enumerate(MODULE.TEMPERATURE_RISE_NODES):
            samples.append(
                {
                    "temperature_rise": rise,
                    "time": float(index),
                    "temperature": 1200.0 + rise,
                    "heat_release_rate": float(index + 1),
                    "mass_fractions": {name: 0.01 * index for name in MODULE.OBSERVABLE_SPECIES},
                }
            )
        trajectory = MODULE.Trajectory(True, "", 7.0, samples, [], 10, 1500.0)
        tolerances = MODULE.Tolerances(1e-12, 1e-12, 1e-12, 1e-12)
        result = MODULE.case_metrics(trajectory, trajectory, tolerances)
        self.assertTrue(result["pass"])
        self.assertEqual(result["idt_relative"], 0.0)
        self.assertEqual(result["inverse_temperature_time_relative"], 0.0)
        self.assertEqual(result["heat_release_nrmse"], 0.0)
        self.assertEqual(result["species_mass_fraction_absolute_maximum"], 0.0)

    def test_semantic_digest_is_order_and_value_sensitive(self):
        class Input:
            def __init__(self, data):
                self.input_data = data

        left = MODULE.semantic_digest([Input({"rate": {"A": 1.0, "b": 2.0}})])
        reordered = MODULE.semantic_digest([Input({"rate": {"b": 2.0, "A": 1.0}})])
        changed = MODULE.semantic_digest([Input({"rate": {"A": 1.0, "b": 3.0}})])
        self.assertEqual(left, reordered)
        self.assertNotEqual(left, changed)

    def test_semantic_comparison_ignores_notes_but_rejects_rate_drift(self):
        class Input:
            def __init__(self, data):
                self.input_data = data

        left = [Input({"note": "source ", "rate": {"A": 1.0, "b": 2.0}})]
        roundtrip = [Input({"note": "source", "rate": {"A": 1.0 + 1e-15, "b": 2.0}})]
        self.assertLess(MODULE.compare_semantic_items(left, roundtrip, "reaction"), 5e-14)
        changed = [Input({"rate": {"A": 1.01, "b": 2.0}})]
        with self.assertRaises(MODULE.CandidateError):
            MODULE.compare_semantic_items(left, changed, "reaction")

    def test_yaml_canonicalization_removes_generated_date(self):
        try:
            import yaml  # noqa: F401
        except ImportError:
            self.skipTest("PyYAML unavailable")
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.yaml"
            path.write_text("generator: YamlWriter\ndate: now\nphases: []\n", encoding="ascii")
            first, _ = MODULE.canonicalize_cantera_yaml(path)
            path.write_text("generator: YamlWriter\ndate: later\nphases: []\n", encoding="ascii")
            second, _ = MODULE.canonicalize_cantera_yaml(path)
            self.assertEqual(first, second)
            self.assertIn(b"date: '1970-01-01T00:00:00Z'", first)


if __name__ == "__main__":
    unittest.main()
