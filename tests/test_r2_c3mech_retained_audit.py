import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))
PATH = TOOLS / "compare_r2_c3mech_retained_candidate.py"
SPEC = importlib.util.spec_from_file_location("retained_audit", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FakeCandidate:
    def __init__(self, species, reactions):
        self.species_names = species
        self.n_species = len(species)
        self.n_reactions = reactions


class RetainedAuditContractTests(unittest.TestCase):
    def test_exact_candidate_identity_is_required(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.yaml"
            path.write_bytes(b"candidate")
            expected = MODULE.sha256_file(path)
            candidate = FakeCandidate(list(MODULE.OBSERVABLE_SPECIES), 7)
            MODULE.validate_candidate_identity(
                path, expected, candidate, len(candidate.species_names), 7
            )
            with self.assertRaises(MODULE.CandidateError):
                MODULE.validate_candidate_identity(
                    path, "0" * 64, candidate, len(candidate.species_names), 7
                )

    def test_dimensions_and_observable_species_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.yaml"
            path.write_bytes(b"candidate")
            expected = MODULE.sha256_file(path)
            candidate = FakeCandidate(list(MODULE.OBSERVABLE_SPECIES), 7)
            with self.assertRaises(MODULE.CandidateError):
                MODULE.validate_candidate_identity(path, expected, candidate, 99, 7)
            missing = FakeCandidate(["C3H8", "O2"], 7)
            with self.assertRaises(MODULE.CandidateError):
                MODULE.validate_candidate_identity(path, expected, missing, 2, 7)


if __name__ == "__main__":
    unittest.main()
