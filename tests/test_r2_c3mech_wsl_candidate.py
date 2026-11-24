#!/usr/bin/env python3
"""Focused tests for the expanded-grid C3Mech candidate evaluation."""

from __future__ import annotations

import copy
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPOSITORY_ROOT.parent
TOOL = REPOSITORY_ROOT / "tools/evaluate_r2_c3mech_wsl_candidate.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("rfam_r2_c3mech_candidate", TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError("tool-import-spec-failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class CandidateEvaluationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tool = load_tool()

    def test_frozen_candidate_reference_and_runner(self) -> None:
        for artifact_id in ("candidate-result", "reference", "runner"):
            with self.subTest(artifact_id=artifact_id):
                self.assertTrue(self.tool.require_file(WORKSPACE_ROOT, artifact_id).is_file())
        self.tool.base.verify_authorities(WORKSPACE_ROOT)

    def test_candidate_selection_and_nonproduction_boundary(self) -> None:
        path = self.tool.require_file(WORKSPACE_ROOT, "candidate-result")
        candidate = self.tool.load_json(path, "candidate")
        self.tool.validate_candidate(candidate)
        holdout = copy.deepcopy(candidate)
        holdout["runner_accounting"]["holdout_access_attempts_during_selection"] = 1
        with self.assertRaisesRegex(self.tool.CandidateError, "holdout-accounting"):
            self.tool.validate_candidate(holdout)
        production = copy.deepcopy(candidate)
        production["approvals"]["production_ready"] = True
        with self.assertRaisesRegex(self.tool.CandidateError, "approval-mismatch"):
            self.tool.validate_candidate(production)

    def test_candidate_yaml_is_separate_and_contains_exact_parameters(self) -> None:
        payload = self.tool.candidate_yaml_bytes().decode("ascii")
        self.assertIn("does not replace RFAM defaults", payload)
        self.assertIn("A: " + format(self.tool.PARAMETERS[0], ".17g"), payload)
        self.assertIn("Ea: " + format(self.tool.PARAMETERS[5], ".17g"), payload)
        self.assertIn("O2: " + format(self.tool.PARAMETERS[3], ".17g"), payload)
        self.assertNotIn("production-ready: true", payload)

    def test_duplicate_keys_and_nonfinite_values_fail_closed(self) -> None:
        with self.assertRaisesRegex(self.tool.CandidateError, "duplicate-json-key"):
            self.tool.strict_pairs([("a", 1), ("a", 2)])
        with self.assertRaisesRegex(self.tool.CandidateError, "nonfinite-json"):
            self.tool.reject_constant("Infinity")

    def test_publication_does_not_replace_existing_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "candidate.yaml"
            output.write_bytes(b"existing\n")
            with self.assertRaisesRegex(self.tool.CandidateError, "already-exists"):
                self.tool.publish_no_replace(output, b"replacement\n")
            self.assertEqual(output.read_bytes(), b"existing\n")


if __name__ == "__main__":
    unittest.main()
