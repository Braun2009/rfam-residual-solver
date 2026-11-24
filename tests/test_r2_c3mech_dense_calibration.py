#!/usr/bin/env python3
"""Focused tests for the dense C3Mech train/holdout calibration route."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPOSITORY_ROOT.parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(relative: str, name: str):
    path = REPOSITORY_ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError("module-load-spec-failed:" + relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class DenseCalibrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.generator = load_module(
            "tools/generate_r2_c3mech_dense_reference.py", "test_bound_dense_generator"
        )
        cls.comparator = load_module(
            "tools/compare_r2_c3mech_dense_reference.py", "test_bound_dense_comparator"
        )

    def test_frozen_identities(self) -> None:
        expected = {
            REPOSITORY_ROOT / "tools/generate_r2_c3mech_dense_reference.py":
                "7b417d7d77f1a99690e560658dbf54cee890becf2bf5afc454535ab7e23da386",
            REPOSITORY_ROOT / "tools/compare_r2_c3mech_dense_reference.py":
                "a1c2798038655132657855f450a4ce0bbdd2f2f14abee15394d86c0d0b141a41",
            REPOSITORY_ROOT / "tests/reference/r2_c3mech_dense_training_cases_v1.json":
                "de56596277f253b748cf4fda168e64489ca47436f9045b1f7d22ad0f46649cdf",
            REPOSITORY_ROOT / "tests/reference/r2_c3mech_dense_holdout_cases_v1.json":
                "92496aeb3fdd38c6d5006919249ebf5698058b93595515c92a5edc204a2d7a6e",
            REPOSITORY_ROOT / "tests/reference/r2_c3mech_coarse_selection_v1.json":
                "2cd3842e2dac60cd82900abad760ac46a720ec051c3a5d32f50fccb6b5d6a562",
            WORKSPACE_ROOT / "tmp/a408-r2-c3mech-44i0-dense-training-v1.json":
                "ae629e1501531d803eec8a74db37dc1e840c5d6679528918e9f8f2edd3862247",
            WORKSPACE_ROOT / "tmp/a408-r2-c3mech-44i0-dense-holdout-v1.json":
                "e1cb70bee2ea084e7d1ed79ebe4481bc87c7c98b1df4a03e0b4fdd3cfa24b65e",
            WORKSPACE_ROOT / "tmp/a408-r2-c3mech-dense-comparison-v1-final.json":
                "b6216844128c36de8971fe2fe1c6b7e50ecc8106aa750d5fffe3c1f5738b041b",
        }
        for path, digest in expected.items():
            with self.subTest(path=path.name):
                self.assertTrue(path.is_file())
                self.assertEqual(sha256(path), digest)

    def test_predeclared_manifests_are_distinct_and_dense(self) -> None:
        training_path, training = self.generator.validate_manifest(WORKSPACE_ROOT, "training")
        holdout_path, holdout = self.generator.validate_manifest(WORKSPACE_ROOT, "holdout")
        self.assertNotEqual(training_path, holdout_path)
        self.assertEqual(training["partition"], "protocol-training")
        self.assertEqual(holdout["partition"], "protocol-holdout")
        self.assertEqual(len(training["sample_times_s"]), 42)
        self.assertEqual(training["sample_times_s"], holdout["sample_times_s"])
        self.assertTrue(set(case["id"] for case in training["cases"]).isdisjoint(
            case["id"] for case in holdout["cases"]
        ))

    def test_generated_references_have_longer_eligible_prefixes(self) -> None:
        expected = {
            "training_reference": {
                "train_phi0p8_diluted": 27,
                "train_phi1p0": 23,
                "train_phi1p2": 20,
            },
            "holdout_reference": {
                "holdout_phi0p6_diluted": 28,
                "holdout_phi1p1": 19,
                "holdout_phi1p4": 20,
            },
        }
        for artifact_id, counts in expected.items():
            path = self.comparator.require_file(WORKSPACE_ROOT, artifact_id)
            document = self.comparator.load_json(path, artifact_id)
            partition = "protocol-training" if artifact_id.startswith("training") else "protocol-holdout"
            self.comparator.reference_view(document, partition)
            self.assertTrue(all(value is False for value in document["approvals"].values()))
            for case in document["cases"]:
                eligible = [
                    sample for sample in case["samples"]
                    if sample["status"] == "eligible-within-mechanism-thermo-intersection"
                ]
                self.assertEqual(len(eligible), counts[case["id"]])
                self.assertGreaterEqual(len(eligible), 19)
                self.assertEqual(
                    case["thermo_envelope_exit"]["status"],
                    "thermo-envelope-exit-detected-no-endpoint-state-published",
                )

    def test_dense_comparison_records_training_before_holdout(self) -> None:
        path = WORKSPACE_ROOT / "tmp/a408-r2-c3mech-dense-comparison-v1-final.json"
        document = self.comparator.load_json(path, "dense-comparison")
        self.assertEqual(document["schema"], "rfam.r2-c3mech-dense-comparison/v1")
        self.assertEqual(
            document["access_order"],
            ["selection", "training-reference", "holdout-reference"],
        )
        self.assertFalse(document["selection"]["holdout_influenced_selection"])
        self.assertFalse(document["selection"]["parameter_vector_publishable"])
        self.assertAlmostEqual(document["training"]["objective"], 6.38898327395609, places=13)
        self.assertAlmostEqual(document["holdout"]["objective"], 1.07690588375512, places=13)
        self.assertGreater(
            document["training"]["objective"],
            document["selection"]["coarse_training_objective"] * 3.0,
        )
        self.assertTrue(all(value is False for value in document["approvals"].values()))

    def test_reference_policy_and_source_drift_fail_closed(self) -> None:
        path = self.comparator.require_file(WORKSPACE_ROOT, "training_reference")
        document = self.comparator.load_json(path, "training-reference")
        adaptive = copy.deepcopy(document)
        adaptive["sampling_policy"]["candidate_or_reference_adaptive"] = True
        with self.assertRaisesRegex(Exception, "dense-sampling-policy-mismatch"):
            self.comparator.reference_view(adaptive, "protocol-training")
        rebound = copy.deepcopy(document)
        rebound["source"]["dense_generator_sha256"] = "0" * 64
        with self.assertRaisesRegex(Exception, "dense-source-authority-mismatch"):
            self.comparator.reference_view(rebound, "protocol-training")

    def test_selection_unknown_field_fails_even_with_matching_file_hash(self) -> None:
        canonical_path = REPOSITORY_ROOT / "tests/reference/r2_c3mech_coarse_selection_v1.json"
        document = json.loads(canonical_path.read_text(encoding="utf-8"))
        document["unreviewed"] = True
        with tempfile.TemporaryDirectory(dir=WORKSPACE_ROOT / "tmp") as directory:
            root = Path(directory)
            path = root / "selection.json"
            raw = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode("utf-8")
            path.write_bytes(raw)
            original = self.comparator.ARTIFACTS["selection"]
            self.comparator.ARTIFACTS["selection"] = (
                Path("selection.json"), hashlib.sha256(raw).hexdigest(), len(raw)
            )
            try:
                with self.assertRaisesRegex(Exception, "selection-schema-open"):
                    self.comparator.load_selection(root)
            finally:
                self.comparator.ARTIFACTS["selection"] = original

    def test_duplicate_keys_and_preexisting_outputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(Exception, "duplicate-json-key"):
            self.generator.strict_pairs([("a", 1), ("a", 2)])
        with self.assertRaisesRegex(Exception, "duplicate-json-key"):
            self.comparator.strict_pairs([("a", 1), ("a", 2)])
        with tempfile.TemporaryDirectory(dir=WORKSPACE_ROOT / "tmp") as directory:
            output = Path(directory) / "sentinel.json"
            output.write_bytes(b"preserve-me\n")
            before = output.read_bytes()
            with self.assertRaisesRegex(Exception, "output-path-already-exists"):
                self.comparator.publish_no_replace(output, b"replacement\n")
            self.assertEqual(output.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
