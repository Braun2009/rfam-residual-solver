#!/usr/bin/env python3
"""Focused fail-close tests for the offline C3Mech candidate screen."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
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
TOOL = REPO / "tools" / "validate_r2_c3mech_candidate_screen.py"
CONTRACT = REPO / "tests" / "reference" / "r2_c3mech_candidate_screen_contract_v2.json"
BASE_TEXT = CONTRACT.read_text(encoding="utf-8")
BASE = json.loads(BASE_TEXT)
GATE_SPEC = importlib.util.spec_from_file_location("r2_c3mech_gate", TOOL)
if GATE_SPEC is None or GATE_SPEC.loader is None:
    raise RuntimeError("unable-to-load-r2-c3mech-gate")
GATE = importlib.util.module_from_spec(GATE_SPEC)
GATE_SPEC.loader.exec_module(GATE)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(document: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    return next(item for item in document["artifacts"] if item["id"] == artifact_id)


class R2C3MechCandidateScreenTests(unittest.TestCase):
    def run_gate(
        self,
        contract: Path = CONTRACT,
        root: Path = WORKSPACE,
        *extra: str,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(TOOL),
                "--workspace-root",
                str(root),
                "--contract",
                str(contract),
                *extra,
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    def require_failure(
        self, result: subprocess.CompletedProcess[str], needle: str
    ) -> None:
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("R2_C3MECH_CANDIDATE_SCREEN_GATE_FAIL", result.stderr)
        self.assertIn(needle, result.stderr)
        self.assertNotIn("R2_C3MECH_CANDIDATE_SCREEN_GATE_BLOCKED", result.stdout + result.stderr)
        self.assertNotIn("Traceback", result.stdout + result.stderr)

    def mutated_contract(
        self, mutation: Callable[[dict[str, Any]], None]
    ) -> subprocess.CompletedProcess[str]:
        document = copy.deepcopy(BASE)
        mutation(document)
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-contract-") as directory:
            path = Path(directory) / "contract.json"
            path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            return self.run_gate(path)

    def raw_contract(self, text: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-contract-") as directory:
            path = Path(directory) / "contract.json"
            path.write_text(text, encoding="utf-8")
            return self.run_gate(path)

    def fixture(self, root: Path) -> Path:
        for item in BASE["artifacts"]:
            source = WORKSPACE / item["path"]
            target = root / item["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        contract = root / "realFluidAllMach/tests/reference/r2_c3mech_candidate_screen_contract_v2.json"
        contract.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(CONTRACT, contract)
        return contract

    def canonical_paths(self) -> dict[str, Path]:
        return {
            item["id"]: WORKSPACE / item["path"]
            for item in BASE["artifacts"]
        }

    def test_canonical_gate_is_blocked(self) -> None:
        result = self.run_gate()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_C3MECH_CANDIDATE_SCREEN_GATE_BLOCKED checks=25", result.stdout)

    def test_fast_screen_regeneration_is_byte_identical(self) -> None:
        result = self.run_gate(CONTRACT, WORKSPACE, "--verify-screen-regeneration")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("R2_C3MECH_CANDIDATE_SCREEN_GATE_BLOCKED checks=26", result.stdout)

    def test_canonical_trajectory_metrics_are_null_without_noninitial_pairs(self) -> None:
        output = json.loads(
            self.canonical_paths()["screen_output"].read_text(encoding="utf-8")
        )
        self.assertEqual(
            output["full_fixed_protocol"]["status"],
            "OUT_OF_MECHANISM_THERMO_ENVELOPE",
        )
        self.assertEqual(output["full_fixed_protocol"]["case_count"], 12)
        results = output["prefix_diagnostic"]["results"]
        for model_id, model in results.items():
            for partition in ("training", "holdout"):
                with self.subTest(model=model_id, partition=partition):
                    aggregate = model[partition + "_aggregate"]
                    self.assertEqual(aggregate["noninitial_eligible_sample_count"], 0)
                    self.assertFalse(aggregate["trajectory_error_metrics_available"])
                    self.assertIsNone(aggregate["maximum_temperature_absolute_error_k"])
                    self.assertIsNone(aggregate["maximum_major_species_l1"])
                    for case in model[partition]:
                        self.assertEqual(case["paired_noninitial_sample_count"], 0)
                        self.assertIsNone(case["maximum_temperature_absolute_error_k"])
                        self.assertIsNone(case["maximum_major_species_l1"])
                        self.assertIsNone(case["last_eligible_temperature_delta_k"])

    def test_numeric_zero_cannot_replace_unavailable_trajectory_metrics(self) -> None:
        mutations = [
            ("case-temperature", lambda document: document["prefix_diagnostic"]["results"]["44I0"]["training"][0].__setitem__("maximum_temperature_absolute_error_k", 0.0)),
            ("case-species", lambda document: document["prefix_diagnostic"]["results"]["44I0"]["training"][0].__setitem__("maximum_major_species_l1", 0.0)),
            ("case-last-delta", lambda document: document["prefix_diagnostic"]["results"]["44I0"]["training"][0].__setitem__("last_eligible_temperature_delta_k", 0.0)),
            ("aggregate-temperature", lambda document: document["prefix_diagnostic"]["results"]["44I0"]["training_aggregate"].__setitem__("maximum_temperature_absolute_error_k", 0.0)),
            ("aggregate-species", lambda document: document["prefix_diagnostic"]["results"]["44I0"]["training_aggregate"].__setitem__("maximum_major_species_l1", 0.0)),
        ]
        for label, mutation in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory(
                prefix="rfam-c3mech-semantic-"
            ) as directory:
                paths = self.canonical_paths()
                document = json.loads(paths["screen_output"].read_text(encoding="utf-8"))
                mutation(document)
                output = Path(directory) / "screen.json"
                output.write_text(json.dumps(document), encoding="utf-8")
                paths["screen_output"] = output
                with self.assertRaisesRegex(GATE.GateError, "screen-output-"):
                    GATE.verify_screen_output(BASE, paths)

    def test_wheel_installed_tree_equivalence_claim_is_rejected_semantically(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-semantic-") as directory:
            paths = self.canonical_paths()
            document = json.loads(paths["screen_plan"].read_text(encoding="utf-8"))
            document["runtime_identity"][
                "wheel_to_installed_tree_equivalence_proven"
            ] = True
            plan = Path(directory) / "plan.json"
            plan.write_text(json.dumps(document), encoding="utf-8")
            paths["screen_plan"] = plan
            with self.assertRaisesRegex(
                GATE.GateError, "plan-wheel-installed-tree-equivalence-enabled"
            ):
                GATE.verify_plan(BASE, paths)

    def test_crossing_evidence_is_bound_and_nonpublishable(self) -> None:
        paths = self.canonical_paths()
        reference = json.loads(
            paths["reference_3r00_training"].read_text(encoding="utf-8")
        )
        for case in reference["cases"]:
            crossing = case["thermo_envelope_exit"]["first_crossing_estimate"]
            self.assertEqual(
                set(crossing),
                {
                    "method", "pressure_pa", "publishable_chemistry_state",
                    "temperature_k", "time_s",
                },
            )
            self.assertFalse(crossing["publishable_chemistry_state"])

        mutations = [
            ("publishable", lambda crossing: crossing.__setitem__("publishable_chemistry_state", True)),
            ("extra-state", lambda crossing: crossing.__setitem__("density_kg_m3", 1.0)),
        ]
        for label, mutation in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory(
                prefix="rfam-c3mech-semantic-"
            ) as directory:
                mutated = copy.deepcopy(reference)
                crossing = mutated["cases"][0]["thermo_envelope_exit"][
                    "first_crossing_estimate"
                ]
                mutated["cases"][0]["samples"][
                    next(
                        index for index, sample in enumerate(mutated["cases"][0]["samples"])
                        if sample["status"] == "outside-thermo-envelope-unscored"
                    )
                ]["exit_detection"]["first_crossing_estimate"] = crossing
                mutation(crossing)
                candidate = Path(directory) / "reference.json"
                candidate.write_text(json.dumps(mutated), encoding="utf-8")
                semantic_paths = self.canonical_paths()
                semantic_paths["reference_3r00_training"] = candidate
                with self.assertRaisesRegex(GATE.GateError, "reference-crossing-"):
                    GATE.verify_references(BASE, semantic_paths)

    def test_every_approval_remains_false(self) -> None:
        for key in BASE["approvals"]:
            with self.subTest(key=key):
                self.require_failure(
                    self.mutated_contract(
                        lambda document, key=key: document["approvals"].__setitem__(key, True)
                    ),
                    "contract-",
                )

    def test_unknown_root_authorization_is_rejected(self) -> None:
        self.require_failure(
            self.mutated_contract(
                lambda document: document.__setitem__("production_authorization", True)
            ),
            "contract-",
        )

    def test_duplicate_nested_approval_is_rejected(self) -> None:
        needle = '    "production_ready": false,'
        replacement = '    "production_ready": true,\n' + needle
        self.assertEqual(BASE_TEXT.count(needle), 1)
        self.require_failure(
            self.raw_contract(BASE_TEXT.replace(needle, replacement, 1)),
            "contract-",
        )

    def test_self_consistent_yaml_rehash_is_rejected_by_frozen_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            mechanism_item = artifact(BASE, "mechanism_44i0")
            mechanism = root / mechanism_item["path"]
            mechanism.write_bytes(mechanism.read_bytes() + b"\n")
            document = json.loads(contract_path.read_text(encoding="utf-8"))
            item = artifact(document, "mechanism_44i0")
            item["sha256"] = sha256(mechanism)
            item["size"] = mechanism.stat().st_size
            contract_path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            self.require_failure(
                self.run_gate(contract_path, root), "contract-"
            )

    def test_live_wheel_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            wheel = root / artifact(BASE, "cantera_wheel")["path"]
            wheel.write_bytes(wheel.read_bytes() + b"\n")
            self.require_failure(
                self.run_gate(contract_path, root),
                "artifact-size-mismatch:cantera_wheel",
            )

    def test_self_consistent_wheel_rehash_is_rejected_by_frozen_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            wheel = root / artifact(BASE, "cantera_wheel")["path"]
            wheel.write_bytes(wheel.read_bytes() + b"\n")
            document = json.loads(contract_path.read_text(encoding="utf-8"))
            wheel_item = artifact(document, "cantera_wheel")
            wheel_item["sha256"] = sha256(wheel)
            wheel_item["size"] = wheel.stat().st_size
            document["runtime_identity"]["cantera_wheel_sha256"] = sha256(wheel)
            document["runtime_identity"]["cantera_wheel_size"] = wheel.stat().st_size
            contract_path.write_text(
                json.dumps(document, indent=2) + "\n", encoding="utf-8"
            )
            self.require_failure(self.run_gate(contract_path, root), "contract-")

    def test_live_screen_output_tamper_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            output = root / artifact(BASE, "screen_output")["path"]
            output.write_bytes(output.read_bytes() + b"\n")
            self.require_failure(
                self.run_gate(contract_path, root), "artifact-size-mismatch:screen_output"
            )

    def test_missing_runner_is_rejected_before_regeneration(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            (root / artifact(BASE, "rfam_runner_executable")["path"]).unlink()
            self.require_failure(
                self.run_gate(contract_path, root, "--verify-screen-regeneration"),
                "artifact-missing:rfam_runner_executable",
            )

    def test_generator_substitution_is_rejected_before_launch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            generator = root / artifact(BASE, "reference_generator")["path"]
            generator.write_text("raise SystemExit('substitute')\n", encoding="utf-8")
            self.require_failure(
                self.run_gate(
                    contract_path,
                    root,
                    "--verify-reference-regeneration",
                    "3R00-training",
                ),
                "artifact-size-mismatch:reference_generator",
            )

    def test_role_lock_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            plan = root / artifact(BASE, "screen_plan")["path"]
            document = json.loads(plan.read_text(encoding="utf-8"))
            document["role_lock"]["selected_model_id"] = "44I0"
            plan.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            self.require_failure(
                self.run_gate(contract_path, root), "artifact-size-mismatch:screen_plan"
            )

    def test_post_exit_payload_reintroduction_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="rfam-c3mech-fixture-") as directory:
            root = Path(directory)
            contract_path = self.fixture(root)
            reference = root / artifact(BASE, "reference_3r00_holdout")["path"]
            document = json.loads(reference.read_text(encoding="utf-8"))
            document["cases"][0]["samples"][-1]["temperature_k"] = 3000.0
            reference.write_text(
                json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )
            self.require_failure(
                self.run_gate(contract_path, root), "artifact-size-mismatch:reference_3r00_holdout"
            )


if __name__ == "__main__":
    unittest.main()
