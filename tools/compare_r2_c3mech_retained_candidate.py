#!/usr/bin/env python3
"""Compare one already-frozen retained mechanism against full C3Mech."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Optional, Sequence

from build_r2_c3mech_retained_candidate import (
    CandidateError,
    OBSERVABLE_SPECIES,
    Tolerances,
    canonical_bytes,
    case_metrics,
    fresh_output,
    load_cases,
    rate_type_counts,
    run_trajectory,
    sha256_file,
)


def validate_candidate_identity(
    candidate_path: Path,
    expected_sha256: str,
    candidate: Any,
    expected_species: int,
    expected_reactions: int,
) -> None:
    if len(expected_sha256) != 64 or any(character not in "0123456789abcdef" for character in expected_sha256):
        raise CandidateError("candidate-expected-sha256-format")
    if sha256_file(candidate_path) != expected_sha256:
        raise CandidateError("candidate-sha256-mismatch")
    if candidate.n_species != expected_species:
        raise CandidateError("candidate-species-count-mismatch")
    if candidate.n_reactions != expected_reactions:
        raise CandidateError("candidate-reaction-count-mismatch")
    missing = [name for name in OBSERVABLE_SPECIES if name not in candidate.species_names]
    if missing:
        raise CandidateError("candidate-observable-species-missing:" + ",".join(missing))


def write_fresh(path: Path, raw: bytes) -> None:
    try:
        with path.open("xb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise CandidateError("result-output-already-exists") from error


def execute(arguments: argparse.Namespace) -> dict[str, Any]:
    import cantera as ct

    source_path = arguments.mechanism.resolve(strict=True)
    candidate_path = arguments.candidate.resolve(strict=True)
    cases_path = arguments.cases.resolve(strict=True)
    source = ct.Solution(str(source_path))
    candidate = ct.Solution(str(candidate_path))
    validate_candidate_identity(
        candidate_path,
        arguments.expected_candidate_sha256,
        candidate,
        arguments.expected_species,
        arguments.expected_reactions,
    )
    cases = load_cases(cases_path, "retained-candidate-audit-v1")
    tolerances = Tolerances(
        arguments.idt_relative_tolerance,
        arguments.inverse_temperature_time_tolerance,
        arguments.heat_release_nrmse_tolerance,
        arguments.species_absolute_tolerance,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in tolerances.__dict__.values()):
        raise CandidateError("tolerance-range")

    rows: dict[str, Any] = {}
    passed = True
    for case in cases:
        reference = run_trajectory(
            lambda: ct.Solution(str(source_path)), case, False, arguments.maximum_time
        )
        if not reference.complete:
            raise CandidateError("audit-reference-incomplete:" + case.case_id + ":" + reference.failure)
        trial = run_trajectory(
            lambda: ct.Solution(str(candidate_path)), case, False, arguments.maximum_time
        )
        metrics = case_metrics(reference, trial, tolerances)
        rows[case.case_id] = metrics
        passed = passed and bool(metrics["pass"])

    return {
        "schema": "rfam.r2-c3mech-retained-audit/v1",
        "status": "internal-detailed-mechanism-surrogate-audit",
        "source": {
            "mechanism_path": str(source_path),
            "mechanism_sha256": sha256_file(source_path),
            "candidate_path": str(candidate_path),
            "candidate_sha256": sha256_file(candidate_path),
            "cases_path": str(cases_path),
            "cases_sha256": sha256_file(cases_path),
            "cantera_version": ct.__version__,
            "python_version": sys.version.split()[0],
        },
        "candidate": {
            "species_count": candidate.n_species,
            "reaction_count": candidate.n_reactions,
            "reaction_types": rate_type_counts(candidate.reactions()),
        },
        "protocol": {
            "reactor": "IdealGasConstPressureReactor",
            "energy": "on",
            "ignition_definition": "linear-interpolated-T-minus-T0-equals-300-K",
            "maximum_time_s": arguments.maximum_time,
            "tolerances": tolerances.__dict__,
            "candidate_reselection_permitted": False,
        },
        "case_metrics": rows,
        "audit_pass": passed,
        "claims": {
            "internal_full_vs_retained_audit_complete": True,
            "experimental_validation_complete": False,
            "mechanism_selected_for_rfam_runtime": False,
            "phase_aware_chemistry_enabled": False,
            "production_ready": False,
        },
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--expected-candidate-sha256", required=True)
    parser.add_argument("--expected-species", type=int, required=True)
    parser.add_argument("--expected-reactions", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum-time", type=float, default=0.005)
    parser.add_argument("--idt-relative-tolerance", type=float, default=0.20)
    parser.add_argument("--inverse-temperature-time-tolerance", type=float, default=0.25)
    parser.add_argument("--heat-release-nrmse-tolerance", type=float, default=0.25)
    parser.add_argument("--species-absolute-tolerance", type=float, default=0.02)
    arguments = parser.parse_args(argv)
    try:
        output = fresh_output(arguments.output, "result-output")
        result = execute(arguments)
        write_fresh(output, canonical_bytes(result))
        print(
            "R2_C3MECH_RETAINED_AUDIT_BLOCKED"
            f" cases={len(result['case_metrics'])}"
            f" audit_pass={int(result['audit_pass'])}"
            " experimental=0 production=0"
        )
        return 0
    except Exception as error:
        reason = str(error).replace("\r", "\\r").replace("\n", "\\n")
        print("R2_C3MECH_RETAINED_AUDIT_FAIL " + reason, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
