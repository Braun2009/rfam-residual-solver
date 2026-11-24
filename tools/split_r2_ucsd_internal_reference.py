#!/usr/bin/env python3
"""Split the frozen internal UCSD screen into isolated protocol artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

from screen_r2_detailed_candidate import validate_reference


EXPECTED_SOURCE_SHA256 = (
    "ffcf5fb0c88b1cf2a4be465bfe2b53b9b306012154fc74c02b45a2b5a1c3984f"
)
EXPECTED_CASES = {
    "protocol-training": [
        "train_phi0p8_diluted",
        "train_phi1p0",
        "train_phi1p2",
    ],
    "protocol-holdout": [
        "holdout_phi0p6_diluted",
        "holdout_phi1p1",
        "holdout_phi1p4",
    ],
}


class SplitError(RuntimeError):
    pass


def strict_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise SplitError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(token: str) -> None:
    raise SplitError("nonfinite-json-token:" + token)


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (
        json.dumps(
            document,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def write_atomic(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def split(source: Path) -> dict[str, bytes]:
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SOURCE_SHA256:
        raise SplitError("source-sha256-mismatch")
    try:
        document = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except (UnicodeError, json.JSONDecodeError) as error:
        raise SplitError("source-json-invalid:" + str(error)) from error
    if not isinstance(document, dict):
        raise SplitError("source-root-not-object")
    validate_reference(document)
    actual = {
        partition: [
            case["id"]
            for case in document["cases"]
            if case.get("partition") == partition
        ]
        for partition in EXPECTED_CASES
    }
    if actual != EXPECTED_CASES or len(document["cases"]) != 6:
        raise SplitError("source-partition-mismatch")

    outputs: dict[str, bytes] = {}
    for partition, case_ids in EXPECTED_CASES.items():
        artifact = dict(document)
        artifact["source_reference_sha256"] = EXPECTED_SOURCE_SHA256
        artifact["partition_artifact"] = partition
        by_id = {case["id"]: case for case in document["cases"]}
        artifact["cases"] = [by_id[case_id] for case_id in case_ids]
        outputs[partition] = canonical_bytes(artifact)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--training-output", type=Path, required=True)
    parser.add_argument("--holdout-output", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        outputs = split(arguments.source.resolve())
        write_atomic(
            arguments.training_output.resolve(), outputs["protocol-training"]
        )
        write_atomic(
            arguments.holdout_output.resolve(), outputs["protocol-holdout"]
        )
    except Exception as error:
        print("R2_UCSD_REFERENCE_SPLIT_FAIL " + str(error), file=sys.stderr)
        return 2
    print(
        "R2_UCSD_REFERENCE_SPLIT_PASS training_cases=3 holdout_cases=3 "
        "rights_clear=0 fit_eligible=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
