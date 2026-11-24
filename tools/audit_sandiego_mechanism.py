#!/usr/bin/env python3
"""Audit the checked-in San Diego Cantera mechanism against RFAM runtime gates.

This is an offline evidence tool. It does not convert the mechanism, mutate
RFAM source files, or claim that the 59-species mechanism is runtime-ready.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import warnings
from collections import Counter
from pathlib import Path
from typing import Any


DEFAULT_MECHANISM = Path("tests/reference/sandiego.yaml")
DEFAULT_OUTPUT = Path("tmp/sandiego_mechanism_audit.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rfam_registry_species(source: Path) -> list[str]:
    text = source.read_text(encoding="utf-8")
    names = re.findall(r'entries_\.push_back\(\{"([^\"]+)"', text)
    if not names:
        raise RuntimeError(f"cannot find RFAM registry entries in {source}")
    return names


def rfam_nasa_species(source: Path) -> list[str]:
    text = source.read_text(encoding="utf-8")
    names = re.findall(r'^\s*add\("([^\"]+)"\s*,', text, flags=re.MULTILINE)
    if not names:
        raise RuntimeError(f"cannot find RFAM NASA7 entries in {source}")
    return names


def species_temperature_ranges(mechanism: Any) -> dict[str, list[float]]:
    ranges: dict[str, list[float]] = {}
    for species in mechanism.species():
        data = species.thermo.input_data
        values = data.get("temperature-ranges", [])
        if values:
            ranges[species.name] = [float(values[0]), float(values[-1])]
    return ranges


def audit(mechanism_path: Path, root: Path) -> dict[str, Any]:
    try:
        import cantera as ct
    except ImportError as error:
        raise RuntimeError("CANTERA_UNAVAILABLE") from error

    # Cantera reports known NASA7 midpoint discontinuities while loading this
    # source. They are evidence to preserve in a detailed thermo audit, but do
    # not affect the structural coverage result emitted by this tool.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        mechanism = ct.Solution(str(mechanism_path))
    species_names = list(mechanism.species_names)
    reactions = list(mechanism.reactions())
    reaction_types = Counter(reaction.reaction_type for reaction in reactions)
    falloff_types = Counter(
        reaction.reaction_type for reaction in reactions
        if reaction.reaction_type.startswith("falloff-")
    )
    third_body_count = sum(
        reaction.third_body is not None for reaction in reactions
    )
    duplicate_count = sum(
        bool(reaction.input_data.get("duplicate", False))
        for reaction in reactions
    )
    reversible_count = sum(
        bool(getattr(reaction, "reversible", False)) for reaction in reactions
    )

    registry = rfam_registry_species(root / "src/thermo/SpeciesRegistry.H")
    nasa7 = rfam_nasa_species(root / "src/thermo/NASA7Thermo.H")
    mechanism_set = set(species_names)
    registry_set = set(registry)
    nasa7_set = set(nasa7)
    ranges = species_temperature_ranges(mechanism)

    missing_registry = sorted(mechanism_set - registry_set)
    missing_nasa7 = sorted(mechanism_set - nasa7_set)
    registry_extra = sorted(registry_set - mechanism_set)
    nasa7_extra = sorted(nasa7_set - mechanism_set)
    minimum_temperature = min(
        values[0] for values in ranges.values() if values
    )
    maximum_temperature = max(
        values[1] for values in ranges.values() if values
    )

    result: dict[str, Any] = {
        "schema": "rfam.sandiego-mechanism-audit/v1",
        "status": "blocked-runtime-integration",
        "mechanism": {
            "path": str(mechanism_path),
            "sha256": sha256(mechanism_path),
            "cantera_version": ct.__version__,
            "species_count": len(species_names),
            "reaction_count": len(reactions),
            "species_names": species_names,
            "reaction_type_counts": dict(sorted(reaction_types.items())),
            "falloff_type_counts": dict(sorted(falloff_types.items())),
            "third_body_reaction_count": third_body_count,
            "duplicate_reaction_count": duplicate_count,
            "reversible_reaction_count": reversible_count,
            "nasa7_temperature_range_K": [
                minimum_temperature, maximum_temperature
            ],
        },
        "rfam_runtime_coverage": {
            "registry_species_count": len(registry),
            "registry_species": registry,
            "nasa7_species_count": len(nasa7),
            "nasa7_species": nasa7,
            "mechanism_species_missing_from_registry": missing_registry,
            "mechanism_species_missing_from_nasa7": missing_nasa7,
            "registry_species_not_in_mechanism": registry_extra,
            "nasa7_species_not_in_mechanism": nasa7_extra,
        },
        "runtime_blockers": [
            "SpeciesRegistry and NASA7Thermo do not cover the full mechanism species set.",
            "DynamicChemistry accepts the RFAM custom YAML schema, not Cantera ck2yaml reaction entries.",
            "Three-body and falloff/Troe pressure dependence are not represented by DynamicChemistry.",
            "A 59-species Domain requires a matching thermodynamic model and EOS species basis.",
        ],
        "runtime_ready": False,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", type=Path, default=DEFAULT_MECHANISM)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    mechanism = args.mechanism.resolve()
    root = args.root.resolve()
    if not mechanism.is_file():
        print(f"MECHANISM_NOT_FOUND: {mechanism}", file=sys.stderr)
        return 2
    try:
        result = audit(mechanism, root)
    except RuntimeError as error:
        print(str(error), file=sys.stderr)
        return 2

    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        "SANDIEGO_MECHANISM_AUDIT_BLOCKED "
        f"species={result['mechanism']['species_count']} "
        f"reactions={result['mechanism']['reaction_count']} "
        f"falloff={sum(result['mechanism']['falloff_type_counts'].values())} "
        f"missing_registry={len(result['rfam_runtime_coverage']['mechanism_species_missing_from_registry'])} "
        f"missing_nasa7={len(result['rfam_runtime_coverage']['mechanism_species_missing_from_nasa7'])}"
    )
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
