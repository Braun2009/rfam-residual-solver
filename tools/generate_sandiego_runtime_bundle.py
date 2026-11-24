#!/usr/bin/env python3
"""Generate deterministic RFAM runtime data and Cantera rate references.

Cantera is an offline authority only. The generated JSON contains the full
59-species/271-reaction San Diego mechanism semantics needed by the future C++
loader, plus per-reaction rates at representative chamber states.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import warnings
from pathlib import Path
from typing import Any


SCHEMA = "rfam.sandiego-runtime-bundle/v1"
DEFAULT_INPUT = Path("tests/reference/sandiego.yaml")
DEFAULT_OUTPUT = Path("tests/reference/sandiego_runtime_bundle_v1.json")
REFERENCE_STATES = (
    ("lean-hot", 1800.0, 4.0e6, 0.8),
    ("stoich-hot", 2500.0, 4.0e6, 1.0),
    ("rich-hot", 3800.0, 4.0e6, 1.35),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: Any) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("non-finite mechanism value")
    return result


def normalize(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): normalize(item) for key, item in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [normalize(item) for item in value]
    if hasattr(value, "tolist"):
        return normalize(value.tolist())
    if isinstance(value, float):
        return finite(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return str(value)


def species_row(species: Any) -> dict[str, Any]:
    thermo = normalize(species.thermo.input_data)
    if thermo.get("model") != "NASA7":
        raise ValueError(f"unsupported thermo model for {species.name}")
    return {
        "name": species.name,
        "molecular_weight_kg_per_mol": finite(species.molecular_weight) / 1000.0,
        "elements": {key: finite(value) for key, value in sorted(species.composition.items())},
        "thermo": thermo,
    }


def reaction_row(index: int, reaction: Any) -> dict[str, Any]:
    return {
        "index": index,
        "equation": reaction.equation,
        "type": reaction.reaction_type,
        "reversible": bool(reaction.reversible),
        "duplicate": bool(reaction.duplicate),
        "reactants": {key: finite(value) for key, value in sorted(reaction.reactants.items())},
        "products": {key: finite(value) for key, value in sorted(reaction.products.items())},
        "orders": {key: finite(value) for key, value in sorted((reaction.orders or {}).items())},
        "input": normalize(reaction.input_data),
    }


def reference_state(gas: Any, reactions: list[Any], state_id: str,
                    temperature: float, pressure: float,
                    equivalence_ratio: float) -> dict[str, Any]:
    gas.TP = temperature, pressure
    gas.set_equivalence_ratio(equivalence_ratio, "C3H8:1", "O2:1")
    # Equilibrate at fixed T/P so the reference composition contains radicals
    # (H/OH/O/...) and every reaction has a non-trivial forward/reverse rate
    # to cross-check, while net production rates vanish (equilibrium).
    gas.equilibrate("TP")
    high: list[float] = []
    low: list[float | None] = []
    for reaction in reactions:
        rate = reaction.rate
        if reaction.reaction_type == "falloff-Troe":
            high.append(finite(rate.high_rate(temperature)))
            low.append(finite(rate.low_rate(temperature)))
        else:
            high.append(finite(rate(temperature)))
            low.append(None)
    return {
        "id": state_id,
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "equivalence_ratio": equivalence_ratio,
        "mole_fractions": {
            name: finite(gas.X[i]) for i, name in enumerate(gas.species_names)
            if gas.X[i] != 0.0
        },
        "forward_rates_of_progress_kmol_m3_s": [finite(v) for v in gas.forward_rates_of_progress],
        "reverse_rates_of_progress_kmol_m3_s": [finite(v) for v in gas.reverse_rates_of_progress],
        "net_rates_of_progress_kmol_m3_s": [finite(v) for v in gas.net_rates_of_progress],
        "net_production_rates_kmol_m3_s": [finite(v) for v in gas.net_production_rates],
        "heat_release_rate_W_m3": finite(gas.heat_release_rate),
        # High/low forward rate constants [kmol,m,s] at this state. These are
        # independent of the mixture composition (low = k0 for falloff, null
        # for non-falloff), so the C++ side can cross-check its own Arrhenius
        # evaluation without re-deriving concentrations.
        "high_rate_constants_kmol_m_s": high,
        "low_rate_constants_kmol_m_s": low,
        # Concentration-based equilibrium constants Kc [kmol,m,s] (Cantera
        # equilibrium_constants) at this state, for reverse-rate checks.
        "equilibrium_constants": [finite(v) for v in gas.equilibrium_constants],
    }


def generate(source: Path) -> dict[str, Any]:
    try:
        import cantera as ct
    except ImportError as error:
        raise RuntimeError("CANTERA_UNAVAILABLE") from error
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        gas = ct.Solution(str(source))
    species = [species_row(item) for item in gas.species()]
    cantera_reactions = list(gas.reactions())
    reactions = [reaction_row(i, item) for i, item in enumerate(cantera_reactions)]
    references = [
        reference_state(gas, cantera_reactions, *state)
        for state in REFERENCE_STATES
    ]
    types: dict[str, int] = {}
    for reaction in reactions:
        types[reaction["type"]] = types.get(reaction["type"], 0) + 1
    return {
        "schema": SCHEMA,
        "source": {
            "path": "tests/reference/" + source.name,
            "sha256": sha256(source),
            "cantera_version": ct.__version__,
        },
        "species_count": len(species),
        "reaction_count": len(reactions),
        "reaction_type_counts": dict(sorted(types.items())),
        "species": species,
        "reactions": reactions,
        "rate_reference_states": references,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    document = generate(args.input.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(document, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n",
        encoding="ascii",
    )
    print(
        f"SANDIEGO_RUNTIME_BUNDLE_PASS species={document['species_count']} "
        f"reactions={document['reaction_count']} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
