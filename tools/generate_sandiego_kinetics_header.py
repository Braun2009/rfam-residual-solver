#!/usr/bin/env python3
"""Generate the RFAM C++ static kinetics table for the San Diego mechanism.

Reads tests/reference/sandiego_runtime_bundle_v1.json (Cantera-generated) and
emits src/thermo/SandiegoKineticsData.H with all 271 reactions:
  - per-reaction type (Arrhenius / three-body / falloff-Troe), reversibility
  - high/low Arrhenius parameters converted to the mol,m,s basis (Cantera
    emits kmol,m,s; A *= 1000^(order-1), Ea /= 1000)
  - explicit orders, third-body efficiencies, Troe falloff coefficients
  - Cantera high/low rate constants at the three bundle reference states
    (kmol basis) embedded for C++ cross-checking without Cantera.

The emitted header is pure C++ (no third-party dependency). Species indices
refer to the order of tests/reference/sandiego_runtime_bundle_v1.json, which
matches SandiegoNasa7Data.H.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


DEFAULT_BUNDLE = Path("tests/reference/sandiego_runtime_bundle_v1.json")
DEFAULT_OUTPUT = Path("src/thermo/SandiegoKineticsData.H")


def esc(value: float) -> str:
    return format(value, ".17g")


def order_of(input_data: dict[str, Any]) -> float:
    orders = input_data.get("orders")
    if orders:
        return float(sum(orders.values()))
    reactants = input_data.get("reactants", {})
    return float(sum(reactants.values()))


def order_of(reaction: dict[str, Any]) -> float:
    orders = reaction.get("orders") or {}
    nu = sum(orders.values()) if orders else sum(reaction["reactants"].values())
    # Cantera counts the third body M in the three-body rate law as an extra
    # order for unit conversion: k has units (kmol/m3)^(1-nu-1)/s.
    if reaction["type"] == "three-body-Arrhenius":
        nu += 1.0
    return float(nu)


def arrhenius_params(spec: dict[str, float], order: float) -> str:
    a_kmol = float(spec["A"])
    b = float(spec["b"])
    ea_kmol = float(spec["Ea"])
    # kmol,m,s -> mol,m,s: the rate constant carries (kmol/m3)^(1-n), i.e.
    # (1000 mol/m3)^(1-n) = 1000^(1-n) * (mol/m3)^(1-n), so A_mol =
    # A_kmol * 1000^(1-order).  (The earlier 1000^(order-1) sign error made
    # rates 1000^(2*order-2) too large, which the symmetric test conversion
    # happened to mask; heat release exposed it.)
    a_mol = a_kmol * math.pow(1000.0, 1.0 - order)
    ea_mol = ea_kmol / 1000.0
    return f"{esc(a_mol)}, {esc(b)}, {esc(ea_mol)}, {esc(order)}"


def coeff_list(entries: dict[str, float], species_index: dict[str, int],
               capacity: int) -> list[str]:
    ordered = sorted(entries.items())
    if len(ordered) > capacity:
        raise RuntimeError("coefficient list exceeds static capacity")
    return [f"{{{species_index[name]}, {esc(value)}}}" for name, value in ordered]


def emit_reaction(index: int, reaction: dict[str, Any],
                  species_index: dict[str, int],
                  references: list[dict[str, Any]]) -> list[str]:
    rtype = reaction["type"]
    input_data = reaction["input"]
    order = order_of(reaction)
    if rtype == "Arrhenius":
        type_id = "0"
        high = arrhenius_params(input_data["rate-constant"], order)
        low = "0, 0, 0, 0"
        troe = "0, 0, 0, 0"
        has_low = "false"
    elif rtype == "three-body-Arrhenius":
        type_id = "1"
        high = arrhenius_params(input_data["rate-constant"], order)
        low = "0, 0, 0, 0"
        troe = "0, 0, 0, 0"
        has_low = "false"
    elif rtype == "falloff-Troe":
        type_id = "2"
        high = arrhenius_params(
            input_data["high-P-rate-constant"], order)
        low = arrhenius_params(
            input_data["low-P-rate-constant"], order)
        troe_spec = input_data["Troe"]
        t2 = troe_spec.get("T2", 0.0)
        troe = (f"{esc(troe_spec['A'])}, {esc(troe_spec['T3'])}, "
                f"{esc(troe_spec['T1'])}, {esc(t2)}")
        has_low = "true"
    else:
        raise RuntimeError(f"unsupported reaction type {rtype}")

    reactants = coeff_list(reaction["reactants"], species_index, 6)
    products = coeff_list(reaction["products"], species_index, 6)
    orders = coeff_list(reaction["orders"] or {}, species_index, 6)
    efficiencies: list[str] = []
    third_body = input_data.get("efficiencies")
    if third_body:
        efficiencies = coeff_list(third_body, species_index, 8)

    high_ref = [esc(float(references[i]["high_rate_constants_kmol_m_s"][index]))
                for i in range(len(references))]
    low_ref: list[str] = []
    for i in range(len(references)):
        value = references[i]["low_rate_constants_kmol_m_s"][index]
        low_ref.append(
            esc(float(value)) if value is not None
            else "std::numeric_limits<double>::quiet_NaN()")
    reversible = "true" if reaction["reversible"] else "false"
    duplicate = "true" if reaction["duplicate"] else "false"

    return [
        "    SandiegoReactionData{",
        f'        "{reaction["equation"]}",',
        f"        {type_id}, {reversible}, {duplicate}, {has_low},",
        f"        {{ {high} }},",
        f"        {{ {low} }},",
        f"        {{ {troe} }},",
        f"        {len(reactants)}, {{ {', '.join(reactants)} }},",
        f"        {len(products)}, {{ {', '.join(products)} }},",
        f"        {len(orders)}, {{ {', '.join(orders)} }},",
        f"        {len(efficiencies)}, {{ {', '.join(efficiencies)} }},",
        "        { " + ", ".join(high_ref) + " },",
        "        { " + ", ".join(low_ref) + " }",
        "    },",
    ]


def generate(bundle_path: Path, output_path: Path) -> None:
    document = json.loads(bundle_path.read_text(encoding="ascii"))
    if document["schema"] != "rfam.sandiego-runtime-bundle/v1":
        raise RuntimeError("unsupported bundle schema")
    species = document["species"]
    reactions = document["reactions"]
    references = document["rate_reference_states"]
    if len(species) != 59 or len(reactions) != 271:
        raise RuntimeError("unexpected species/reaction counts")
    species_index = {row["name"]: i for i, row in enumerate(species)}

    lines: list[str] = []
    lines.append("/* Auto-generated by tools/generate_sandiego_kinetics_header.py. */")
    lines.append("/* San Diego 271-reaction kinetics; DO NOT EDIT BY HAND.       */")
    lines.append("")
    lines.append("#ifndef SandiegoKineticsData_H")
    lines.append("#define SandiegoKineticsData_H")
    lines.append("")
    lines.append("#include <cmath>")
    lines.append("#include <limits>")
    lines.append("#include <string>")
    lines.append("")
    lines.append("namespace Foam {")
    lines.append("namespace SandiegoKinetics {")
    lines.append("")
    lines.append("// A/n/Ea are on the mol,m,s basis (concentrations in mol/m3);")
    lines.append("// order is the reaction order used for the kmol<->mol unit")
    lines.append("// conversion (order = sum of explicit orders, or reactants).")
    lines.append("struct ArrheniusParams")
    lines.append("{")
    lines.append("    double A, n, Ea, order;")
    lines.append("};")
    lines.append("")
    lines.append("struct SpeciesCoeff")
    lines.append("{")
    lines.append("    int species;")
    lines.append("    double coeff;")
    lines.append("};")
    lines.append("")
    lines.append("struct SandiegoReactionData")
    lines.append("{")
    lines.append("    const char* equation;")
    lines.append("    int type;              // 0 Arrhenius, 1 three-body, 2 falloff-Troe")
    lines.append("    bool reversible, duplicate, hasLow;")
    lines.append("    ArrheniusParams high;  // k_inf (falloff), k (Arrhenius/three-body)")
    lines.append("    ArrheniusParams low;   // k0 (falloff only)")
    lines.append("    double troe[4];       // falloff Troe: A, T3, T1, T2 (T2=0 if absent)")
    lines.append("    int nReactants;")
    lines.append("    SpeciesCoeff reactants[6];")
    lines.append("    int nProducts;")
    lines.append("    SpeciesCoeff products[6];")
    lines.append("    int nOrders;")
    lines.append("    SpeciesCoeff orders[6];")
    lines.append("    int nEfficiencies;")
    lines.append("    SpeciesCoeff efficiencies[8];")
    lines.append("    double highKRef[3];   // Cantera high k [kmol,m,s] at reference states")
    lines.append("    double lowKRef[3];    // Cantera low k (falloff), NaN otherwise")
    lines.append("};")
    lines.append("")
    lines.append("constexpr int nReactions = 271;")
    lines.append("")
    lines.append("const SandiegoReactionData reactions[nReactions] = {")
    for index, reaction in enumerate(reactions):
        lines.extend(emit_reaction(index, reaction, species_index, references))
    lines.append("};")
    lines.append("")
    lines.append("} // namespace SandiegoKinetics")
    lines.append("} // namespace Foam")
    lines.append("")
    lines.append("#endif")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="ascii")
    print(f"SANDIEGO_KINETICS_HEADER_PASS reactions={len(reactions)} "
          f"output={output_path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        generate(args.bundle.resolve(), args.output.resolve())
    except (OSError, KeyError, ValueError, json.JSONDecodeError, RuntimeError) as error:
        print(f"SANDIEGO_KINETICS_HEADER_FAIL {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
