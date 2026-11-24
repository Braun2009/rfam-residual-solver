#!/usr/bin/env python3
"""Cross-check the generated C++ NASA-7 table against Cantera (offline).

Reads the generated header's embedded coefficients and evaluates Cp/h for every
species at probe temperatures, comparing against the Cantera 3.1.0 NASA-7
polynomial values. This locks the coefficient transcription to the authoritative
source without requiring Cantera at RFAM build time.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


DEFAULT_HEADER = Path("src/thermo/SandiegoNasa7Data.H")
PROBES_K = (300.0, 600.0, 1100.0, 1500.0, 2500.0)
REL_TOL = 1.0e-9
EDGE_GUARD = 1.0  # K, keep probes off range/mid boundaries (source NASA-7
                  # data has small discontinuities at Tmid and range edges)
RUNIV = 8.314462618


def parse_header(path: Path) -> list[dict]:
    text = path.read_text(encoding="ascii")
    records: list[dict] = []
    # each record block: name, Tlow, Tcommon, Thigh, {lo...}, {hi...}, mw, C H O N
    pattern = re.compile(
        r'"(?P<name>[^"]+)"\s*,\s*'
        r"(?P<tlow>[-\d.eE+]+)\s*,\s*(?P<tcom>[-\d.eE+]+)\s*,\s*(?P<thigh>[-\d.eE+]+)\s*,\s*"
        r"\{\s*(?P<lo>[^}]+)\}\s*,\s*\{\s*(?P<hi>[^}]+)\}\s*,\s*"
        r"(?P<mw>[-\d.eE+]+)\s*,\s*(?P<nc>\d+)\s*,\s*(?P<nh>\d+)\s*,\s*(?P<no>\d+)\s*,\s*(?P<nn>\d+)"
    )
    for match in pattern.finditer(text):
        records.append({key: match.group(key) for key in match.groupdict()})
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    args = parser.parse_args()
    try:
        import cantera as ct
    except ImportError as error:
        print("CANTERA_UNAVAILABLE", file=sys.stderr)
        return 2

    records = parse_header(args.header.resolve())
    if len(records) != 59:
        print(f"SANDIEGO_NASA7_AUDIT_FAIL header_species={len(records)}", file=sys.stderr)
        return 2

    gas = ct.Solution(str(Path("tests/reference/sandiego.yaml").resolve()))
    worst_cp = 0.0
    worst_h = 0.0
    for record in records:
        name = record["name"]
        if name not in gas.species_names:
            print(f"SANDIEGO_NASA7_AUDIT_FAIL unknown={name}", file=sys.stderr)
            return 2
        index = gas.species_index(name)
        low = [float(v) for v in re.split(r",\s*", record["lo"])]
        high = [float(v) for v in re.split(r",\s*", record["hi"])]
        tlow = float(record["tlow"])
        tcom = float(record["tcom"])
        thigh = float(record["thigh"])
        if tlow != gas.species(index).thermo.min_temp or thigh != gas.species(index).thermo.max_temp:
            print(
                f"SANDIEGO_NASA7_AUDIT_FAIL range {name} "
                f"{tlow}/{thigh} vs {gas.species(index).thermo.min_temp}/"
                f"{gas.species(index).thermo.max_temp}",
                file=sys.stderr,
            )
            return 2
        for T in PROBES_K:
            # Pure polynomial evaluation (no clamping), matching Cantera's
            # NasaPoly2; only keep probes off the Tmid split boundary where
            # the source data itself has small discontinuities.
            Tuse = T
            if abs(Tuse - tcom) < EDGE_GUARD:
                Tuse = tcom + EDGE_GUARD if Tuse >= tcom else tcom - EDGE_GUARD
            coeffs = low if Tuse < tcom else high
            cp = RUNIV * (coeffs[0] + Tuse * (coeffs[1] + Tuse * (coeffs[2] + Tuse * (coeffs[3] + Tuse * coeffs[4]))))
            h = RUNIV * (Tuse * (coeffs[0] + Tuse * (coeffs[1] / 2 + Tuse * (coeffs[2] / 3 + Tuse * (coeffs[3] / 4 + Tuse * (coeffs[4] / 5))))) + coeffs[5])
            gas.TP = Tuse, 101325.0
            ref_cp = gas.standard_cp_R[index] * RUNIV
            ref_h = gas.standard_enthalpies_RT[index] * RUNIV * Tuse
            worst_cp = max(worst_cp, abs(cp - ref_cp) / max(abs(ref_cp), 1.0))
            worst_h = max(worst_h, abs(h - ref_h) / max(abs(ref_h), 1.0))
            if worst_cp > REL_TOL or worst_h > REL_TOL:
                print(
                    f"SANDIEGO_NASA7_AUDIT_FAIL {name} T={T} "
                    f"cp={cp} vs {ref_cp} h={h} vs {ref_h}",
                    file=sys.stderr,
                )
                return 2
    print(f"SANDIEGO_NASA7_AUDIT_PASS species={len(records)} "
          f"max_rel_cp={worst_cp:.3e} max_rel_h={worst_h:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
