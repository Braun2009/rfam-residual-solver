#!/usr/bin/env python3
"""Generate the GERG-2008 residual-Helmholtz reference table.

Uses torch-flash (independent Helmholtz multifluid backend, parameter set
multiparameter.gerg-2008) to tabulate the mass-basis residual Helmholtz
energy a_res(T,rho) and its first/second derivatives for pure O2 and C3H8
on a (T, rho) grid.  The GERG residual is defined for EVERY rho (including
the dense liquid branch where the PR isochoric path fails), which is the
motivation for the GERG adoption.

Derivatives are obtained by torch autograd on alpha_residual:
  a_res = R*T*alpha_residual(T, rho_mol, x) / MW
  da_res/dT, d2a_res/dT2, da_res/drho, d2a_res/drho2, d2a_res/dTdrho

A 1e-20 mole-fraction epsilon substitutes for exact zeros in the binary
composition (torch-flash computes 0/0 -> NaN for zero components).

Schema: rfam.gerg-residual-table/v1
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

import torch
import torch_flash as tf


SCHEMA_VERSION = "rfam.gerg-residual-table/v1"
R = 8.314462618
EPS = 1e-20

FLUIDS = {
    "O2": {
        "MW": 0.0319988,
        "temperatures": list(range(60, 501, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0,
                      900.0, 1000.0, 1100.0, 1200.0, 1300.0],
    },
    "C3H8": {
        "MW": 0.044097,
        "temperatures": list(range(90, 651, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0],
    },
}

# Binary mixture tables at fixed O2 mole fractions (GERG-2008 reducing
# functions with the torch-flash parameter set; departure terms are zero
# in that parameter set, so the mixture residual is the standard reducing-
# function mix without interaction correction).  xO2 in mole fraction.
MIXTURES = {
    "O2C3H8_x0.05": {"xO2": 0.05, "MW": 0.0319988*0.05 + 0.044097*0.95,
        "temperatures": list(range(90, 501, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0]},
    "O2C3H8_x0.25": {"xO2": 0.25, "MW": 0.0319988*0.25 + 0.044097*0.75,
        "temperatures": list(range(90, 501, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0]},
    "O2C3H8_x0.5": {"xO2": 0.5, "MW": 0.0319988*0.5 + 0.044097*0.5,
        "temperatures": list(range(90, 501, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0]},
    "O2C3H8_x0.75": {"xO2": 0.75, "MW": 0.0319988*0.75 + 0.044097*0.25,
        "temperatures": list(range(90, 501, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0,
                      800.0, 900.0, 1000.0]},
    "O2C3H8_x0.95": {"xO2": 0.95, "MW": 0.0319988*0.95 + 0.044097*0.05,
        "temperatures": list(range(90, 501, 10)),
        "densities": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0,
                      200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0,
                      900.0, 1000.0, 1100.0, 1200.0, 1300.0]},
}


def make_eos():
    return tf.gerg2008(("oxygen", "propane"))


def residual_at(eos, xO2: float, T: float, rho: float, MW: float):
    """a_res and derivatives (mass basis) at (T, rho) for a fixed xO2."""
    x = torch.tensor([xO2, max(EPS, 1.0 - xO2)])
    Tt = torch.tensor(float(T), requires_grad=True)
    rho_mol = torch.tensor(float(rho) / MW, requires_grad=True)
    ar = eos.alpha_residual(Tt, rho_mol, x)
    a = R * Tt * ar / MW
    gT, grho = torch.autograd.grad(a, (Tt, rho_mol), create_graph=True)
    gTT, = torch.autograd.grad(gT, Tt, create_graph=True)
    gTrho, = torch.autograd.grad(gT, rho_mol, create_graph=True)
    grrho, = torch.autograd.grad(grho, rho_mol, create_graph=True)
    return {
        "a_res": float(a.detach()),
        "da_dT": float(gT.detach()),
        "d2a_dT2": float(gTT.detach()),
        "da_drho_mol": float(grho.detach()),
        "d2a_drho2_mol": float(grrho.detach()),
        "d2a_dTdrho_mol": float(gTrho.detach()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
        default=Path(__file__).resolve().parents[1] /
                 "tests/reference/gerg_residual_table_v1.csv")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        with open(args.output, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        if not rows:
            print("EMPTY")
            return 2
        bad = [r for r in rows if not math.isfinite(float(r["a_res"]))]
        print(f"GERG_RESIDUAL_TABLE_CHECK rows={len(rows)} finite={len(rows)-len(bad)}")
        return 0 if not bad else 2

    eos = make_eos()
    rows = []
    for species, spec in FLUIDS.items():
        idx = 0 if species == "O2" else 1
        for T in spec["temperatures"]:
            for rho in spec["densities"]:
                v = residual_at(eos, 1.0 if idx == 0 else EPS, T, rho,
                                spec["MW"])
                rows.append({
                    "schema_version": SCHEMA_VERSION,
                    "species": species,
                    "temperature_k": T,
                    "rho_kg_m3": rho,
                    "a_res_J_kg": v["a_res"],
                    "da_dT_J_kgK": v["da_dT"],
                    "d2a_dT2_J_kgK2": v["d2a_dT2"],
                    "da_drho_J_kg_per_kgm3": v["da_drho_mol"] / spec["MW"],
                    "d2a_drho2_J_kg_per_kgm32":
                        v["d2a_drho2_mol"] / spec["MW"]**2,
                    "d2a_dTdrho_J_kgK_per_kgm3":
                        v["d2a_dTdrho_mol"] / spec["MW"],
                })
                print(f"{species} {T}K {rho}kg/m3: a_res={v['a_res']:.1f}",
                      flush=True)
    for mix, spec in MIXTURES.items():
        for T in spec["temperatures"]:
            for rho in spec["densities"]:
                v = residual_at(eos, spec["xO2"], T, rho, spec["MW"])
                rows.append({
                    "schema_version": SCHEMA_VERSION,
                    "species": mix,
                    "temperature_k": T,
                    "rho_kg_m3": rho,
                    "a_res_J_kg": v["a_res"],
                    "da_dT_J_kgK": v["da_dT"],
                    "d2a_dT2_J_kgK2": v["d2a_dT2"],
                    "da_drho_J_kg_per_kgm3": v["da_drho_mol"] / spec["MW"],
                    "d2a_drho2_J_kg_per_kgm32":
                        v["d2a_drho2_mol"] / spec["MW"]**2,
                    "d2a_dTdrho_J_kgK_per_kgm3":
                        v["d2a_dTdrho_mol"] / spec["MW"],
                })
                print(f"{mix} {T}K {rho}kg/m3: a_res={v['a_res']:.1f}",
                      flush=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    data = Path(args.output).read_bytes()
    print(f"wrote {len(rows)} rows, sha256={hashlib.sha256(data).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
