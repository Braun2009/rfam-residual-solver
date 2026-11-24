"""Greedy elimination: 59 -> smallest passing mechanism.

Anchors: adiabatic T (<1%), IDT @1300/1600/2000K (<15%).
PROTECTED species never removed.
Parallel over candidate deletions.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cantera as ct
import yaml as _yaml
import re
import numpy as np
from concurrent.futures import ProcessPoolExecutor

SRC = "tests/reference/sandiego.yaml"
PROTECTED = set("C3H8 C3H6 N-C3H7 I-C3H7 C2H5 C2H4 C2H3 C2H2 CH4 CH3 CH2O HCO CO CO2 H2O H2 H O OH HO2 H2O2 O2 N2".split())

_tok = re.compile(r"[A-Za-z][A-Za-z0-9-]*")

def build(keep):
    doc = _yaml.safe_load(open(SRC, encoding="utf-8"))
    ks = set(keep)
    doc["phases"][0]["species"] = [s for s in doc["phases"][0]["species"] if s in ks]
    out = []
    for r in doc.get("reactions", []):
        toks = [t for t in _tok.findall(str(r.get("equation", ""))) if t != "M"]
        if any(t not in ks for t in toks):
            continue
        if "efficiencies" in r:
            r["efficiencies"] = {k: v for k, v in r["efficiencies"].items() if k in ks}
        out.append(r)
    doc["reactions"] = out
    return ct.Solution(yaml=_yaml.dump(doc, sort_keys=False))

def idt(g, t0):
    g.TPX = t0, 4e6, {"C3H8": 1.0, "O2": 5.0 * 2.667 / 3.6276}
    r = ct.IdealGasConstPressureReactor(g)
    net = ct.ReactorNet([r])
    while net.time < 0.05:
        net.step()
        if g.T > t0 + 400.0:
            return net.time
    return None

def adia(g):
    g.TPX = 500.0, 4e6, {"C3H8": 1.0, "O2": 5.0 * 2.667 / 3.6276}
    g.equilibrate("HP")
    return g.T

def _fill_full():
    global FULL_IDT, FULL_ADIA
    if FULL_IDT:
        return
    for t0 in (1300, 1600, 2000):
        FULL_IDT[t0] = idt(ct.Solution(SRC), t0)
    FULL_ADIA = adia(ct.Solution(SRC))

FULL_IDT = {}
FULL_ADIA = None
_fill_full()

def check(keep):
    """Return True if all anchors pass."""
    try:
        g = build(keep)
    except Exception:
        return False
    T = adia(g)
    if abs(T - FULL_ADIA) / FULL_ADIA >= 0.01:
        return False
    for t0 in (1300, 1600, 2000):
        t = idt(g, t0)
        if t is None:
            return False
        if abs(t - FULL_IDT[t0]) / FULL_IDT[t0] >= 0.15:
            return False
    return True

def _try_delete(c):
    full_names = set(ct.Solution(SRC).species_names)
    keep = (full_names - {c}) | PROTECTED
    return c, check(keep)

def main():
    _fill_full()
    full_names = set(ct.Solution(SRC).species_names)
    keep = set(full_names)
    print(f"start: {len(keep)} species")
    removable = sorted(full_names - PROTECTED)
    removed = []
    while removable:
        # parallel test all candidates
        results = {}
        with ProcessPoolExecutor(max_workers=14) as ex:
            for c, ok in ex.map(_try_delete, removable):
                results[c] = ok
        deletable = [c for c in removable if results.get(c)]
        if not deletable:
            break
        # remove the one with... all equally good; pick first
        c = deletable[0]
        keep.remove(c)
        removable.remove(c)
        removed.append(c)
        print(f"removed {c}: {len(keep)} left")
    print(f"\nfinal: {len(keep)} species")
    print("kept:", ", ".join(sorted(keep)))
    print("removed:", ", ".join(removed))

if __name__ == "__main__":
    main()
