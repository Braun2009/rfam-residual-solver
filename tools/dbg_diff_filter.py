"""Diff the two filtering implementations on the 26-species elimination set."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "tools")
import cantera as ct
import yaml as _yaml
import re

SRC = "tests/reference/sandiego.yaml"
ELIM_26 = set("C2H2 C2H3 C2H4 C2H5 C2H6 C3H6 C3H8 CH2CHO CH2O CH3 CH3O CH4 CO CO2 H H2 H2O H2O2 HCO HO2 I-C3H7 N-C3H7 N2 O O2 OH".split())

def idt(g, t0=1300.0):
    g.TPX = t0, 4e6, {"C3H8": 1.0, "O2": 5.0 * 2.667 / 3.6276}
    r = ct.IdealGasConstPressureReactor(g)
    net = ct.ReactorNet([r])
    while net.time < 0.05:
        net.step()
        if g.T > t0 + 400.0:
            return net.time
    return None

def build_a(keep):  # dbg_eliminate style
    doc = _yaml.safe_load(open(SRC, encoding="utf-8"))
    ks = set(keep)
    doc["phases"][0]["species"] = [s for s in doc["phases"][0]["species"] if s in ks]
    tok = re.compile(r"[A-Za-z][A-Za-z0-9-]*")
    out = []
    for r in doc.get("reactions", []):
        toks = [t for t in tok.findall(str(r.get("equation", ""))) if t != "M"]
        if any(t not in ks for t in toks):
            continue
        if "efficiencies" in r:
            r["efficiencies"] = {k: v for k, v in r["efficiencies"].items() if k in ks}
        out.append(r)
    doc["reactions"] = out
    return ct.Solution(yaml=_yaml.dump(doc, sort_keys=False))

from reduce_sandiego_mechanism import _yaml_dump_reduced, reduced_gas

g_a = build_a(ELIM_26)
g_b = reduced_gas(ct.Solution(SRC), ELIM_26)
print("A rx:", g_a.n_reactions, " B rx:", g_b.n_reactions)
print("A idt1300:", idt(g_a))
print("B idt1300:", idt(g_b))
# diff reaction lists
ya = _yaml.safe_load(open(SRC, encoding="utf-8"))
tok = re.compile(r"[A-Za-z][A-Za-z0-9-]*")
def rxset(gas, keep):
    doc = _yaml.safe_load(open(SRC, encoding="utf-8"))
    ks = set(keep)
    out = set()
    for r in doc.get("reactions", []):
        toks = [t for t in tok.findall(str(r.get("equation", ""))) if t != "M"]
        if any(t not in ks for t in toks):
            continue
        out.add(r["equation"])
    return out
ra = rxset(None, ELIM_26)
# B's equation set from its yaml
yb = _yaml.safe_load(_yaml_dump_reduced(ct.Solution(SRC), ELIM_26))
rb = set(r["equation"] for r in yb["reactions"])
print("only in A:", ra - rb)
print("only in B:", rb - ra)
