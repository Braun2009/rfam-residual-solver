import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import cantera as ct
import yaml as _yaml
import re

SRC = "tests/reference/sandiego.yaml"
PROTECTED = {"C3H8", "C3H6", "N-C3H7", "I-C3H7", "C2H5", "C2H4", "C2H3", "C2H2",
             "CH4", "CH3", "CH2O", "HCO", "CO", "CO2", "H2O", "H2",
             "H", "O", "OH", "HO2", "H2O2", "O2", "N2"}

def filt(keep):
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
    return ct.Solution(yaml=_yaml.dump(doc, sort_keys=False)), len(out)

def idt(g, t0=1300.0):
    g.TPX = t0, 4e6, {"C3H8": 1.0, "O2": 5.0 * 2.667 / 3.6276}
    r = ct.IdealGasConstPressureReactor(g)
    net = ct.ReactorNet([r])
    while net.time < 0.05:
        net.step()
        if g.T > t0 + 400.0:
            return net.time
    return None

g23, n23 = filt(PROTECTED)
print("PROTECTED-23: rx=%d idt1300=%s" % (n23, idt(g23)))

full = ct.Solution(SRC)
candidates = [s for s in full.species_names if s not in PROTECTED]
print("candidates: %d" % len(candidates))
for c in candidates:
    g, n = filt(PROTECTED | {c})
    t = idt(g)
    if t is not None:
        print("+%s: rx=%d idt1300=%.2e IGNITES" % (c, n, t))
