import cantera as ct, yaml, re
SRC = "tests/reference/sandiego.yaml"
PROTECTED = set("C3H8 C3H6 N-C3H7 I-C3H7 C2H5 C2H4 C2H3 C2H2 CH4 CH3 CH2O HCO CO CO2 H2O H2 H O OH HO2 H2O2 O2 N2".split())
doc = yaml.safe_load(open(SRC, encoding="utf-8"))
tok = re.compile(r"[A-Za-z][A-Za-z0-9-]*")
out = []
for r in doc.get("reactions", []):
    toks = [t for t in tok.findall(str(r.get("equation", ""))) if t != "M"]
    if any(t not in PROTECTED for t in toks):
        continue
    out.append(r["equation"])
print("kept %d reactions" % len(out))
for e in out:
    print(" ", e)
