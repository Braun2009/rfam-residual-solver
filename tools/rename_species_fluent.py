"""rename_species_fluent.py - rename S-CH2/I-C3H7/N-C3H7 to Fluent-safe
names at the yaml level, then re-export Chemkin via Cantera."""
import yaml
import cantera as ct

RENAME = {"S-CH2": "CH2S", "I-C3H7": "IC3H7", "N-C3H7": "NC3H7"}

doc = yaml.safe_load(open("tests/reference/sandiego30.yaml", encoding="utf-8"))
# rename in phase species list
doc["phases"][0]["species"] = [RENAME.get(s, s) for s in doc["phases"][0]["species"]]
# rename species definitions
for sp in doc["species"]:
    if sp["name"] in RENAME:
        sp["name"] = RENAME[sp["name"]]
# rename inside reactions
for r in doc.get("reactions", []):
    for f in ("equation",):
        if f in r:
            eq = r[f]
            for old, new in RENAME.items():
                eq = eq.replace(old, new)
            r[f] = eq
    if "efficiencies" in r:
        r["efficiencies"] = {RENAME.get(k, k): v
                             for k, v in r["efficiencies"].items()}
    if "orders" in r:
        r["orders"] = {RENAME.get(k, k): v for k, v in r["orders"].items()}

yaml_text = yaml.dump(doc, sort_keys=False)
g = ct.Solution(yaml=yaml_text)
open("tests/reference/sandiego30_fluent.yaml", "w", encoding="utf-8").write(yaml_text)
import os
out = "tests/reference/sandiego30_chem.inp"
if os.path.exists(out):
    os.remove(out)
g.write_chemkin(out, quiet=True)
print("renamed + re-exported:", g.n_species, "sp", g.n_reactions, "rx")
print([s for s in g.species_names if "CH2" in s or "C3H7" in s])
