"""Find the 4 margin species that fix IDT@2000K."""
import sys
sys.path.insert(0, "tools")
import cantera as ct
from reduce_sandiego_mechanism import FINAL_30, reduced_gas, idt, SRC

base26 = set("C2H2 C2H3 C2H4 C2H5 C2H6 C3H6 C3H8 CH2CHO CH2O CH3 CH3O CH4 CO CO2 H H2 H2O H2O2 HCO HO2 I-C3H7 N-C3H7 N2 O O2 OH".split())
candidates = [s for s in ct.Solution(SRC).species_names if s not in base26]
full = ct.Solution(SRC)
tf = {t0: idt(ct.Solution(SRC), 4e6, 2.667, t0) for t0 in (1300, 1600, 2000)}
print("full:", {k: f"{v:.2e}" for k, v in tf.items()})

best = None
for c in candidates:
    keep = base26 | {c}
    g = reduced_gas(full, keep)
    score = 0.0
    details = []
    for t0 in (1300, 1600, 2000):
        t = idt(g, 4e6, 2.667, t0)
        if t is None:
            score += 10
            details.append("none")
        else:
            rel = abs(t - tf[t0]) / tf[t0]
            score += rel
            details.append(f"{rel:.1%}")
    print(f"+{c}: {' '.join(details)}")
    if best is None or score < best[1]:
        best = (c, score)
print("best single:", best)
