import cantera as ct
import re

g = ct.Solution("tests/reference/sandiego30_fluent.yaml")
hdr = open("src/fluent/GergUdrgmNasa7.H", encoding="utf-8").read()
RU = 8.314462618

def nasa_cp(a, T):
    return RU * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)

def nasa_h(a, T):
    return RU * T * (a[0] + a[1]*T/2 + a[2]*T**2/3 + a[3]*T**3/4
                     + a[4]*T**4/5 + a[5]/T)

pat = re.compile(r'\{"([^"]+)", ([0-9.]+), ([0-9.]+), ([0-9.]+), ([0-9.]+),\s*\{([^}]*)\}, \{([^}]*)\}\}')
ok = True
for m in pat.finditer(hdr):
    name = m.group(1)
    if name not in g.species_names:
        continue
    lo = [float(x) for x in m.group(6).split(",")]
    hi = [float(x) for x in m.group(7).split(",")]
    i = g.species_index(name)
    mw = g.molecular_weights[i] / 1000.0
    for T, a in ((500.0, lo), (2000.0, hi), (3000.0, hi)):
        g.TPX = T, 4e6, {name: 1.0}
        cp_ref = g.cp_mass * mw
        h_ref = g.enthalpy_mass * mw
        dcp = abs(nasa_cp(a, T) - cp_ref) / cp_ref
        dh = abs(nasa_h(a, T) - h_ref) / max(abs(h_ref), 1.0)
        if dcp > 1e-4 or dh > 1e-4:
            print("MISMATCH %s @%gK cp=%g ref=%g h=%g ref=%g" % (
                name, T, nasa_cp(a, T), cp_ref, nasa_h(a, T), h_ref))
            ok = False
print("NASA7_HEADER_" + ("PASS" if ok else "FAIL"))
