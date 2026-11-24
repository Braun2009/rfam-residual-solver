#!/usr/bin/env python3
"""vortex_geom_gen.py - generate vortex_full360.geo (gmsh OpenCASCADE).

Verified gmsh 4.8 idioms:
  x = newv; Box(x) = {...};
  u = newv; BooleanUnion(u) = { Volume{a}; Delete; }{ Volume{b}; Delete; };
  r = newv; BooleanDifference(r) = { Volume{a}; Delete; }{ Volume{b}; Delete; };
  out[] = Extrude{ {0,0,1}, {0,0,0}, angle }{ Surface{s}; };

CFD frame: x axial (0=plate face, downstream +), plate centre at origin,
units mm. Fluid = (chamber + plume) - plate - 3 pintle heads.
"""
import json, math

PROFILE = "chamber_profile_3pintle.json"
OUT = "vortex_full360.geo"

def main():
    doc = json.load(open(PROFILE, encoding="utf-8"))
    prof = doc["profile"]
    step = max(1, len(prof) // 120)
    pts = [(p["x_mm"], p["r_mm"]) for p in prof[::step]]
    if pts[-1][0] < 406.5:
        pts.append((406.5, prof[-1]["r_mm"]))
    axes = [(12.5*math.cos(math.radians(a)), 12.5*math.sin(math.radians(a)))
            for a in (-150.7, -30.1, 90.4)]

    L = []
    A = L.append
    A('SetFactory("OpenCASCADE");')
    A('')
    # ---------- chamber ----------
    A('// chamber: meridian section revolved 2pi')
    ptags = []
    for i, (x, r) in enumerate(pts):
        A(f'pt{i} = newp; Point(pt{i}) = {{{x}, 0, {r}, 5}};')
        ptags.append(f'pt{i}')
    A('ax0 = newp; Point(ax0) = {0, 0, 0, 5};')
    A('ax1 = newp; Point(ax1) = {406.5, 0, 0, 5};')
    A(f'spl = newl; BSpline(spl) = {{{", ".join(ptags)}}};')
    A('lax = newl; Line(lax) = {ax0, ax1};')
    A(f'lstart = newl; Line(lstart) = {{ax0, {ptags[0]}}};')
    A(f'lend = newl; Line(lend) = {{{ptags[-1]}, ax1}};')
    A('clch = newll; Curve Loop(clch) = {lax, -lend, -spl, -lstart};')
    A('sch = news; Plane Surface(sch) = {clch};')
    A('outcham[] = Extrude{ {0,0,1}, {0,0,0}, 2*Pi }{ Surface{sch}; };')
    A('cham = outcham[0];')
    A('')
    # ---------- plume ----------
    A('// plume: rectangle section revolved 2pi')
    A('p1 = newp; Point(p1) = {406.5, 60, 0, 10};')
    A('p2 = newp; Point(p2) = {1406.5, 60, 0, 30};')
    A('p3 = newp; Point(p3) = {1406.5, 300, 0, 60};')
    A('p4 = newp; Point(p4) = {406.5, 300, 0, 60};')
    A('qax0 = newp; Point(qax0) = {406.5, 0, 0, 10};')
    A('qax1 = newp; Point(qax1) = {1406.5, 0, 0, 30};')
    A('l1 = newl; Line(l1) = {p1, p2};')
    A('l2 = newl; Line(l2) = {p2, p3};')
    A('l3 = newl; Line(l3) = {p3, p4};')
    A('l4 = newl; Line(l4) = {p4, p1};')
    A('lax2 = newl; Line(lax2) = {qax0, qax1};')
    A('lc1 = newl; Line(lc1) = {qax0, p1};')
    A('lc2 = newl; Line(lc2) = {p2, qax1};')
    A('ll2 = newll; Curve Loop(ll2) = {lax2, -lc2, l2, l3, l4, -lc1};')
    A('sf = news; Plane Surface(sf) = {ll2};')
    A('outplume[] = Extrude{ {0,0,1}, {0,0,0}, 2*Pi }{ Surface{sf}; };')
    A('plume = outplume[0];')
    A('')
    # ---------- plate solid ----------
    A('vplate = newv; Cylinder(vplate) = {0, 0, -8, 0, 0, 8, 24};')
    A('')
    # ---------- pintle heads ----------
    H = 21.8
    head_vars = ['vplate']
    for idx, (yc, zc) in enumerate(axes):
        t = idx + 1
        A(f'// pintle {t}: LOX + C3H8 frusta + hemi cap')
        A('q0 = newp; Point(q0) = {0, 0, 0, 0.6};')
        A('q1 = newp; Point(q1) = {0, 7.5, 0, 0.6};')
        A('q2 = newp; Point(q2) = {26, 17.3, 0, 0.6};')
        A('m1 = newl; Line(m1) = {q0, q1};')
        A('m2 = newl; Line(m2) = {q1, q2};')
        A('m3 = newl; Line(m3) = {q2, q0};')
        A('c1 = newll; Curve Loop(c1) = {m1, m2, m3};')
        A('s1 = news; Plane Surface(s1) = {c1};')
        A(f'vfrust[] = Extrude{{ {{0,0,1}}, {{0,0,0}}, {2*H}*Pi/180 }}{{ Surface{{s1}}; }};')
        A('frust = vfrust[0];')
        A('head = frust;')
        ang = math.degrees(math.atan2(zc, yc))
        A(f'Rotate {{ {{1,0,0}}, {{0,0,0}}, ({ang} - {H})*Pi/180 }}{{ Volume{{head}}; }};')
        A(f'Translate {{ Volume{{head}}; }}{{ 0, {yc}, {zc} }};')
        head_vars.append('head')
    A('')
    # ---------- fluid ----------
    A('// fluid = (chamber + plume) - plate - pintle heads')
    A('dom = newv; BooleanUnion(dom) = { Volume{cham}; Delete; }{ Volume{plume}; Delete; };')
    A('fluid = newv; BooleanDifference(fluid) = { Volume{dom}; Delete; }'
      f'{{ Volume{{{", ".join(head_vars)}}}; Delete; }};')
    A('')
    A('Physical Volume("fluid", 1) = { fluid };')
    A('Printf("GEO_GEN_DONE");')
    open(OUT, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"wrote {OUT}")

if __name__ == "__main__":
    main()
