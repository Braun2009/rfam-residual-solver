#!/usr/bin/env python3
"""vortex_geom_api.py - full-360 chamber+plume+injector fluid domain + mesh
(gmsh Python API 4.15.2, Anaconda on A406).

CFD frame: x axial (0=plate face, downstream +), plate centre at origin, mm.
Fluid = (chamber + plume) - plate - 3 pintle frusta.
Size fields: wall 0.08mm first layer, cone 0.4mm, throat 0.6mm,
shock region 1mm, far field 14-45mm.
"""
import json, math
import gmsh

PROFILE = "chamber_profile_3pintle.json"

def main():
    doc = json.load(open(PROFILE, encoding="utf-8"))
    prof = doc["profile"]
    step = max(1, len(prof) // 120)
    pts = [(p["x_mm"], p["r_mm"]) for p in prof[::step]]
    if pts[-1][0] < 406.5:
        pts.append((406.5, prof[-1]["r_mm"]))
    axes = [(12.5*math.cos(math.radians(a)), 12.5*math.sin(math.radians(a)))
            for a in (-150.7, -30.1, 90.4)]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("vortex360")
    occ = gmsh.model.occ

    # ---------- chamber ----------
    pids = []
    for (x, r) in pts:
        pids.append(occ.addPoint(x, 0.0, r))
    ax0 = occ.addPoint(0, 0, 0)
    ax1 = occ.addPoint(406.5, 0, 0)
    spl = occ.addBSpline(pids)
    lax = occ.addLine(ax0, ax1)
    lstart = occ.addLine(ax0, pids[0])
    lend = occ.addLine(pids[-1], ax1)
    cl = occ.addCurveLoop([lax, -lend, -spl, -lstart])
    sch = occ.addPlaneSurface([cl])
    cham = occ.revolve([(2, sch)], 0, 0, 0, 1, 0, 0, 2*math.pi)
    cham_vol = [t for d, t in cham if d == 3]

    # ---------- plume ----------
    p1 = occ.addPoint(406.5, 60, 0); p2 = occ.addPoint(1406.5, 60, 0)
    p3 = occ.addPoint(1406.5, 300, 0); p4 = occ.addPoint(406.5, 300, 0)
    q0 = occ.addPoint(406.5, 0, 0); q1 = occ.addPoint(1406.5, 0, 0)
    l1 = occ.addLine(p1, p2); l2 = occ.addLine(p2, p3)
    l3 = occ.addLine(p3, p4); l4 = occ.addLine(p4, p1)
    la = occ.addLine(q0, q1); lc1 = occ.addLine(q0, p1); lc2 = occ.addLine(p2, q1)
    cl2 = occ.addCurveLoop([la, -lc2, l2, l3, l4, -lc1])
    sf = occ.addPlaneSurface([cl2])
    plume = occ.revolve([(2, sf)], 0, 0, 0, 1, 0, 0, 2*math.pi)
    plume_vol = [t for d, t in plume if d == 3]

    # ---------- plate solid ----------
    plate = occ.addCylinder(0, 0, -8, 0, 0, 8, 24)

    # ---------- pintle frusta ----------
    heads = [plate]
    for i, (yc, zc) in enumerate(axes):
        ang = math.degrees(math.atan2(zc, yc))
        t0 = occ.addPoint(0, 0, 0)
        t1 = occ.addPoint(0, 7.5, 0)
        t2 = occ.addPoint(26, 17.3, 0)
        m1 = occ.addLine(t0, t1); m2 = occ.addLine(t1, t2); m3 = occ.addLine(t2, t0)
        c1 = occ.addCurveLoop([m1, m2, m3])
        s1 = occ.addPlaneSurface([c1])
        frust = occ.revolve([(2, s1)], 0, 0, 0, 1, 0, 0, 43.6*math.pi/180)
        frust_vol = [t for d, t in frust if d == 3][0]
        occ.rotate([(3, frust_vol)], 0, 0, 0, 1, 0, 0, (ang - 21.8)*math.pi/180)
        occ.translate([(3, frust_vol)], 0, yc, zc)
        heads.append(frust_vol)
        print(f"pintle {i+1}: angle {ang:.1f} deg, vol {frust_vol}")

    # ---------- boolean ----------
    dom = occ.fuse([(3, cham_vol[0])], [(3, plume_vol[0])])[0]
    fluid = occ.cut(dom, [(3, int(h)) for h in heads], removeTool=False)[0]
    print("fluid:", fluid)
    occ.synchronize()
    # delete everything except the fluid volume for a clean CAD export
    fluid_vol = [t for d, t in fluid if d == 3]
    all3d = [t for d, t in gmsh.model.getEntities(3)]
    to_del = [t for t in all3d if t not in fluid_vol]
    if to_del:
        occ.remove([(3, t) for t in to_del], recursive=True)
    occ.synchronize()

    gmsh.model.addPhysicalGroup(3, [t for d, t in fluid if d == 3], 1)
    gmsh.model.setPhysicalName(3, 1, "fluid")

    bnd = gmsh.model.getBoundary(fluid, oriented=False)
    print("boundary surfaces:", len(bnd))

    # classify + physical groups
    wall_dt, cone_dt, far_dt = [], [], []
    phy = {}
    for (d, tag) in bnd:
        com = occ.getCenterOfMass(d, tag)
        x, y, z = com
        r = math.hypot(y, z)
        name = None
        if r > 299.0 or x > 1400.0:
            name = "farfield"
            far_dt.append(tag)
        elif x < -0.5:
            name = "wall_plate_back"
            wall_dt.append(tag)
        elif r <= 17.5 and x <= 26.5:
            if r < 2.5 and x > 24.0:
                name = "inlet_torch"
            elif x < 13.0:
                name = "inlet_lox"
            else:
                name = "inlet_fuel"
            cone_dt.append(tag)
        elif r > 17.0 and r < 25.0 and x < 0.5:
            name = "wall_plate"
            wall_dt.append(tag)
        elif x < 407.0 and r > 17.0:
            name = "wall_chamber"
            wall_dt.append(tag)
        if name:
            phy.setdefault(name, []).append(tag)
    for name, tags in phy.items():
        g = gmsh.model.addPhysicalGroup(2, tags)
        gmsh.model.setPhysicalName(2, g, name)
    print("physics:", {k: len(v) for k, v in phy.items()})

    # ---------- per-surface mesh sizes (no background field) ----------
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 4.0)
    for tag in wall_dt:
        gmsh.model.mesh.setSize([(2, tag)], 0.5)
    for tag in cone_dt:
        gmsh.model.mesh.setSize([(2, tag)], 1.0)
    for tag in far_dt:
        gmsh.model.mesh.setSize([(2, tag)], 25.0)
    # throat region: refine surfaces near the throat by size on wall bands
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 6)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    if __import__("os").environ.get("VORTEX_GEOM_ONLY"):
        # fine surface mesh for STL export
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 3.0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 15.0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
        occ.synchronize()
        bnd3 = gmsh.model.getBoundary([(3, t) for d, t in fluid if d == 3], oriented=False)
        for (d, tag) in bnd3:
            com = gmsh.model.occ.getCenterOfMass(d, tag)
            x, y, z = com
            r = math.hypot(y, z)
            if r > 299.0 or x > 1400.0:
                gmsh.model.mesh.setSize([(2, tag)], 15.0)
            elif r <= 17.5 and x <= 26.5:
                gmsh.model.mesh.setSize([(2, tag)], 0.5)
            else:
                gmsh.model.mesh.setSize([(2, tag)], 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write("vortex_full360_surface.msh")
        gmsh.write("vortex_full360_fluid_fine.stl")
        gmsh.finalize()
        print("GEOM_EXPORT_OK")
        return
    occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 4.0)
    bnd2 = gmsh.model.getBoundary([(3, t) for d, t in fluid if d == 3], oriented=False)
    for (d, tag) in bnd2:
        com = gmsh.model.occ.getCenterOfMass(d, tag)
        x, y, z = com
        r = math.hypot(y, z)
        if r > 299.0 or x > 1400.0 or x < -0.5:
            gmsh.model.mesh.setSize([(2, tag)], 20.0)
        elif r <= 17.5 and x <= 26.5:
            gmsh.model.mesh.setSize([(2, tag)], 1.0)
        else:
            gmsh.model.mesh.setSize([(2, tag)], 0.5)
    print("meshing ...")
    print("BND2_N:", len(bnd2), "first:", bnd2[:4])
    gmsh.model.mesh.setSize([bnd2[0]], 1.5)
    print("SIZE_PROBE:", gmsh.model.mesh.getSizes(bnd2)[:4])
    gmsh.model.mesh.generate(3)
    ntet = len(gmsh.model.mesh.getElementsByType(4)[1])
    print("MESH_DONE tets=", ntet)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("vortex_full360.msh")
    gmsh.write("vortex_full360.step")
    gmsh.finalize()
    print("GEOM_OK")

if __name__ == "__main__":
    main()
