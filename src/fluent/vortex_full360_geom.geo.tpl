/* vortex_full360_geom.geo - full-360 chamber + plume + injector-head booleans
   CFD frame: x = axial (0 = injector plate face, downstream +x)
              y,z = transverse, plate centre at origin
   Units: mm (gmsh) -> m on export. */
SetFactory("OpenCASCADE");

// ---- chamber profile (v5) ----
// points from chamber_profile_3pintle.json, decimated
// (generated at runtime by vortex_geom_gen.py; this file is a template)
DefineConstant[ NPROF = {0, Name "profile point count"} ];

// ---- three pintle axes (CFD transverse coords, from STL analysis) ----
// theta = -150.7, -30.1, 90.4 deg on bolt circle R=12.5
yc1 = 12.5*Cos(-150.7*Pi/180); zc1 = 12.5*Sin(-150.7*Pi/180);
yc2 = 12.5*Cos(-30.1*Pi/180);  zc2 = 12.5*Sin(-30.1*Pi/180);
yc3 = 12.5*Cos(90.4*Pi/180);   zc3 = 12.5*Sin(90.4*Pi/180);
Printf("pintle axes: (%g,%g) (%g,%g) (%g,%g)", yc1, zc1, yc2, zc2, yc3, zc3);

// ---- profile spline points (placeholders replaced by generator) ----
// pprof[i] = newp; Point(pprof[i]) = {x, 0, r, lc};
