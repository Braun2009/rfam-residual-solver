#!/usr/bin/env python3
"""Remesh a recovered fluid boundary, retaining and measuring topology.

This experiment records all geometry changes. No CFD quality acceptance is
inferred from success of a remesher.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import pymeshlab as ml
import trimesh
import manifold3d as md

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--mesh',required=True);p.add_argument('--geometry',required=True);p.add_argument('--out',required=True)
p.add_argument('--base-mm',type=float,default=1.5);p.add_argument('--injector-mm',type=float,default=.3)
p.add_argument('--throat-mm',type=float,default=.5);p.add_argument('--max-dist-mm',type=float,default=.015)
a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
d=np.load(Path(a.mesh)/'fluid_volume_mm.npz');points=d['points'];faces=d['boundary']
used=np.unique(faces);ids=np.full(len(points),-1);ids[used]=np.arange(len(used));points=points[used];faces=ids[faces]
manifest=json.loads((Path(a.geometry)/'geometry_manifest.json').read_text())
board=manifest['board_cad_z_mm'];centre=np.array(manifest['centre_xy_mm']);start=time.time()
ms=ml.MeshSet();ms.add_mesh(ml.Mesh(points,faces),'fluid')

def check(stage):
    m=ms.current_mesh();v=m.vertex_matrix();f=m.face_matrix()
    t=trimesh.Trimesh(v,f,process=False)
    report={'stage':stage,'vertices':len(v),'triangles':len(f),'closed':bool(t.is_watertight),
            'euler':int(t.euler_number),'volume_mm3':float(t.volume),'elapsed_s':time.time()-start}
    print(report,flush=True)
    np.savez_compressed(out/f'{stage}.npz',vertices=v,faces=f)
    if not t.is_watertight or not t.is_winding_consistent:raise ValueError('Remeshing opened the surface')
    return report

report=[check('before')]
# Explicit absolute tolerance: 0.0001 mm, about three float32 ULPs at
# the injector. Abort if merging changes the closed-surface topology.
ms.apply_filter('meshing_merge_close_vertices',threshold=ml.PureValue(0.0001))
ms.apply_filter('meshing_remove_duplicate_faces')
ms.apply_filter('meshing_remove_null_faces')
report.append(check('merge'))
if report[-1]['euler']!=report[0]['euler']:raise ValueError('Merge changed genus')
for name,length,selection in [('base',a.base_mm,None),('injector',a.injector_mm,f'(z0 < {board+24} && z1 < {board+24} && z2 < {board+24})'),('throat',a.throat_mm,f'(z0 > {board+254} && z1 > {board+254} && z2 > {board+254} && z0 < {board+310} && z1 < {board+310} && z2 < {board+310})')]:
    if selection:ms.apply_filter('compute_selection_by_condition_per_face',condselect=selection)
    ms.apply_filter('meshing_isotropic_explicit_remeshing',iterations=6,adaptive=False,
                    selectedonly=selection is not None,targetlen=ml.PureValue(length),featuredeg=45,
                    checksurfdist=True,maxsurfdist=ml.PureValue(a.max_dist_mm),
                    splitflag=True,collapseflag=True,swapflag=True,smoothflag=True,reprojectflag=True)
    report.append(check(name))
    if report[-1]['euler']!=report[0]['euler']:raise ValueError(f'{name} changed genus')
m=ms.current_mesh();v=m.vertex_matrix();f=m.face_matrix()
s=md.Manifold(md.Mesh64(v,f.astype(np.uint64)))
if s.status()!=md.Error.NoError:raise ValueError(s.status())
patches=s.slice(manifest['inlet_cad_z_mm']+.02).decompose()
if len(patches)!=64:raise ValueError(f'Inlet connectivity changed: {len(patches)}')
np.savez_compressed(out/'remeshed_boundary_mm.npz',vertices=v,faces=f)
(out/'remesh_report.json').write_text(json.dumps({'stages':report,'inlet_apertures':len(patches),'max_projection_distance_mm':a.max_dist_mm,'targets_mm':{'base':a.base_mm,'injector':a.injector_mm,'throat':a.throat_mm}},indent=2))
print('SURFACE_REMESH_FINISHED',flush=True)
