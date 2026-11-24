#!/usr/bin/env python3
"""Prepare remeshed boundary planes and measure sampled CAD deviation.

Input is a recovered tetrahedral boundary or a vertices/faces NPZ. All units
are millimetres. Positive merge tolerance is opt-in and always audited.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pymeshlab as ml
import trimesh


def load(path):
    d=np.load(path)
    v=d['vertices'] if 'vertices' in d else d['points']
    f=d['faces'] if 'faces' in d else d['boundary']
    used=np.unique(f);ids=np.full(len(v),-1,dtype=np.int64);ids[used]=np.arange(len(used))
    return trimesh.Trimesh(v[used],ids[f],process=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True);p.add_argument('--reference',required=True)
    p.add_argument('--geometry',required=True);p.add_argument('--out',required=True)
    p.add_argument('--merge-mm',type=float,default=0.0)
    a=p.parse_args()
    m=load(a.input);ref=load(a.reference)
    manifest=json.loads((Path(a.geometry)/'geometry_manifest.json').read_text())
    before=m.volume;start_euler=m.euler_number;snaps=[]
    if a.merge_mm:
        if not 0<a.merge_mm<=.001:raise ValueError('Merge tolerance must be <=0.001 mm')
        ms=ml.MeshSet();ms.add_mesh(ml.Mesh(m.vertices,m.faces))
        ms.apply_filter('meshing_merge_close_vertices',threshold=ml.PureValue(a.merge_mm))
        ms.apply_filter('meshing_remove_duplicate_faces');ms.apply_filter('meshing_remove_null_faces')
        mm=ms.current_mesh();m=trimesh.Trimesh(mm.vertex_matrix(),mm.face_matrix(),process=False)
    for z in [manifest['inlet_cad_z_mm'],ref.bounds[1,2]]:
        use=np.abs(m.vertices[:,2]-z)<1e-5
        snaps.append({'z_mm':float(z),'vertices':int(use.sum()),'max_displacement_mm':float(np.abs(m.vertices[use,2]-z).max()) if use.any() else 0})
        m.vertices[use,2]=z
    if not m.is_watertight or not m.is_winding_consistent or m.euler_number!=start_euler:
        raise ValueError('Boundary preparation changed closed topology')
    report={'merge_tolerance_mm':a.merge_mm,'plane_snaps':snaps,'closed':True,'euler':m.euler_number,
            'volume_change_during_preparation_mm3':m.volume-before,
            'relative_volume_difference_from_reference':(m.volume-ref.volume)/ref.volume,
            'sample_count_per_direction':5000,'sampling_is_not_a_Hausdorff_bound':True,'distances':{}}
    rng=np.random.default_rng(147)
    for name,src,dst in [('reference_to_prepared',ref,m),('prepared_to_reference',m,ref)]:
        idx=rng.choice(len(src.faces),size=5000,p=src.area_faces/src.area)
        w=rng.random((len(idx),2));flip=w.sum(1)>1;w[flip]=1-w[flip]
        tri=src.triangles[idx];pts=tri[:,0]+w[:,:1]*(tri[:,1]-tri[:,0])+w[:,1:]*(tri[:,2]-tri[:,0])
        _,dist,_=trimesh.proximity.closest_point(dst,pts)
        report['distances'][name]=np.percentile(dist,[0,50,95,99,100]).tolist()
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,vertices=m.vertices,faces=m.faces)
    out.with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('BOUNDARY_PREPARED',out,report,flush=True)

if __name__=='__main__':main()
