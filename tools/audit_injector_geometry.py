#!/usr/bin/env python3
"""Independent sampled geometry/region audit of a tetrahedral candidate.

Connectivity and exact recovered-boundary volume closure are checked by
mesh_injector_tetra.py. This tool reports sampled geometric error, cell sizes
and normals. Sampled geometry results are not a Hausdorff or intersection proof.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import trimesh
from scipy.optimize import linear_sum_assignment


def to_mesh(npz):
    d=np.load(npz)
    v=d['vertices'] if 'vertices' in d else d['points']
    f=d['faces'] if 'faces' in d else d['boundary']
    used=np.unique(f);ids=np.full(len(v),-1,dtype=np.int64);ids[used]=np.arange(len(used))
    return trimesh.Trimesh(v[used],ids[f],process=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',required=True);parser.add_argument('--geometry',required=True)
    parser.add_argument('--reference-mesh',required=True);parser.add_argument('--out',required=True)
    a=parser.parse_args();root=Path(a.mesh);geom=Path(a.geometry)
    d=np.load(root/'fluid_volume_mm.npz');p=d['points'];c=d['cells'];b=d['boundary']
    manifest=json.loads((geom/'geometry_manifest.json').read_text())
    target=json.loads((root/'mesh_audit.json').read_text())
    reference=json.loads((Path(a.reference_mesh)/'mesh_audit.json').read_text())
    # Check every physical boundary normal against its owner cell centre.
    x=p[c];centres=x.mean(1)
    normals=np.cross(p[b[:,1]]-p[b[:,0]],p[b[:,2]]-p[b[:,0]])
    face_c=p[b].mean(1);dot=np.einsum('ij,ij->i',normals,face_c-centres[d['boundary_owner']])
    if np.any(dot<=0):raise ValueError('Boundary face normal does not point out of its owner cell')
    before=to_mesh(geom/'fluid_boundary_mm.npz');after=to_mesh(root/'fluid_volume_mm.npz')
    rng=np.random.default_rng(514)
    if not after.is_watertight or not after.is_winding_consistent:
        raise ValueError('Exported tetrahedral boundary is not closed and oriented')
    report={'sampling_is_not_hausdorff_or_intersection_proof':True,
            'closed_boundary':bool(after.is_watertight),'euler':int(after.euler_number),
            'outward_boundary_faces':len(b),'wrong_boundary_normals':0,
            'volume_change_from_exact_boolean_relative':(after.volume-before.volume)/before.volume,
            'surface_distance_samples':{},'regional_edges':{},'inlet_area_changes':[]}
    for name,src,dst in [('cad_to_mesh',before,after),('mesh_to_cad',after,before)]:
        ids=rng.choice(len(src.faces),size=8000,p=src.area_faces/src.area)
        w=rng.random((len(ids),2));flip=w.sum(1)>1;w[flip]=1-w[flip]
        tri=src.triangles[ids];pts=tri[:,0]+w[:,:1]*(tri[:,1]-tri[:,0])+w[:,1:]*(tri[:,2]-tri[:,0])
        _,distance,_=trimesh.proximity.closest_point(dst,pts)
        if not np.isfinite(distance).all():
            raise ValueError('Surface-distance query produced nonfinite values')
        report['surface_distance_samples'][name]={'count':len(pts),'distance_mm_percentiles':np.percentile(distance,[0,50,95,99,100]).tolist()}
        print(name,report['surface_distance_samples'][name],flush=True)
    # Compare each aperture with the unsmoothed CAD-derived mesh by position.
    old=reference['inlets'];new=target['inlets']
    dist=np.linalg.norm(np.array([i['centroid_cad_mm'] for i in old])[:,None,:]-np.array([i['centroid_cad_mm'] for i in new])[None,:,:],axis=2)
    i,j=linear_sum_assignment(dist)
    for oi,ni in zip(i,j):
        o,n=old[oi],new[ni]
        report['inlet_area_changes'].append({'candidate_id':n['id'],'family':n['design_family'],
            'area_reference_mm2':o['area_mm2'],'area_candidate_mm2':n['area_mm2'],
            'relative_area_change':n['area_mm2']/o['area_mm2']-1,
            'centroid_distance_mm':float(dist[oi,ni])})
    lengths=np.linalg.norm(x[:,[0,0,0,1,1,2]]-x[:,[1,2,3,2,3,3]],axis=2)
    axial=centres[:,2]-manifest['board_cad_z_mm']
    for name,lo,hi in [('nozzle',-2,15),('near_field',15,30),('chamber',30,254),('throat',254,286),('downstream_refinement',286,310),('nozzle_exit',310,406.5)]:
        use=(axial>=lo)&(axial<hi)
        report['regional_edges'][name]={'cells':int(use.sum()),'edge_length_mm_p01_p50_p99':np.percentile(lengths[use],[1,50,99]).tolist()}
    solid=to_mesh(geom/'assembled_head_mm.npz')
    # Probe the most at-risk near-head cells; allow only the measured boundary
    # approximation, reporting intrusions rather than reclassifying the CAD.
    use=np.flatnonzero(axial<17)
    sample=rng.choice(use,size=min(5000,len(use)),replace=False)
    signed=trimesh.proximity.signed_distance(solid,centres[sample])
    if not np.isfinite(signed).all():
        raise ValueError('Solid-distance query produced nonfinite values')
    report['sampled_cell_centres_vs_metal']={'count':len(sample),'inside_metal':int(np.sum(signed>0)),
            'deeper_than_0p04_mm':int(np.sum(signed>.04)),'maximum_signed_depth_mm':float(signed.max())}
    if np.any(signed>.04):raise ValueError('Sampled fluid-cell centre penetrates metal beyond measured remesh scale')
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('GEOMETRY_AUDIT_WRITTEN',out,flush=True)

if __name__=='__main__':main()
