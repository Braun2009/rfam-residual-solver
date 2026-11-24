#!/usr/bin/env python3
"""Tetrahedralize a Boolean fluid boundary and audit real cell connectivity.

Diagnostic mesh only: no claim of RFAM import, prism layers, or CFD quality.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np
import manifold3d as md
import tetgen
import meshio
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',required=True)
    parser.add_argument('--out',required=True)
    parser.add_argument('--simplify-mm',type=float,default=0.0)
    parser.add_argument('--quality',action='store_true')
    parser.add_argument('--min-ratio',type=float,default=1.5)
    parser.add_argument('--steiner-limit',type=int,default=500000)
    parser.add_argument('--max-volume-mm3',type=float,default=-1.0)
    parser.add_argument('--epsilon',type=float,default=1e-12)
    parser.add_argument('--input-tets',help='Audit an existing NPZ points/cells mesh without rerunning TetGen')
    parser.add_argument('--boundary',help='Boundary vertices/faces NPZ or an audited mesh points/boundary NPZ')
    parser.add_argument('--preserve-boundary',action='store_true',help='Do not insert Steiner points on the boundary')
    parser.add_argument('--relaxed-inlets',action='store_true',
                        help='Allow merged inlet apertures (coarse meshes cannot resolve every 2 mm hole)')
    args=parser.parse_args()
    root=Path(args.geometry);out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((root/'geometry_manifest.json').read_text())
    d=np.load(root/'fluid_boundary_mm.npz')
    original=md.Manifold(md.Mesh64(d['vertices'].astype(np.float64),d['faces'].astype(np.uint64)))
    if args.boundary:
        remesh=np.load(args.boundary)
        vertices=remesh['vertices'] if 'vertices' in remesh else remesh['points']
        triangles=remesh['faces'] if 'faces' in remesh else remesh['boundary']
        used=np.unique(triangles);mapping=np.full(len(vertices),-1,dtype=np.int64);mapping[used]=np.arange(len(used))
        surface=md.Manifold(md.Mesh64(vertices[used].astype(np.float64),mapping[triangles].astype(np.uint64)))
    else:
        surface=original.simplify(args.simplify_mm)
    if surface.status()!=md.Error.NoError:raise ValueError(surface.status())
    sm=surface.to_mesh64()
    points=np.asarray(sm.vert_properties)[:,:3].copy()
    faces=np.asarray(sm.tri_verts,dtype=np.int32)
    print('SURFACE',len(points),len(faces),'volume_mm3',surface.volume(),flush=True)
    np.savez_compressed(out/'mesher_boundary_mm.npz',vertices=points,faces=faces)
    start=time.time()
    if args.input_tets:
        raw=np.load(args.input_tets)
        nodes=raw['points'].copy();cells=raw['cells'].copy()
    else:
        t=tetgen.TetGen(points,faces)
        t.tetrahedralize(plc=True,quality=args.quality,minratio=args.min_ratio,
                        steinerleft=args.steiner_limit,neighout=True,quiet=False,verbose=0,
                        epsilon=args.epsilon,maxvolume=args.max_volume_mm3,
                        fixedvolume=args.max_volume_mm3>0,nobisect=args.preserve_boundary)
        nodes=t.node;cells=t.elem
        np.savez_compressed(out/'tetgen_raw_mm.npz',points=nodes,cells=cells)
    if nodes.ndim!=2 or nodes.shape[1]!=3 or not np.isfinite(nodes).all():
        raise ValueError('Invalid mesh point coordinates')
    if cells.ndim!=2 or cells.shape[1]!=4 or not np.issubdtype(cells.dtype,np.integer):
        raise ValueError('Expected integer linear-tetrahedron connectivity')
    if len(cells)==0 or np.any(cells<0) or np.any(cells>=len(nodes)):
        raise ValueError('Tetrahedral connectivity is empty or out of bounds')
    canonical=np.sort(cells,axis=1)
    if np.any(np.diff(canonical,axis=1)==0):raise ValueError('Repeated tetrahedron corner')
    if len(np.unique(canonical,axis=0))!=len(cells):raise ValueError('Duplicate tetrahedra')
    del canonical
    # Compact unused optimizer nodes; orphan vertices are not flow cells.
    used=np.unique(cells);removed_orphans=len(nodes)-len(used)
    mapping=np.full(len(nodes),-1,dtype=np.int64);mapping[used]=np.arange(len(used))
    cells=mapping[cells];nodes=nodes[used]
    print('TETRAHEDRA_READY','nodes',len(nodes),'tets',len(cells),'seconds',time.time()-start,flush=True)
    # Normalize orientation; negative input ordering is not silently treated
    # as a positive-volume overlap. Shared-face side tests follow below.
    xyz=nodes[cells]
    signed=np.einsum('ij,ij->i',xyz[:,1]-xyz[:,0],np.cross(xyz[:,2]-xyz[:,0],xyz[:,3]-xyz[:,0]))/6
    if not np.isfinite(signed).all() or np.any(signed==0):raise ValueError('Invalid or zero volume tetrahedra')
    negative=signed<0
    cells[negative,0],cells[negative,1]=cells[negative,1].copy(),cells[negative,0].copy()
    volumes=np.abs(signed)
    rel=abs(volumes.sum()-surface.volume())/surface.volume()
    # STL input coordinates have float32 precision (about 3e-5 mm at the
    # injector). Boundary recovery may collapse sub-resolution slivers.
    # Bound this input/output difference by one input-coordinate ULP times
    # boundary area; test exact closure against the recovered faces below.
    coordinate_ulp=float(np.spacing(np.float32(np.abs(points).max())))
    input_volume_error_bound=coordinate_ulp*surface.surface_area()
    if abs(volumes.sum()-surface.volume())>input_volume_error_bound:
        raise ValueError(f'Tetrahedral fill differs from input beyond STL precision: {rel}')
    # A shared face, not merely a touching point, defines cell connectivity.
    all_faces=cells[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
    keys=np.sort(all_faces,axis=1)
    _,first,inverse,counts=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
    if np.any(counts>2):raise ValueError('Nonmanifold tetrahedron face')
    order=np.argsort(inverse,kind='stable')
    starts=np.r_[0,np.cumsum(counts)[:-1]]
    pair_start=starts[counts==2]
    owners=order[pair_start]//4;neighbours=order[pair_start+1]//4
    # Owners and neighbours must be on opposite sides of every shared face.
    centres=nodes[cells].mean(axis=1)
    internal=all_faces[order[pair_start]]
    for lo in range(0,len(owners),200000):
        hi=min(lo+200000,len(owners));x=nodes[internal[lo:hi]]
        normal=np.cross(x[:,1]-x[:,0],x[:,2]-x[:,0]);fc=x.mean(axis=1)
        left=np.einsum('ij,ij->i',normal,fc-centres[owners[lo:hi]])
        right=np.einsum('ij,ij->i',normal,centres[neighbours[lo:hi]]-fc)
        if np.any(left<=0) or np.any(right<=0):raise ValueError('Inverted/overlapping cells across shared face')
    graph=coo_matrix((np.ones(len(owners),dtype=np.uint8),(owners,neighbours)),shape=(len(cells),len(cells))).tocsr()
    ncomponents,labels=connected_components(graph,directed=False)
    if ncomponents!=1:raise ValueError(f'Fluid cell domain has {ncomponents} disconnected components')
    boundary_indices=first[counts==1]
    boundary=all_faces[boundary_indices]
    bpoints=nodes[boundary]
    # Translation before summation improves cancellation for offset CAD data.
    shifted=bpoints-nodes.mean(axis=0)
    boundary_volume=float(np.einsum('ij,ij->i',shifted[:,0],np.cross(shifted[:,1],shifted[:,2])).sum()/6)
    closure=abs(volumes.sum()-boundary_volume)/volumes.sum()
    if closure>1e-10:raise ValueError(f'Recovered cell/boundary volume closure failed: {closure}')
    inlet=(np.max(np.abs(bpoints[:,:,2]-manifest['inlet_cad_z_mm']),axis=1)<1e-7)
    outlet=(np.max(np.abs(bpoints[:,:,2]-nodes[:,2].max()),axis=1)<1e-7)
    boundary_owner=boundary_indices//4
    if not np.any(inlet) or not np.any(outlet):raise ValueError('Missing inlet or outlet faces')
    # Count inlet patches by shared surface edges, preserving each aperture.
    inlet_faces=boundary[inlet]
    edge=np.sort(inlet_faces[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
    _,idx,inv,cnt=np.unique(edge,axis=0,return_index=True,return_inverse=True,return_counts=True)
    order=np.argsort(inv,kind='stable');starts=np.r_[0,np.cumsum(cnt)[:-1]];st=starts[cnt==2]
    a=order[st]//3;b=order[st+1]//3
    ng,group=connected_components(coo_matrix((np.ones(len(a)),(a,b)),shape=(len(inlet_faces),len(inlet_faces))).tocsr(),directed=False)
    ip=nodes[inlet_faces];areas=.5*np.linalg.norm(np.cross(ip[:,1]-ip[:,0],ip[:,2]-ip[:,0]),axis=1)
    ic=ip.mean(axis=1)
    inlet_groups=[]
    inlet_ids=np.zeros(len(boundary),dtype=np.int32)
    for g in range(ng):
        use=group==g;area=float(areas[use].sum());cent=(ic[use]*areas[use,None]).sum(axis=0)/area
        inlet_groups.append({'id':g+1,'area_mm2':area,'centroid_cad_mm':cent.tolist(),'boundary_triangles':int(use.sum()),'reaches_outlet':True})
    inlet_ids[inlet]=group+1
    if ng!=manifest.get('inlet_apertures') and not args.relaxed_inlets:
        raise ValueError(f'Inlet patch count changed: {ng} vs {manifest.get("inlet_apertures")}')
    # Area and location classify the deliberately truncated mass-flow inlets.
    # Physical names come from the design: 20 small propane feeds per pintle,
    # one main LOX passage per pintle, one central torch passage.
    axes=np.array([x['xy_mm'] for x in manifest['pintles']])
    centre=np.array(manifest['centre_xy_mm'])
    group_counts={name:[0,0,0] for name in ['fuel','lox']}
    torch_count=0
    boundary_zone=np.zeros(len(boundary),dtype=np.int32)  # 0 = wall
    boundary_zone[outlet]=1
    for g in inlet_groups:
        xy=np.array(g['centroid_cad_mm'][:2]);dist=np.linalg.norm(axes-xy,axis=1);p=int(np.argmin(dist))
        axis_distance=float(np.linalg.norm(xy-centre))
        if args.relaxed_inlets:
            # A coarse surface merges the 2 mm propane holes, so classify by
            # geometry instead of by exact patch area: the torch sits on the
            # domain axis, the LOX bore at a pintle centre, fuel farther out.
            if axis_distance<3.0:
                family='torch';zone=30;torch_count+=1
            elif dist[p]<5.0:
                family='lox';zone=10+p;group_counts['lox'][p]+=1
            else:
                family='fuel';zone=20+p;group_counts['fuel'][p]+=1
        elif 3.0<g['area_mm2']<3.3 and 5<dist[p]<10:
            family='fuel';zone=20+p;group_counts['fuel'][p]+=1
        elif 20<g['area_mm2']<70 and dist[p]<.2:
            family='lox';zone=10+p;group_counts['lox'][p]+=1
        elif np.linalg.norm(xy-centre)<.2 and 8<g['area_mm2']<20:
            family='torch';zone=30;torch_count+=1
        else:raise ValueError(f'Unclassified inlet aperture: {g}')
        g['design_family']=family;g['pintle_index']=p if family!='torch' else None
        g['boundary_zone']=zone
        boundary_zone[inlet_ids==g['id']]=zone
    if group_counts!={'fuel':[20,20,20],'lox':[1,1,1]} or torch_count!=1:
        if not args.relaxed_inlets:
            raise ValueError(f'Unexpected inlet distribution: {group_counts}, torch {torch_count}')
        print('RELAXED_INLETS distribution:',group_counts,'torch',torch_count,flush=True)
    edge2=np.sum((xyz[:,[0,0,0,1,1,2]]-xyz[:,[1,2,3,2,3,3]])**2,axis=(1,2))
    mean_ratio=12*np.power(3*volumes,2/3)/edge2
    report={
        'status':'TOPOLOGY_PASS_CFD_QUALITY_NOT_ACCEPTED', 'coordinate_system':'CAD Z downstream','units':'mm',
        'nodes':len(nodes),'tetrahedra':len(cells),'shared_faces':len(owners),'boundary_faces':len(boundary),
        'fluid_components_by_shared_faces':int(ncomponents),'nonmanifold_faces':int((counts>2).sum()),
        'duplicate_cells':0,'shared_face_side_check':True,'orphan_vertices_removed':removed_orphans,
        'input_negative_orderings_reoriented':int(negative.sum()),
        'volume_min_mm3':float(volumes.min()),'volume_sum_mm3':float(volumes.sum()),
        'surface_volume_mm3':surface.volume(),'fill_relative_volume_error':float(rel),
        'input_coordinate_ulp_mm':coordinate_ulp,'input_volume_error_bound_mm3':input_volume_error_bound,
        'recovered_boundary_volume_mm3':boundary_volume,'recovered_boundary_volume_closure':closure,
        'original_surface_volume_mm3':original.volume(),'surface_simplification_mm':args.simplify_mm,
        'surface_simplification_relative_volume_change':float((surface.volume()-original.volume())/original.volume()),
        'inlet_connected_patches':int(ng),'inlets':inlet_groups,
        'design_inlet_counts':group_counts,'torch_inlets':torch_count,
        'quality_mean_ratio_percentiles':np.percentile(mean_ratio,[0,.1,1,5,50,95,100]).tolist(),
        'quality_mean_ratio_below_0p01':int(np.sum(mean_ratio<.01)),
        'input_tetrahedra':args.input_tets,'input_boundary':args.boundary,
        'boundary_zones':{'0':'wall','1':'outlet','10':'lox_0','11':'lox_1','12':'lox_2',
                          '20':'fuel_0','21':'fuel_1','22':'fuel_2','30':'torch'},
        'generation_quality_refinement':args.quality,'min_ratio_requested':args.min_ratio,
        'steiner_limit':args.steiner_limit,'max_volume_mm3':args.max_volume_mm3,
        'preserve_boundary_requested':args.preserve_boundary,
        'elapsed_s':time.time()-start,
        'limitations':['No prism boundary layers or verified production spacing.',
                       'No RFAM unstructured mesh reader/face operators are implemented.',
                       'Assembly uses mating-plane reconstruction; inlet labels need geometry review.']}
    np.savez_compressed(out/'fluid_volume_mm.npz',points=nodes,cells=cells,boundary=boundary,
                        boundary_owner=boundary_owner,inlet_patch=inlet_ids,boundary_zone=boundary_zone,
                        internal_face_nodes=internal,internal_face_owner=owners,internal_face_neighbour=neighbours)
    meshio.write(out/'fluid_volume_mm.vtu',meshio.Mesh(nodes,[('tetra',cells)],cell_data={'fluid_component':[labels.astype(np.int32)],'volume_mm3':[volumes]}))
    meshio.write(out/'boundary_patches_mm.vtu',meshio.Mesh(nodes,[('triangle',boundary)],cell_data={'inlet_patch':[inlet_ids],'owner_cell':[boundary_owner],'boundary_zone':[boundary_zone]}))
    (out/'mesh_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('TOPOLOGY_AUDIT_PASS',len(cells),'inlet_patches',ng,'rel_volume_error',rel,flush=True)

if __name__=='__main__':main()
