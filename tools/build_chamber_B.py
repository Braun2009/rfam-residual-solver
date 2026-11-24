#!/usr/bin/env python3
"""Route B equivalent-exit plane: chamber/nozzle only, no injector passages.

CAD axial Z and mm. Equivalent injection rings are a declared model, not CAD
extracted outlets. The remaining inlet-plane area is an impermeable plate.
"""
import argparse,json,math,struct,hashlib
from pathlib import Path
import numpy as np
import gmsh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile',required=True);p.add_argument('--out',required=True)
    p.add_argument('--h-mm',type=float,default=6.)
    p.add_argument('--jet-mm',type=float,default=1.0)
    a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    doc=json.loads(Path(a.profile).read_text());prof=np.array([[q['x_mm'],q['r_mm']] for q in doc['profile']])
    # Use actual design profile knots at important transitions.
    zboard=270.9527893066406;cut=14.;axis=np.array([49.215348,48.350049])
    sample=np.unique(np.r_[cut,np.arange(20,406,10),123,273,276.5,406.5])
    radius=np.interp(sample,prof[:,0],prof[:,1])
    gmsh.initialize();gmsh.option.setNumber('General.Terminal',1)
    try:
        gmsh.model.add('B_equivalent_injection');o=gmsh.model.occ
        pts=[o.addPoint(0,0,zboard+cut)]+[o.addPoint(r,0,zboard+z) for z,r in zip(sample,radius)]+[o.addPoint(0,0,zboard+sample[-1])]
        curves=[o.addLine(pts[i],pts[(i+1)%len(pts)]) for i in range(len(pts))]
        surface=o.addPlaneSurface([o.addCurveLoop(curves)])
        rev=o.revolve([(2,surface)],0,0,0,0,0,1,2*math.pi)
        volume=[x for x in rev if x[0]==3]
        z=zboard+cut
        rings=[];spec=[]
        for i,ang in enumerate([210,330,90]):
            cx=12.5*math.cos(math.radians(ang));cy=12.5*math.sin(math.radians(ang))
            for family,zone,rin,area in [('lox',10+i,1.5,20.990),('fuel',20+i,5.5,59.922)]:
                rout=math.sqrt(rin**2+area/math.pi)
                outer=o.addDisk(cx,cy,z,rout,rout);inner=o.addDisk(cx,cy,z,rin,rin)
                ann=o.cut([(2,outer)],[(2,inner)])[0]
                rings.extend(ann);spec.append({'zone':zone,'family':family,'xy_from_axis_mm':[cx,cy],'r_inner_mm':rin,'r_outer_mm':rout,'target_area_mm2':area})
        torch=o.addDisk(0,0,z,1.7,1.7);rings.append((2,torch))
        spec.append({'zone':30,'family':'torch','xy_from_axis_mm':[0,0],'r_inner_mm':0,'r_outer_mm':1.7,'target_area_mm2':math.pi*1.7**2})
        o.fragment(volume,rings);o.synchronize()
        # Remove standalone imprint helpers; volume boundary is the authority.
        vols=gmsh.model.getEntities(3)
        if len(vols)!=1:raise ValueError('Expected one chamber volume')
        bnd=gmsh.model.getBoundary(vols,oriented=False)
        groups={};tagzone={}
        for dim,tag in bnd:
            bb=o.getBoundingBox(dim,tag);cent=o.getCenterOfMass(dim,tag)
            zone=0
            if abs(bb[2]-z)<1e-5 and abs(bb[5]-z)<1e-5:
                area=o.getMass(dim,tag)
                for s in spec:
                    if abs(area-s['target_area_mm2'])<.01 and np.linalg.norm(np.array(cent[:2])-s['xy_from_axis_mm'])<.01:zone=s['zone'];break
            elif abs(bb[2]-(zboard+406.5))<1e-5 and abs(bb[5]-(zboard+406.5))<1e-5:zone=1
            groups.setdefault(zone,[]).append(tag);tagzone[tag]=zone
        if set(groups)!={0,1,10,11,12,20,21,22,30}:raise ValueError(f'Missing boundary groups: {groups.keys()}')
        for zone,tags in groups.items():gmsh.model.addPhysicalGroup(2,tags,tag=zone+100,name=f'zone_{zone}')
        gmsh.model.addPhysicalGroup(3,[vols[0][1]],tag=1,name='fluid')
        field=gmsh.model.mesh.field
        jet=field.add('Box')
        for k,v in dict(VIn=a.jet_mm,VOut=a.h_mm,XMin=-24,XMax=24,YMin=-24,YMax=24,ZMin=z-.1,ZMax=z+12,Thickness=15).items():field.setNumber(jet,k,v)
        throat=field.add('Box')
        for k,v in dict(VIn=a.h_mm*.5,VOut=a.h_mm,XMin=-30,XMax=30,YMin=-30,YMax=30,ZMin=zboard+255,ZMax=zboard+293,Thickness=12).items():field.setNumber(throat,k,v)
        mn=field.add('Min');field.setNumbers(mn,'FieldsList',[jet,throat]);field.setAsBackgroundMesh(mn)
        gmsh.option.setNumber('Mesh.MeshSizeFromPoints',0);gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0);gmsh.option.setNumber('Mesh.Algorithm3D',1)
        gmsh.option.setNumber('General.NumThreads',8)
        gmsh.model.mesh.generate(3);gmsh.model.mesh.optimize('Netgen')
        ids,xyz,_=gmsh.model.mesh.getNodes();points=np.array(xyz).reshape(-1,3);points[:,:2]+=axis
        remap=np.full(int(max(ids))+1,-1,dtype=np.int64);remap[ids]=np.arange(len(ids))
        _,conn=gmsh.model.mesh.getElementsByType(4);cells=remap[np.array(conn).reshape(-1,4)]
        t=points[cells];v=np.einsum('ij,ij->i',t[:,1]-t[:,0],np.cross(t[:,2]-t[:,0],t[:,3]-t[:,0]))/6
        neg=v<0;cells[neg,0],cells[neg,1]=cells[neg,1].copy(),cells[neg,0].copy()
        boundary=[];zones=[]
        for tag,zone in tagzone.items():
            _,c=gmsh.model.mesh.getElementsByType(2,tag);f=remap[np.array(c).reshape(-1,3)]
            boundary.append(f);zones.extend([zone]*len(f))
        boundary=np.vstack(boundary);zones=np.array(zones,dtype=np.int32)
        allfaces=cells[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
        key,first,inv,count=np.unique(np.sort(allfaces,axis=1),axis=0,return_index=True,return_inverse=True,return_counts=True)
        if np.any(count>2):raise ValueError('Nonmanifold faces')
        lookup={tuple(sorted(f)):z for f,z in zip(boundary,zones)}
        ext=first[count==1];boundary=allfaces[ext];zones=np.array([lookup[tuple(sorted(f))] for f in boundary],dtype=np.int32)
        order=np.argsort(inv,kind='stable');start=np.r_[0,np.cumsum(count)[:-1]];pair=start[count==2]
        own=order[pair]//4;nei=order[pair+1]//4
        n,_=connected_components(coo_matrix((np.ones(len(own)),(own,nei)),shape=(len(cells),len(cells))),directed=False)
        if n!=1:raise ValueError('Disconnected chamber')
        np.savez_compressed(out/'fluid_volume_mm.npz',points=points,cells=cells,boundary=boundary,boundary_zone=zones,boundary_owner=ext//4)
        areas=np.linalg.norm(np.cross(points[boundary[:,1]]-points[boundary[:,0]],points[boundary[:,2]]-points[boundary[:,0]]),axis=1)/2
        for s in spec:s['mesh_area_mm2']=float(areas[zones==s['zone']].sum())
        manifest={'route':'B_equivalent_exit_plane','cad_board_z_mm':zboard,'inlet_plane_from_board_mm':cut,'axis_xy_mm':axis.tolist(),
                  'pintle_centres_from_axis_mm':[s['xy_from_axis_mm'] for s in spec if s['family']=='lox'],
                  'patches':spec,'cells':len(cells),'points':len(points),'components':n,'min_volume_mm3':float(abs(v).min()),
                  'volume_mm3':float(abs(v).sum()),'geometry':'axisymmetric chamber profile downstream of 14 mm; equivalent annular inlets; remainder of cut plane is wall',
                  'limits':['No internal injector passages','Equivalent plate excludes upstream recirculation','No ignition/performance conclusion','No prism boundary layers']}
        (out/'B_manifest.json').write_text(json.dumps(manifest,indent=2))
        print('B_MESH_PASS',len(cells),'cells',len(boundary),'boundary faces',flush=True)
    finally:gmsh.finalize()

if __name__=='__main__':main()
