#!/usr/bin/env python3
"""Audit 64 distinct upstream injector passages on an actual tetra slice.

The slice must lie BEFORE passage components split or merge. This is a
cell-centred diagnostic, not integration over CAD injector outlet surfaces.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from tetra_field_to_vtu import read_mesh,read_field
from inspect_injector_volume import cells_from_npz


def passage_cut(points,cells,z):
    grid=cells_from_npz(points,cells)
    ids=numpy_to_vtk(np.arange(len(cells),dtype=np.int32),deep=True)
    ids.SetName('cell_id');grid.GetCellData().AddArray(ids)
    plane=vtk.vtkPlane();plane.SetOrigin(0,0,z);plane.SetNormal(0,0,1)
    cutter=vtk.vtkCutter();cutter.SetInputData(grid);cutter.SetCutFunction(plane);cutter.Update()
    conn=vtk.vtkPolyDataConnectivityFilter();conn.SetInputData(cutter.GetOutput())
    conn.SetExtractionModeToAllRegions();conn.ColorRegionsOn();conn.Update()
    result=conn.GetOutput();count=conn.GetNumberOfExtractedRegions()
    if count!=64:raise ValueError(f'upstream slice requires exactly 64 passage regions, found {count}')
    origin=vtk_to_numpy(result.GetCellData().GetArray('cell_id'))
    region_point=vtk_to_numpy(result.GetPointData().GetArray('RegionId'))
    verts=vtk_to_numpy(result.GetPoints().GetData())
    labels=np.empty(len(origin),dtype=np.int32);area=np.empty(len(origin));centre=np.empty((len(origin),2))
    for i in range(len(origin)):
        poly=result.GetCell(i)
        idx=[poly.GetPointId(j) for j in range(poly.GetNumberOfPoints())]
        labels[i]=region_point[idx[0]]
        xy=verts[idx,:2]
        area[i]=.5*abs(np.dot(xy[:,0],np.roll(xy[:,1],-1))
                     -np.dot(xy[:,1],np.roll(xy[:,0],-1)))
        centre[i]=np.mean(xy,axis=0)
    regions=np.bincount(labels,weights=area,minlength=count)
    xy=np.column_stack([np.bincount(labels,weights=area*centre[:,axis],minlength=count)
                         for axis in (0,1)])/regions[:,None]
    return origin,labels,area,regions,xy


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--inlet-audit',required=True,help='delivery/reports/mesh_audit.json')
    ap.add_argument('--z-mm',type=float,default=270.5)
    ap.add_argument('--field',action='append',required=True,help='LABEL=field.rfamtetf')
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    points,cells,_=read_mesh(args.mesh)
    origin,labels,areas,region_areas,xy=passage_cut(points,cells,args.z_mm/1000)
    inlet=json.loads(Path(args.inlet_audit).read_text(encoding='utf-8'))['inlets']
    if len(inlet)!=64:raise ValueError('inlet audit does not contain 64 patches')
    inlet_xy=np.asarray([row['centroid_cad_mm'][:2] for row in inlet])*1e-3
    distances=np.linalg.norm(xy[:,None,:]-inlet_xy[None,:,:],axis=2)
    assignment=np.argmin(distances,axis=1)
    if len(set(assignment))!=64 or distances[np.arange(64),assignment].max()>0.0001:
        raise ValueError('section-to-inlet mapping is ambiguous')
    passages=[dict(inlet_patch=int(inlet[j]['id']),family=inlet[j]['design_family'],
                    pintle_index=inlet[j].get('pintle_index'),region=int(i),
                    area_mm2=float(region_areas[i]*1e6),
                    inlet_centroid_distance_mm=float(distances[i,j]*1000))
              for i,j in enumerate(assignment)]
    report={'mesh':args.mesh,'inlet_audit':args.inlet_audit,'z_mm':args.z_mm,
            'kind':'upstream distinct passages, NOT geometric injector exits',
            'area_mm2':float(areas.sum()*1e6),'max_match_mm':max(row['inlet_centroid_distance_mm'] for row in passages),
            'passages':passages,'fields':{}}
    family_names=['lox','fuel','torch']
    for spec in args.field:
        label,path=spec.split('=',1)
        cent,f,y,ns=read_field(path)
        if len(f)!=len(cells) or ns!=2:raise ValueError('field/mesh mismatch '+label)
        mass=areas*f[origin,5].astype(float)*f[origin,4].astype(float)
        oxygen=mass*y[origin,0].astype(float);fuel=mass*y[origin,1].astype(float)
        m=np.bincount(labels,weights=mass,minlength=64)
        o=np.bincount(labels,weights=oxygen,minlength=64)
        b=np.bincount(labels,weights=fuel,minlength=64)
        split=[]
        for row in passages:
            i=row['region'];split.append(dict(inlet_patch=row['inlet_patch'],family=row['family'],
                pintle_index=row['pintle_index'],mass_kg_s=float(m[i]),
                O2_kg_s=float(o[i]),fuel_kg_s=float(b[i])))
        grouped={}
        for family in family_names:
            sel=[r for r in split if r['family']==family]
            if len(sel)!={'lox':3,'fuel':60,'torch':1}[family]:
                raise ValueError('incorrect passage-family count')
            grouped[family]=dict(count=len(sel),
                area_mm2=sum(row['area_mm2'] for row in passages if row['family']==family),
                mass_kg_s=sum(r['mass_kg_s'] for r in sel),
                O2_kg_s=sum(r['O2_kg_s'] for r in sel),fuel_kg_s=sum(r['fuel_kg_s'] for r in sel),
                min_patch_kg_s=min(r['mass_kg_s'] for r in sel),
                max_patch_kg_s=max(r['mass_kg_s'] for r in sel))
        report['fields'][label]=dict(source=path,groups=grouped,per_passage=split,
            mass_closure_kg_s=float(np.sum(m-o-b)),
            note='cell-centred rho*u_z*A; axial section in channels before outlets')
    Path(args.out).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for label,field in report['fields'].items():
        print(label,'z_mm',args.z_mm,'area_mm2',report['area_mm2'],
              'LOX',field['groups']['lox']['mass_kg_s'],
              'fuel',field['groups']['fuel']['mass_kg_s'],
              'torch',field['groups']['torch']['mass_kg_s'])
    print('INJECTOR_PASSAGE_CUT_DIAGNOSTIC_WRITTEN',args.out)


if __name__=='__main__':main()
