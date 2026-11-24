#!/usr/bin/env python3
"""Sample independent channel components immediately upstream of their mergers.

An axial cut through an oblique LOX exit is NOT its discharge surface. Values
are cell-centred diagnostics; no post-exit mixing or ignition is inferred.
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


def components(points,cells,z):
    grid=cells_from_npz(points,cells)
    ids=numpy_to_vtk(np.arange(len(cells),dtype=np.int32),deep=True)
    ids.SetName('cell_id');grid.GetCellData().AddArray(ids)
    plane=vtk.vtkPlane();plane.SetOrigin(0,0,z);plane.SetNormal(0,0,1)
    cutter=vtk.vtkCutter();cutter.SetInputData(grid);cutter.SetCutFunction(plane);cutter.Update()
    conn=vtk.vtkPolyDataConnectivityFilter();conn.SetInputData(cutter.GetOutput())
    conn.SetExtractionModeToAllRegions();conn.ColorRegionsOn();conn.Update()
    poly=conn.GetOutput()
    if poly.GetNumberOfCells()==0:raise ValueError('plane misses mesh')
    origin=vtk_to_numpy(poly.GetCellData().GetArray('cell_id'))
    labels_at_points=vtk_to_numpy(poly.GetPointData().GetArray('RegionId'))
    vertices=vtk_to_numpy(poly.GetPoints().GetData())
    labels=np.empty(len(origin),dtype=np.int32)
    area=np.empty(len(origin),dtype=np.float64)
    centre=np.empty((len(origin),2),dtype=np.float64)
    for i in range(len(origin)):
        cell=poly.GetCell(i)
        idx=[cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        labels[i]=labels_at_points[idx[0]]
        xy=vertices[idx,:2]
        area[i]=.5*abs(np.dot(xy[:,0],np.roll(xy[:,1],-1))
                     -np.dot(xy[:,1],np.roll(xy[:,0],-1)))
        centre[i]=xy.mean(axis=0)
    n=conn.GetNumberOfExtractedRegions()
    region_area=np.bincount(labels,weights=area,minlength=n)
    centroid=np.column_stack([np.bincount(labels,weights=area*centre[:,k],minlength=n)
                              for k in range(2)])/region_area[:,None]
    return origin,labels,area,region_area,centroid


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--inlet-audit',required=True)
    ap.add_argument('--z-mm',type=float,required=True)
    ap.add_argument('--field',action='append',required=True,help='LABEL=FIELD')
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    points,cells,_=read_mesh(args.mesh)
    source,region,area,region_area,xy=components(points,cells,args.z_mm/1000)
    inlet=json.loads(Path(args.inlet_audit).read_text(encoding='utf-8'))['inlets']
    inlet_xy=np.array([x['centroid_cad_mm'][:2] for x in inlet])/1000
    largest=int(np.argmax(region_area))
    independent=np.arange(len(region_area))!=largest
    distance=np.linalg.norm(xy[:,None,:]-inlet_xy[None,:,:],axis=2)
    nearest=distance.argmin(axis=1)
    family={int(k):inlet[int(nearest[k])]['design_family'] for k in np.flatnonzero(independent)}
    expected={271.5:{'fuel':60,'lox':9,'torch':1},
              271.8:{'fuel':60,'lox':9,'torch':1},
              276.6:{'fuel':60,'lox':9,'torch':1},
              276.8:{'fuel':60,'lox':9,'torch':1},
              277.68:{'fuel':30,'lox':3,'torch':1},
              280.8:{'fuel':0,'lox':3,'torch':1},
              280.9:{'fuel':0,'lox':3,'torch':1},
              283.3:{'fuel':0,'lox':0,'torch':1}}
    lookup=expected.get(round(args.z_mm,2))
    counts={name:sum(f==name for f in family.values()) for name in ('fuel','lox','torch')}
    if lookup is None or counts!=lookup:
        raise ValueError(f'component families at z={args.z_mm} do not match the audited topology: {counts}')
    report={'kind':'independent internal channel cut, NOT injector outlet',
            'z_mm':args.z_mm,'region_count':len(region_area),
            'external_area_mm2':float(region_area[largest]*1e6),
            'internal_count_by_family':counts,'fields':{}}
    for spec in args.field:
        label,path=spec.split('=',1)
        _,data,y,ns=read_field(path)
        if len(data)!=len(cells) or ns!=2:raise ValueError('field/mesh mismatch')
        items={}
        for name in ('fuel','lox','torch'):
            members=np.array([k for k,v in family.items() if v==name],dtype=int)
            sel=np.isin(region,members)
            a=area[sel];indices=source[sel]
            if not len(a):
                items[name]=dict(components=0,area_mm2=0.,Y_O2_min=None,Y_O2_max=None,
                                 Y_O2_mean=None,rho_mean=None,T_mean=None,axial_mass_kg_s=0.)
                continue
            weights=a/a.sum()
            composition=y[indices,0].astype(float)
            density=data[indices,5].astype(float)
            items[name]=dict(components=len(members),area_mm2=float(a.sum()*1e6),
                Y_O2_min=float(composition.min()),Y_O2_max=float(composition.max()),
                Y_O2_mean=float(np.sum(composition*weights)),
                rho_mean=float(np.sum(density*weights)),
                T_mean=float(np.sum(data[indices,1]*weights)),
                axial_mass_kg_s=float(np.sum(a*density*data[indices,4])))
        external=region==largest
        ext_y=y[source[external],0].astype(float)
        ext_area=area[external]
        outside=dict(area_mm2=float(ext_area.sum()*1e6),
            Y_O2_min=float(ext_y.min()),Y_O2_max=float(ext_y.max()),
            nonbackground_area_frac=float(ext_area[np.abs(ext_y-.784)>.1].sum()/ext_area.sum()),
            pure_O2_polygon_count=int(np.count_nonzero(ext_y>.95)),
            fuel_rich_polygon_count=int(np.count_nonzero(ext_y<.2)))
        report['fields'][label]=dict(source=path,families=items,external=outside)
    Path(args.out).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for label,content in report['fields'].items():
        print(label,'z_mm',args.z_mm,'regions',len(region_area),
              [(name,d['components'],d['Y_O2_mean'],d['axial_mass_kg_s'])
               for name,d in content['families'].items()],
              'external_Y_range',(content['external']['Y_O2_min'],
                                  content['external']['Y_O2_max']))
    print('PREEXIT_COMPONENTS_DIAGNOSTIC_WRITTEN',args.out)


if __name__=='__main__':main()
