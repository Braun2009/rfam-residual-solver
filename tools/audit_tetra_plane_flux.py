#!/usr/bin/env python3
"""Integrate cell-centred fluxes on a VTK tetra/plane cut (diagnostic only).

The cut is a fixed axial slice, not the geometric exit surfaces of the LOX
slots or 60 discrete fuel orifices. No unsteady mass is integrated in time.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from inspect_injector_volume import cells_from_npz
from tetra_field_to_vtu import read_mesh,read_field


def cut_area_by_cell(points,cells,z):
    grid=cells_from_npz(points,cells)
    ids=numpy_to_vtk(np.arange(len(cells),dtype=np.int32),deep=True)
    ids.SetName('cell_id');grid.GetCellData().AddArray(ids)
    plane=vtk.vtkPlane();plane.SetOrigin(0,0,z);plane.SetNormal(0,0,1)
    cutter=vtk.vtkCutter();cutter.SetInputData(grid);cutter.SetCutFunction(plane);cutter.Update()
    cut=cutter.GetOutput()
    if cut.GetNumberOfCells()==0:
        raise ValueError('plane does not intersect the mesh')
    idata=vtk_to_numpy(cut.GetCellData().GetArray('cell_id'))
    vertices=vtk_to_numpy(cut.GetPoints().GetData())
    if not len(idata):raise ValueError('plane does not intersect the mesh')
    area=np.empty(len(idata),dtype=np.float64)
    for i in range(len(idata)):
        cell=cut.GetCell(i)
        p=vertices[[cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())],:2]
        area[i]=abs(np.dot(p[:,0],np.roll(p[:,1],-1))
                     -np.dot(p[:,1],np.roll(p[:,0],-1)))*.5
    return idata,area


def audit(args):
    points,cells,_=read_mesh(args.mesh)
    idata,area=cut_area_by_cell(points,cells,args.z)
    summary={'mesh':str(args.mesh),'z_m':args.z,'polygons':len(idata),
             'total_area_mm2':float(area.sum()*1e6),
             'method':'VTK tetra cut, cell-centred rho*u_z*A; not CAD exits',
             'converged':False,'fields':{}}
    for field in args.field:
        label,path=field.split('=',1)
        cent,f,y,ns=read_field(path)
        if len(cells)!=len(f) or ns!=2:raise ValueError('mesh/field mismatch '+label)
        mass=area*f[idata,5].astype(float)*f[idata,4].astype(float)
        oxygen=mass*y[idata,0].astype(float)
        fuel=mass*y[idata,1].astype(float)
        summary['fields'][label]={'source':path,
            'net_mass_kg_s':float(mass.sum()),
            'net_O2_kg_s':float(oxygen.sum()),
            'net_fuel_kg_s':float(fuel.sum()),
            'positive_mass_kg_s':float(mass[mass>0].sum()),
            'negative_mass_kg_s':float(mass[mass<0].sum()),
            'composition_closure_kg_s':float((mass-oxygen-fuel).sum())}
    return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--z',type=float,required=True,help='SI metres')
    ap.add_argument('--field',action='append',required=True,help='LABEL=field.rfamtetf')
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    summary=audit(args)
    Path(args.out).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for label,data in summary['fields'].items():
        print(label,'z',args.z,'area_mm2',summary['total_area_mm2'],
              'mass',data['net_mass_kg_s'],'O2',data['net_O2_kg_s'],
              'fuel',data['net_fuel_kg_s'])
    print('TETRA_PLANE_FLUX_DIAGNOSTIC_WRITTEN',args.out)


if __name__=='__main__':main()
