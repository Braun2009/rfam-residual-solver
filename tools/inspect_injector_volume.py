#!/usr/bin/env python3
"""Render and export genuine tetrahedral cuts, never floating interior planes."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection


def cells_from_npz(points,cells):
    p=vtk.vtkPoints();p.SetData(numpy_to_vtk(np.asarray(points,dtype=np.float64),deep=True))
    ca=vtk.vtkCellArray()
    ca.SetData(numpy_to_vtkIdTypeArray(np.arange(0,len(cells)*4+1,4,dtype=np.int64),deep=True),
               numpy_to_vtkIdTypeArray(cells.astype(np.int64).ravel(),deep=True))
    g=vtk.vtkUnstructuredGrid();g.SetPoints(p);g.SetCells(vtk.VTK_TETRA,ca)
    return g


def triangulate_surface(poly):
    t=vtk.vtkTriangleFilter();t.SetInputData(poly);t.Update()
    return t.GetOutput()


def save_stl(poly,path):
    w=vtk.vtkSTLWriter();w.SetFileName(str(path));w.SetInputData(triangulate_surface(poly));w.SetFileTypeToBinary()
    if not w.Write():raise IOError(path)


def clip_grid(grid,origin,normal):
    plane=vtk.vtkPlane();plane.SetOrigin(*origin);plane.SetNormal(*normal)
    c=vtk.vtkClipDataSet();c.SetInputData(grid);c.SetClipFunction(plane);c.SetValue(0);c.Update()
    return c.GetOutput()


def volume_cut(grid,origin,normal):
    s=vtk.vtkDataSetSurfaceFilter();s.SetInputData(clip_grid(grid,origin,normal));s.Update()
    # The new cut faces come from clipped volume cells. No caps are invented
    # without cells beneath them, and no original interior faces are exported.
    return s.GetOutput()


def slice_polygons(grid,origin,normal,axes,offset):
    plane=vtk.vtkPlane();plane.SetOrigin(*origin);plane.SetNormal(*normal)
    c=vtk.vtkCutter();c.SetInputData(grid);c.SetCutFunction(plane);c.Update()
    pd=c.GetOutput();pts=vtk_to_numpy(pd.GetPoints().GetData())
    polys=[]
    for i in range(pd.GetNumberOfCells()):
        cell=pd.GetCell(i)
        if cell.GetNumberOfPoints()<3:continue
        coords=pts[[cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]][:,axes]-offset
        polys.append(coords)
    return polys


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--mesh',required=True);a.add_argument('--geometry',required=True);a.add_argument('--out',required=True)
    args=a.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    data=np.load(Path(args.mesh)/'fluid_volume_mm.npz');points=data['points'];cells=data['cells']
    manifest=json.loads((Path(args.geometry)/'geometry_manifest.json').read_text())
    centre=np.asarray(manifest['centre_xy_mm']);board=manifest['board_cad_z_mm']
    grid=cells_from_npz(points,cells)
    # Full-length half, and a nozzle close-up. Both are cut from actual cells.
    half=volume_cut(grid,[centre[0],0,0],[1,0,0])
    save_stl(half,out/'fluid_tetra_half_mm.stl')
    local=clip_grid(grid,[0,0,board+30],[0,0,-1])
    localhalf=volume_cut(local,[centre[0],0,0],[1,0,0])
    save_stl(localhalf,out/'nozzle_tetra_half_mm.stl')
    # Presentation frame only: axial -> X, transverse -> Y, cut normal -> +Z.
    # This lets BCAD light and frame the actual cut face without editing CAD.
    matrix=vtk.vtkMatrix4x4();matrix.Identity()
    transform=np.array([[0,0,1,-board],[0,1,0,-centre[1]],[-1,0,0,centre[0]],[0,0,0,1]])
    for i in range(4):
        for j in range(4):matrix.SetElement(i,j,transform[i,j])
    vt=vtk.vtkTransform();vt.SetMatrix(matrix)
    tf=vtk.vtkTransformPolyDataFilter();tf.SetInputData(localhalf);tf.SetTransform(vt);tf.Update()
    save_stl(tf.GetOutput(),out/'nozzle_cut_view_mm.stl')
    (out/'view_transform.json').write_text(json.dumps({'cad_to_view_4x4':transform.tolist(),'units':'mm'},indent=2))
    fig,axs=plt.subplots(2,3,figsize=(15,10),constrained_layout=True)
    for ax,z in zip(axs.ravel(),[manifest['inlet_cad_z_mm']+.01,273,276.5,279.5,282,287]):
        polys=slice_polygons(grid,[0,0,z],[0,0,1],[0,1],centre)
        ax.add_collection(PolyCollection(polys,facecolors='#72becb',edgecolors='#244651',linewidths=.18))
        width=max(25,float(np.linalg.norm(points[np.abs(points[:,2]-z)<.5,:2]-centre,axis=1).max()) if np.any(np.abs(points[:,2]-z)<.5) else 25)
        ax.set(xlim=(-width,width),ylim=(-width,width),xlabel='X from axis [mm]',ylabel='Y from axis [mm]',title=f'CAD z={z:.2f} mm: {len(polys)} cut cell polygons')
        ax.set_aspect('equal')
    fig.suptitle('Actual tetrahedral cell cuts: blue = fluid cells; white = excluded solid')
    fig.savefig(out/'tetra_cross_sections.png',dpi=170);plt.close(fig)
    # Longitudinal through central axis and the upper pintle.
    polys=slice_polygons(grid,[centre[0],0,0],[1,0,0],[2,1],np.array([board,centre[1]]))
    fig,axs=plt.subplots(2,1,figsize=(14,10),constrained_layout=True)
    for ax in axs:
        ax.add_collection(PolyCollection(polys,facecolors='#72becb',edgecolors='#244651',linewidths=.2))
        ax.set_aspect('equal');ax.set(xlabel='axial distance from board [mm]',ylabel='transverse distance [mm]')
    axs[0].set(xlim=(-2,408),ylim=(-67,67),title='Entire chamber: slice through volume cells')
    axs[1].set(xlim=(-2,30),ylim=(-5,24),title='Nozzle and torch: actual fluid passages, blank metal')
    fig.savefig(out/'tetra_longitudinal_sections.png',dpi=170);plt.close(fig)
    # Dimensionless tetrahedron mean ratio: 1=regular, tends to 0 for slivers.
    v=points[cells];vol=np.abs(np.einsum('ij,ij->i',v[:,1]-v[:,0],np.cross(v[:,2]-v[:,0],v[:,3]-v[:,0])))/6
    edge2=np.sum((v[:,[0,0,0,1,1,2]]-v[:,[1,2,3,2,3,3]])**2,axis=(1,2))
    quality=12*np.power(3*vol,2/3)/edge2
    report={'real_cell_count':len(cells),'cut_from_volume_cells':True,
            'mean_ratio_percentiles':np.percentile(quality,[0,.1,1,5,50,95,100]).tolist(),
            'mean_ratio_below_0p01':int(np.sum(quality<.01)),
            'cell_volume_min_mm3':float(vol.min()),
            'production_quality_accepted':False,
            'inspection_stl_units':'mm','inspection_stl_is_a_surface_of_clipped_cells':True}
    (out/'inspection_report.json').write_text(json.dumps(report,indent=2))
    print('CELL_CUTS_EXPORTED',out,report,flush=True)

if __name__=='__main__':main()
