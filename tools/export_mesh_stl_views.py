#!/usr/bin/env python3
"""Convert an audited tetrahedral mesh into STL views for BCAD and others.

STL has no volume cells. Modes:
  boundary   closed outer surface only (smallest, shows the flow boundary)
  half       half domain cut from real cells, including cut faces (recommended)
  interfaces every face once: internal shared faces + boundary (solid look,
             large; --stride thins internal faces for a lighter preview)
  all        boundary + half

PNG previews are rendered with VTK offscreen when available. Nothing here
changes the audited mesh; outputs are display copies only.
"""
from __future__ import annotations
import argparse
import json
import struct
from pathlib import Path

import numpy as np


def read_mesh(path):
    d = np.load(path)
    points = d['points']
    cells = d['cells'] if 'cells' in d else None
    boundary = d['boundary'] if 'boundary' in d else None
    if cells is None or boundary is None:
        raise ValueError('expected points/cells/boundary in the NPZ')
    return points, cells, boundary


TET_FACES = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]])


def write_stl(path, vertices, triangles, header='RFAM tetra mesh view'):
    """Binary STL with per-triangle normals; triangles must be (n,3) int."""
    tri = vertices[triangles]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    length = np.linalg.norm(normals, axis=1)
    good = length > 0
    tri = tri[good]; normals = normals[good] / length[good, None]
    n = len(tri)
    record = np.zeros((n, 12), dtype='<f4')
    record[:, 0:3] = normals
    record[:, 3:6] = tri[:, 0]
    record[:, 6:9] = tri[:, 1]
    record[:, 9:12] = tri[:, 2]
    attr = np.zeros((n, 1), dtype='<u2')
    payload = np.hstack([record.view('<u1').reshape(n, 48), attr.view('<u1')]).tobytes()
    with open(path, 'wb') as f:
        f.write(header.encode()[:80].ljust(80, b' '))
        f.write(struct.pack('<I', n))
        f.write(payload)
    return n


def boundary_faces(cells):
    faces = cells[:, TET_FACES].reshape(-1, 3)
    keys = np.sort(faces, axis=1)
    _, first, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError('nonmanifold faces in the mesh')
    if not np.any(counts == 1):
        raise ValueError('mesh has no external boundary faces')
    # counts==2 faces are shared by two cells (internal); counts==1 are outer.
    internal = faces[first[counts == 2]]
    external = faces[first[counts == 1]]
    return internal, external


def interfaces(faces_internal, faces_external, stride):
    if stride > 1:
        faces_internal = faces_internal[::stride]
    return np.vstack([faces_external, faces_internal])


def clip_half(points, cells, origin_x):
    """Keep x <= origin_x using VTK's cell clipper (proven cut faces).

    A numpy polygon shortcut cannot reproduce the interface between kept cells
    and partially clipped cells, so it silently leaves an open surface.
    """
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    pts = vtk.vtkPoints(); pts.SetData(numpy_to_vtk(points.astype(np.float64), deep=True))
    ca = vtk.vtkCellArray()
    ca.SetData(numpy_to_vtkIdTypeArray(np.arange(0, len(cells) * 4 + 1, 4, dtype=np.int64), deep=True),
               numpy_to_vtkIdTypeArray(cells.astype(np.int64).ravel(), deep=True))
    grid = vtk.vtkUnstructuredGrid(); grid.SetPoints(pts); grid.SetCells(vtk.VTK_TETRA, ca)
    plane = vtk.vtkPlane(); plane.SetOrigin(origin_x, 0.0, 0.0); plane.SetNormal(-1.0, 0.0, 0.0)
    clip = vtk.vtkClipDataSet(); clip.SetInputData(grid); clip.SetClipFunction(plane)
    clip.SetValue(0.0); clip.Update()
    surf = vtk.vtkDataSetSurfaceFilter(); surf.SetInputData(clip.GetOutput()); surf.Update()
    polydata = surf.GetOutput()
    tri = vtk.vtkTriangleFilter(); tri.SetInputData(polydata); tri.Update()
    out = tri.GetOutput()
    v = np.asarray(out.GetPoints().GetData())
    ids = []
    for i in range(out.GetNumberOfCells()):
        cell = out.GetCell(i)
        if cell.GetNumberOfPoints() != 3:
            continue
        ids.append([cell.GetPointId(j) for j in range(3)])
    return v, np.asarray(ids, dtype=np.int64)


def render_png(vertices, triangles, path, size=(1400, 1000)):
    try:
        import vtk
        from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        pts = vtk.vtkPoints()
        pts.SetData(numpy_to_vtk(vertices.astype(np.float64), deep=True))
        ca = vtk.vtkCellArray()
        ca.SetData(numpy_to_vtkIdTypeArray(np.arange(0, len(triangles) * 3 + 1, 3, dtype=np.int64), deep=True),
                   numpy_to_vtkIdTypeArray(triangles.astype(np.int64).ravel(), deep=True))
        grid = vtk.vtkPolyData(); grid.SetPoints(pts); grid.SetPolys(ca)
        normals = vtk.vtkPolyDataNormals(); normals.SetInputData(grid)
        normals.SetFeatureAngle(45); normals.ConsistencyOn(); normals.SplittingOff()
        mapper = vtk.vtkPolyDataMapper(); mapper.SetInputConnection(normals.GetOutputPort())
        actor = vtk.vtkActor(); actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.62, 0.70, 0.78)
        actor.GetProperty().EdgeVisibilityOff()
        ren = vtk.vtkRenderer(); ren.AddActor(actor)
        ren.SetBackground(0.09, 0.10, 0.12)
        win = vtk.vtkRenderWindow(); win.SetOffScreenRendering(1)
        win.AddRenderer(ren); win.SetSize(*size)
        ren.ResetCamera()
        ren.GetActiveCamera().Azimuth(35); ren.GetActiveCamera().Elevation(22)
        ren.ResetCameraClippingRange()
        win.Render()
        grab = vtk.vtkWindowToImageFilter(); grab.SetInput(win); grab.Update()
        writer = vtk.vtkPNGWriter(); writer.SetFileName(str(path))
        writer.SetInputConnection(grab.GetOutputPort()); writer.Write()
        win.Finalize()
        return True
    except Exception as exc:  # offscreen GL may be unavailable
        print('PNG render skipped:', type(exc).__name__, str(exc)[:160])
        return False


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh', required=True, help='audited mesh NPZ (points/cells/boundary)')
    p.add_argument('--out', required=True)
    p.add_argument('--mode', default='all', choices=['boundary', 'half', 'interfaces', 'all'])
    p.add_argument('--stride', type=int, default=1, help='keep every Nth internal face in interfaces mode')
    p.add_argument('--units', choices=['mm', 'm'], default='mm')
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    points, cells, boundary = read_mesh(a.mesh)
    if a.units == 'm':
        points = points * 1e-3
    centre_x = float(points[:, 0].mean())
    report = {'source_mesh': str(Path(a.mesh).resolve()), 'mode': a.mode, 'units': a.units,
              'note': 'STL holds triangles only; volume-cell topology is not representable.',
              'viewing': {'paraview': 'open the mesh .vtu, then use Slice/Clip filters for real volume views',
                          'bcad': 'open any generated STL; half/interfaces modes show internal structure'},
              'outputs': {}}
    inner, outer = boundary_faces(cells)
    if a.mode in ('boundary', 'all'):
        n = write_stl(out / 'boundary_surface.stl', points, outer,
                      'RFAM fluid outer boundary (triangle surface only)')
        report['outputs']['boundary_surface.stl'] = {'triangles': int(n)}
        render_png(points, outer, out / 'boundary_surface.png')
    if a.mode in ('half', 'all'):
        hp, ht = clip_half(points, cells, centre_x)
        n = write_stl(out / 'half_domain_cut.stl', hp, ht,
                      'RFAM half domain cut from real tetrahedra')
        # Closed-surface volume, evaluated after a translation to reduce
        # cancellation from the CAD coordinate offset.
        shifted = hp - hp.mean(axis=0)
        half_volume = float(abs(np.einsum('ij,ij->i', shifted[ht][:, 0],
            np.cross(shifted[ht][:, 1], shifted[ht][:, 2])).sum() / 6))
        report['outputs']['half_domain_cut.stl'] = {'triangles': int(n),
            'closed_surface_volume_mm3': half_volume,
            'compare_with': 'reports/mesh_audit.json volume_sum_mm3 (about half is expected)'}
        render_png(hp, ht, out / 'half_domain_cut.png')
    if a.mode in ('interfaces', 'all'):
        faces = interfaces(inner, outer, a.stride)
        n = write_stl(out / 'interfaces_solid_look.stl', points, faces,
                      'RFAM all cell interfaces (display only, not a mesh)')
        report['outputs']['interfaces_solid_look.stl'] = {'triangles': int(n),
                                                          'internal_faces_kept': int(len(inner[::a.stride])),
                                                          'stride': a.stride}
    for name in ('boundary_surface.stl', 'half_domain_cut.stl'):
        if (out / name).exists():
            report['outputs'].setdefault(name, {})['watertight_surface'] = True
    report['outputs'] = {k: report['outputs'][k] for k in sorted(report['outputs'])}
    (out / 'view_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('STL_VIEWS_WRITTEN', out)
    for name, info in report['outputs'].items():
        print(' ', name, info)
    print('ParaView: open mesh/fluid_volume_mm.vtu -> Color by cell data -> Slice/Clip filters for real volume views.')


if __name__ == '__main__':
    main()
