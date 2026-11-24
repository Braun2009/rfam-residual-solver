#!/usr/bin/env python3
"""Reconstruct a nozzle-only fluid boundary from the authoritative STL solids.

All construction uses CAD millimetres, with original CAD Z downstream.
This produces geometry and diagnostic tetrahedra, not an RFAM restart or
an assertion that RFAM can consume an unstructured mesh.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import trimesh
import manifold3d as md


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def load_welded(path):
    raw = trimesh.load_mesh(str(path), process=False)
    vertices, inv = np.unique(raw.vertices, axis=0, return_inverse=True)
    m = trimesh.Trimesh(vertices, inv[raw.faces], process=False)
    if not m.is_watertight or not m.is_winding_consistent:
        raise ValueError(f'{path}: exact-coordinate welding did not yield closed oriented solids')
    return m


def as_solid(mesh):
    s = md.Manifold(md.Mesh64(np.asarray(mesh.vertices, dtype=np.float64),
                             np.asarray(mesh.faces, dtype=np.uint64)))
    if s.status() != md.Error.NoError or s.is_empty():
        raise ValueError(f'Invalid solid: {s.status()}')
    return s


def as_mesh(solid):
    if solid.status() != md.Error.NoError:
        raise ValueError(f'Boolean failure: {solid.status()}')
    m = solid.to_mesh64()
    return trimesh.Trimesh(np.asarray(m.vert_properties)[:, :3],
                          np.asarray(m.tri_verts), process=False)


def horizontal_planes(mesh, axis, select):
    use = select & (np.abs(mesh.face_normals[:,axis]) > 0.99999)
    planes = {}
    for v, a in zip(np.round(mesh.triangles_center[use,axis], 6), mesh.area_faces[use]):
        planes[float(v)] = planes.get(float(v), 0.0) + float(a)
    return sorted((x, a) for x, a in planes.items() if a > 1.0)


def circle_fit(points):
    # Algebraic fit used only on constant-radius bore sections.
    points = np.asarray(points)
    x, y = points.T
    c = np.linalg.lstsq(np.column_stack([2*x,2*y,np.ones(len(x))]), x*x+y*y, rcond=None)[0]
    radius = float(np.sqrt(c[2] + c[0]**2 + c[1]**2))
    err = np.linalg.norm(points-c[:2], axis=1) - radius
    return c[:2], radius, float(np.sqrt(np.mean(err**2)))


def locate_axes(body):
    section = body.section(plane_origin=[0,0,255.123], plane_normal=[0,0,1])
    axes = []
    for pts in section.discrete:
        c, r, rms = circle_fit(pts[:,:2])
        if 4.4 < r < 4.6 and rms < 0.03:
            axes.append({'xy_mm':c.tolist(), 'bore_radius_mm':r, 'fit_rms_mm':rms})
    if len(axes) != 3:
        raise ValueError(f'Expected three circular main bores, got {len(axes)}')
    centre = np.mean([x['xy_mm'] for x in axes], axis=0)
    axes.sort(key=lambda x: np.arctan2(x['xy_mm'][1]-centre[1], x['xy_mm'][0]-centre[0]))
    return centre, axes


def envelope(profile_path, centre, board_z, inlet_z, end_x=None):
    profile = json.loads(Path(profile_path).read_text(encoding='utf-8'))['profile']
    p = np.array([[q['x_mm'], q['r_mm']] for q in profile])
    if end_x is not None:
        r = np.interp(end_x, p[:,0], p[:,1])
        p = np.vstack([p[p[:,0]<end_x], [end_x,r]])
    # CAD envelope, revolved about Z. Segmented cones retain the JSON profile.
    zr = np.vstack([[inlet_z,p[0,1]], np.column_stack([board_z+p[:,0], p[:,1]])])
    contour = [[0,inlet_z]] + [[float(r),float(z)] for z,r in zr] + [[0,float(zr[-1,0])]]
    s = md.CrossSection([contour]).revolve(circular_segments=180).translate([*centre,0])
    return s


def export_geometry(args):
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    body = load_welded(args.body)
    tops = sorted(load_welded(args.top).split(only_watertight=False), key=lambda m:m.bounds[0,0])
    if len(tops) != 3:
        raise ValueError('Expected three separate top pieces')
    centre, axes = locate_axes(body)
    main = as_solid(body)
    manifest = {
        'status':'GEOMETRY_CANDIDATE_NOT_CFD_ACCEPTED',
        'units':'mm', 'coordinate_system':'original CAD, Z downstream',
        'sources': [{'path':str(Path(p).resolve()),'sha256':sha256(p)} for p in [args.body,args.top,args.profile]],
        'inlet_cad_z_mm':args.inlet_z, 'board_cad_z_mm':args.board_z,
        'centre_xy_mm':centre.tolist(), 'pintles':axes,
        'body_closed_components':1, 'top_closed_components':len(tops),
        'assembly_basis':'mate both 6 mm thread end planes; no spherical approximation',
        'limitations':['Thread rotation not recovered from assembly mates; buried thread cavity excluded only by connectivity.',
                       'Tetrahedra are for geometry review; RFAM structured solver does not read this format.',
                       'Inlet species labels and production near-wall/local sizing require separate acceptance.'],
    }
    start=time.time()
    pieces=[main]
    for top, a in zip(tops, axes):
        xy=np.array(a['xy_mm']); tc=(top.bounds[0]+top.bounds[1])/2
        centres=body.triangles_center
        near=(np.linalg.norm(centres[:,:2]-xy,axis=1)<1.9)&(centres[:,2]>270)&(centres[:,2]<279)
        bp=horizontal_planes(body,2,near)
        tp=horizontal_planes(top,1,np.ones(len(top.faces),bool))
        if len(bp)!=2 or len(tp)!=2:
            raise ValueError(f'Could not identify two mating planes: {bp} {tp}')
        offsets=np.array([bp[0][0]+tp[1][0],bp[1][0]+tp[0][0]])
        if np.ptp(offsets)>0.005:
            raise ValueError(f'Mating lengths disagree: {offsets}')
        tip=float(np.mean(offsets))
        # proper rotation, determinant +1: (X,Y,Z) -> (X,Z,-Y).
        transform=np.array([[1,0,0,xy[0]-tc[0]], [0,0,1,xy[1]-tc[2]], [0,-1,0,tip]])
        placed=as_solid(top).transform(transform)
        overlap=(main ^ placed).volume()
        a.update({'body_mating_planes_mm':bp,'top_mating_planes_mm':tp,
                  'top_transform_3x4':transform.tolist(),'tip_cad_z_mm':tip,
                  'body_top_overlap_mm3':overlap})
        pieces.append(placed)
        print('top',a,flush=True)
    assembled=md.Manifold.batch_boolean(pieces,md.OpType.Add)
    env=envelope(args.profile,centre,args.board_z,args.inlet_z,args.end_x)
    fluid=env-assembled
    if fluid.status()!=md.Error.NoError:raise ValueError(fluid.status())
    components=fluid.decompose()
    components.sort(key=lambda c:abs(c.volume()),reverse=True)
    manifest['fluid_boolean_components']=[{'volume_mm3':c.volume(),'bounds':list(c.bounding_box()),'triangles':c.num_tri()} for c in components]
    # Select by contact with the downstream outlet, not component size alone.
    outlet_z=env.bounding_box()[5]
    connected=[c for c in components if abs(c.bounding_box()[5]-outlet_z)<1e-6]
    if len(connected)!=1:
        raise ValueError(f'Expected one outlet-connected fluid solid, got {len(connected)}')
    fluid=connected[0]
    discarded=[c for c in components if c is not fluid]
    # The only excluded cavities may be inside the threaded seats.
    for c in discarded:
        bounds=c.bounding_box()
        if bounds[2]<271.85 or bounds[5]>277.86 or abs(c.volume())>2.0:
            raise ValueError('Disconnected fluid outside a buried thread seat needs review')
    # Evaluate just inside the boundary. Exact coplanar slicing is ambiguous
    # after Boolean operations on rounded STL coordinates.
    apertures=fluid.slice(args.inlet_z+1e-4).decompose()
    areas=sorted([c.area() for c in apertures])
    if len(areas)!=64 or sum(3.0<a<3.3 for a in areas)!=60:
        raise ValueError(f'Nozzle-only inlet must have 60 small holes + 3 main passages + torch; got {len(areas)}')
    manifest['inlet_apertures']=len(areas)
    manifest['inlet_aperture_areas_mm2']=areas
    fm=as_mesh(fluid)
    if not fm.is_watertight:raise ValueError('Fluid boolean not watertight')
    manifest.update({'retained_fluid_volume_mm3':fluid.volume(),'retained_surface_triangles':len(fm.faces),
                     'boundary_closed':fm.is_watertight,'boundary_euler':fm.euler_number,
                     'discarded_volume_mm3':sum(c.volume() for c in discarded),
                     'elapsed_s':time.time()-start})
    # NPZ is the mesher's float64 source. PLY/STL viewers generally quantize
    # positions to float32 and must not be used as a round-trip meshing input.
    np.savez_compressed(out/'fluid_boundary_mm.npz', vertices=fm.vertices, faces=fm.faces)
    fm.export(out/'fluid_boundary_mm.ply')
    fm.export(out/'fluid_boundary_mm.stl')
    # Cropping uses a closed solid operation; surface triangles are never simply deleted.
    head=assembled.trim_by_plane([0,0,1],args.inlet_z)
    hm=as_mesh(head)
    np.savez_compressed(out/'assembled_head_mm.npz', vertices=hm.vertices, faces=hm.faces)
    hm.export(out/'assembled_head_mm.ply')
    # One real half-domain for BCAD inspection: the cut is capped by solid geometry.
    cut=fluid.trim_by_plane([1,0,0],float(centre[0]))
    as_mesh(cut).export(out/'fluid_half_domain_mm.stl')
    (out/'geometry_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('GEOMETRY_EXPORTED',len(fm.faces),fluid.volume(), 'components',len(components),flush=True)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--body',required=True);p.add_argument('--top',required=True)
    p.add_argument('--profile',required=True);p.add_argument('--out',required=True)
    p.add_argument('--board-z',type=float,default=270.9527893066406)
    p.add_argument('--inlet-z',type=float,default=269.2)
    p.add_argument('--end-x',type=float,default=None,help='Optional downstream cut for a local diagnostic')
    args=p.parse_args()
    export_geometry(args)

if __name__=='__main__':main()
