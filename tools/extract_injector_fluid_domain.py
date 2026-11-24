#!/usr/bin/env python3
"""Reconstruct the complete injector fluid domain from original STL.

Extracts internal flow passages from the solid body, caps inlet openings,
and generates a CFD-quality tetrahedral mesh with proper boundary labels.
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


AXIS_XY = np.array([49.21534822760748, 48.35004867381059])
BODY_SHA256 = '89a9579104f1ef532c7814fa6c612205e5a69c4ce1fa037163e5c6ba7a67f7d6'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def as_manifold(mesh):
    """Convert trimesh to manifold3d."""
    return md.Manifold(mesh=md.Mesh(
        vert_properties=np.asarray(mesh.vertices, dtype=np.float32),
        tri_verts=np.asarray(mesh.faces, dtype=np.uint32)))


def as_trimesh(manifold):
    """Convert manifold3d to trimesh."""
    m = manifold.to_mesh()
    return trimesh.Trimesh(vertices=np.asarray(m.vert_properties),
                          faces=np.asarray(m.tri_verts), process=False)


def create_lox_inlet_cap(inner_radius=4.55, outer_radius=6.35, z=0):
    """Create a planar cap for the LOX inlet hole at z=0."""
    theta = np.linspace(0, 2*np.pi, 39, endpoint=True)
    inner_loop = np.column_stack([AXIS_XY[0] + inner_radius * np.cos(theta),
                                   AXIS_XY[1] + inner_radius * np.sin(theta),
                                   np.full(len(theta), z)])
    # Create a simple disk cap using triangulation
    center = np.array([[AXIS_XY[0], AXIS_XY[1], z]])
    faces = []
    for i in range(len(inner_loop) - 1):
        faces.append([0, i+1, i+2])
    verts = np.vstack([center, inner_loop])
    cap = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    # Ensure correct normal direction (pointing into fluid domain, -Z)
    if cap.face_normals[0, 2] > 0:
        cap.faces = cap.faces[:, ::-1]
    return cap


def create_propane_hole_cap(angle_deg, radius=23.2, hole_radius=1.4, z=261.0):
    """Create a planar cap for one propane hole."""
    theta_center = np.radians(angle_deg)
    center_xy = AXIS_XY + radius * np.array([np.cos(theta_center), np.sin(theta_center)])
    
    # Simple circular cap
    theta = np.linspace(0, 2*np.pi, 25, endpoint=True)
    loop = np.column_stack([center_xy[0] + hole_radius * np.cos(theta),
                           center_xy[1] + hole_radius * np.sin(theta),
                           np.full(len(theta), z)])
    center = np.array([[center_xy[0], center_xy[1], z]])
    faces = []
    for i in range(len(loop) - 1):
        faces.append([0, i+1, i+2])
    verts = np.vstack([center, loop])
    cap = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    # Normal should point upstream (-Z for inlet at higher z)
    if cap.face_normals[0, 2] > 0:
        cap.faces = cap.faces[:, ::-1]
    return cap


def extract_fluid_domain(body_path, chamber_profile_path, output_dir):
    """Extract fluid domain from solid body and cap inlets."""
    if sha256(body_path) != BODY_SHA256:
        raise ValueError('injector_v3.STL hash mismatch')
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load and validate solid
    raw = trimesh.load_mesh(str(body_path), process=False)
    verts, inv = np.unique(raw.vertices, axis=0, return_inverse=True)
    solid_mesh = trimesh.Trimesh(vertices=verts, faces=inv[raw.faces], process=False)
    
    if not solid_mesh.is_watertight or not solid_mesh.is_winding_consistent:
        raise ValueError('solid body is not closed/oriented')
    
    print('Solid body: {} faces, volume {} mm3'.format(len(solid_mesh.faces), solid_mesh.volume))
    
    # Load chamber profile for downstream envelope
    profile = json.loads(Path(chamber_profile_path).read_text(encoding='utf-8'))
    profile_points = profile['profile']
    
    # Create bounding envelope: upstream cylinder + downstream chamber contour
    # Upstream: cylinder from z=0 to z=269, radius sufficient to contain body
    upstream_r = 50.0  # mm, contains all flow passages
    downstream_z_end = 290.0  # mm, past injector tip
    
    # Build envelope as revolution around Z axis
    # For now, use a simple approach: create outer cylinder + subtract solid
    
    # Step 1: Create inlet caps
    lox_cap = create_lox_inlet_cap()
    
    propane_caps = []
    for k in range(48):
        angle = -90 + 7.5 * k  # Start from user's selected hole
        propane_caps.append(create_propane_hole_cap(angle))
    
    # Step 2: Cap the solid body at inlet openings
    # We need to "plug" the holes so Boolean subtraction creates clean inlets
    # This is complex - the solid has open holes, we need to identify and cap them
    
    # Simplified approach for now: work with the solid as-is
    # and identify boundary faces after Boolean operation
    
    # Export intermediate geometry
    solid_mesh.export(output_dir / 'injector_solid.stl')
    lox_cap.export(output_dir / 'lox_inlet_cap.stl')
    
    combined_caps = trimesh.util.concatenate([lox_cap] + propane_caps)
    combined_caps.export(output_dir / 'propane_inlet_caps.stl')
    
    manifest = {
        'status': 'GEOMETRY_EXTRACTION_IN_PROGRESS',
        'body_sha256': BODY_SHA256,
        'lox_inlet_cap_area_mm2': float(lox_cap.area),
        'propane_hole_count': 48,
        'propane_cap_area_mm2_total': float(sum(c.area for c in propane_caps)),
        'note': 'Full fluid domain extraction requires Boolean operations on closed manifolds; not yet implemented'
    }
    
    (output_dir / 'fluid_extraction_manifest.json').write_text(
        json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    
    print('Inlet caps exported to:', output_dir)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--body', required=True, help='injector_v3.STL path')
    parser.add_argument('--chamber-profile', required=True, help='chamber_profile_3pintle.json')
    parser.add_argument('--out', required=True, help='output directory')
    args = parser.parse_args()
    
    extract_fluid_domain(args.body, args.chamber_profile, args.out)


if __name__ == '__main__':
    main()
