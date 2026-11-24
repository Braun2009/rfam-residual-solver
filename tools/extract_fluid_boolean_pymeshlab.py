#!/usr/bin/env python3
"""Fully automated fluid domain extraction using PyMeshLab Boolean operations.

This bypasses ANSYS/Fluent entirely and creates the fluid domain STL
using pure Python libraries.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pymeshlab
import trimesh


AXIS_XY = np.array([49.21534822760748, 48.35004867381059])
BODY_SHA256 = '89a9579104f1ef532c7814fa6c612205e5a69c4ce1fa037163e5c6ba7a67f7d6'


def create_cylinder_mesh(center_xy, radius, height, segments=64):
    """Create a cylinder mesh using trimesh."""
    # Create cylinder along Z axis
    cyl = trimesh.creation.cylinder(radius=radius, height=height, sections=segments)
    # Translate to position: center XY, start at z=0
    cyl.apply_translation([center_xy[0], center_xy[1], height/2])
    return cyl


def extract_fluid_boolean(body_stl, output_dir, cylinder_radius=35, cylinder_height=280):
    """Extract fluid domain using PyMeshLab Boolean difference."""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print('Loading solid STL...')
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(str(body_stl))
    
    solid = ms.current_mesh()
    print(f'  Solid: {solid.vertex_number()} vertices, {solid.face_number()} faces')
    
    # Create cylinder
    print(f'Creating bounding cylinder (r={cylinder_radius}, h={cylinder_height})...')
    cylinder_mesh = create_cylinder_mesh(AXIS_XY, cylinder_radius, cylinder_height)
    
    # Save cylinder temporarily
    temp_cylinder = output_dir / 'temp_cylinder.stl'
    cylinder_mesh.export(temp_cylinder)
    
    # Load cylinder into MeshSet
    ms.load_new_mesh(str(temp_cylinder))
    cyl = ms.current_mesh()
    print(f'  Cylinder: {cyl.vertex_number()} vertices, {cyl.face_number()} faces')
    
    # PyMeshLab Boolean difference: cylinder - solid
    print('Performing Boolean difference (cylinder - solid)...')
    print('  This may take several minutes...')
    
    start = time.time()
    try:
        # Boolean difference: first mesh (cylinder) - second mesh (solid)
        ms.generate_boolean_difference(first_mesh=1, second_mesh=0)
        elapsed = time.time() - start
        print(f'  Boolean operation completed in {elapsed:.1f}s')
        
        # The result is added as a new mesh
        fluid = ms.current_mesh()
        print(f'  Result: {fluid.vertex_number()} vertices, {fluid.face_number()} faces')
        
        # Clean up
        print('Cleaning result...')
        ms.meshing_remove_duplicate_vertices()
        ms.meshing_remove_duplicate_faces()
        ms.meshing_remove_unreferenced_vertices()
        
        fluid_clean = ms.current_mesh()
        print(f'  After cleanup: {fluid_clean.vertex_number()} vertices, {fluid_clean.face_number()} faces')
        
        # Save
        output_stl = output_dir / 'fluid_domain_boolean.stl'
        print(f'Saving to {output_stl}...')
        ms.save_current_mesh(str(output_stl))
        
        # Verify with trimesh
        print('\nVerifying result...')
        result_mesh = trimesh.load_mesh(str(output_stl), process=False)
        print(f'  Volume: {result_mesh.volume:.1f} mm^3')
        print(f'  Watertight: {result_mesh.is_watertight}')
        
        # Check cross-sections
        sec0 = result_mesh.section(plane_origin=[0,0,1], plane_normal=[0,0,1])
        sec261 = result_mesh.section(plane_origin=[0,0,261], plane_normal=[0,0,1])
        
        if sec0:
            print(f'  z=1: {len(sec0.discrete)} contours')
        if sec261:
            holes = sum(1 for path in sec261.discrete
                       if 3 < abs(np.sum(path[:-1,0]*path[1:,1]-path[1:,0]*path[:-1,1])/2) < 15
                       and 18 < np.linalg.norm(path[:-1,:2].mean(axis=0)-AXIS_XY) < 28)
            print(f'  z=261: {holes} propane hole candidates')
        
        manifest = {
            'method': 'pymeshlab_boolean_difference',
            'cylinder_radius_mm': cylinder_radius,
            'cylinder_height_mm': cylinder_height,
            'result_vertices': fluid_clean.vertex_number(),
            'result_faces': fluid_clean.face_number(),
            'result_volume_mm3': float(result_mesh.volume),
            'watertight': result_mesh.is_watertight,
            'time_seconds': elapsed
        }
        
        (output_dir / 'boolean_extraction_manifest.json').write_text(
            json.dumps(manifest, indent=2) + '\n')
        
        # Cleanup temp file
        temp_cylinder.unlink()
        
        print('\n=== SUCCESS ===')
        print(f'Fluid domain STL: {output_stl}')
        print('Next: Tetrahedralize with tetgen or import into mesher')
        
        return output_stl, manifest
        
    except Exception as e:
        print(f'ERROR during Boolean operation: {e}')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--body', default='D:/Vortex/cad/injector/injector_v3.STL')
    parser.add_argument('--out', default='work/path3_boolean')
    parser.add_argument('--radius', type=float, default=35)
    parser.add_argument('--height', type=float, default=280)
    args = parser.parse_args()
    
    extract_fluid_boolean(args.body, args.out, args.radius, args.height)


if __name__ == '__main__':
    main()
