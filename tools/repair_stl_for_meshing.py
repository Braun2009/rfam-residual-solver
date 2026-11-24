#!/usr/bin/env python3
"""Repair and clean the injector STL for robust meshing.

Uses PyMeshLab to fix common STL issues that cause mesher failures:
- Non-manifold edges
- Self-intersections
- Degenerate faces
- Inconsistent normals
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

try:
    import pymeshlab
except ImportError:
    print("ERROR: pymeshlab not installed")
    print("Install with: pip install pymeshlab")
    raise


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def repair_stl(input_path, output_path, aggressive=False):
    """Repair STL mesh using MeshLab filters."""
    
    print('Loading mesh:', input_path)
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(str(input_path))
    
    m = ms.current_mesh()
    print('  Vertices:', m.vertex_number())
    print('  Faces:', m.face_number())
    
    # Check topology if available
    try:
        topo_init = ms.get_topological_measures()
        print('  Initial topology:', topo_init.get('number_components', 'N/A'), 'components')
    except:
        print('  Topology check skipped (API unavailable)')
    
    print('\nApplying repair filters...')
    
    # 1. Remove duplicate faces
    print('  Removing duplicate faces...')
    ms.meshing_remove_duplicate_faces()
    
    # 2. Remove duplicate vertices
    print('  Removing duplicate vertices...')
    ms.meshing_remove_duplicate_vertices()
    
    # 3. Remove unreferenced vertices
    print('  Removing unreferenced vertices...')
    ms.meshing_remove_unreferenced_vertices()
    
    # 4. Remove null faces (zero area)
    print('  Removing null faces...')
    ms.meshing_remove_null_faces()
    
    # 5. Fix non-manifold edges (split or collapse)
    print('  Fixing non-manifold edges...')
    ms.meshing_repair_non_manifold_edges(method='Remove Faces')
    
    # 6. Fix non-manifold vertices
    print('  Fixing non-manifold vertices...')
    ms.meshing_repair_non_manifold_vertices(vertdispratio=0.0)
    
    # 7. Unify normals
    print('  Unifying normals...')
    ms.meshing_re_orient_faces_coherently()
    
    if aggressive:
        print('\nAggressive repair mode:')
        
        # 8. Close holes
        print('  Closing small holes...')
        ms.meshing_close_holes(maxholesize=30)
        
        # 9. Smooth to remove self-intersections
        print('  Smoothing to resolve self-intersections...')
        ms.apply_coord_laplacian_smoothing(stepsmoothnum=3, cotangentweight=False)
        
        # 10. Remesh if needed
        print('  Remeshing for uniformity...')
        ms.meshing_isotropic_explicit_remeshing(iterations=3, adaptive=True)
    
    # Final cleanup
    print('\nFinal cleanup...')
    ms.meshing_remove_duplicate_faces()
    ms.meshing_remove_unreferenced_vertices()
    
    # Check result
    m_fixed = ms.current_mesh()
    
    print('\nRepaired mesh:')
    print('  Vertices:', m_fixed.vertex_number())
    print('  Faces:', m_fixed.face_number())
    
    # Try to get topology measures
    try:
        topo = ms.get_topological_measures()
        is_closed = topo.get('is_mesh_closed', topo.get('number_holes', -1) == 0)
        print('  Watertight:', is_closed)
        print('  Connected components:', topo.get('number_components', 'N/A'))
    except Exception as e:
        print('  Topology measures unavailable:', str(e)[:50])
        topo = {'is_mesh_closed': None, 'number_components': 1}
    
    # Save
    print('\nSaving to:', output_path)
    ms.save_current_mesh(str(output_path))
    
    return {
        'input_vertices': m.vertex_number(),
        'input_faces': m.face_number(),
        'output_vertices': m_fixed.vertex_number(),
        'output_faces': m_fixed.face_number(),
        'is_watertight': topo.get('is_mesh_closed', None),
        'connected_components': topo.get('number_components', 1)
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='input STL')
    parser.add_argument('--output', required=True, help='output repaired STL')
    parser.add_argument('--aggressive', action='store_true',
                       help='apply aggressive repairs (hole filling, remeshing)')
    parser.add_argument('--verify-hash', 
                       default='89a9579104f1ef532c7814fa6c612205e5a69c4ce1fa037163e5c6ba7a67f7d6',
                       help='expected SHA256 of input')
    args = parser.parse_args()
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    
    # Verify input hash
    input_hash = sha256(input_path)
    if args.verify_hash and input_hash != args.verify_hash:
        print(f'WARNING: Input hash mismatch')
        print(f'  Expected: {args.verify_hash}')
        print(f'  Got:      {input_hash}')
        response = input('Continue anyway? [y/N] ')
        if response.lower() != 'y':
            return
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    result = repair_stl(input_path, output_path, args.aggressive)
    
    print('\n=== Repair Summary ===')
    print(f'Vertices: {result["input_vertices"]} → {result["output_vertices"]}')
    print(f'Faces: {result["input_faces"]} → {result["output_faces"]}')
    print(f'Watertight: {result["is_watertight"]}')
    print(f'Components: {result["connected_components"]}')
    
    if not result['is_watertight']:
        print('\nWARNING: Repaired mesh is still not watertight')
        print('Try --aggressive mode or manual repair in MeshLab GUI')
    else:
        print('\nSUCCESS: Mesh is now watertight and ready for meshing')


if __name__ == '__main__':
    main()
