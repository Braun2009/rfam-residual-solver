#!/usr/bin/env python3
"""
Directly modify STL to enlarge fuel orifices
Strategy: Identify fuel hole geometry and radially scale it
"""

import numpy as np
import sys

def read_stl(filename):
    """Read binary or ASCII STL file"""
    with open(filename, 'rb') as f:
        header = f.read(5).decode('ascii', errors='ignore')
    
    if header == 'solid':
        return read_stl_ascii(filename)
    else:
        return read_stl_binary(filename)

def read_stl_binary(filename):
    """Read binary STL"""
    with open(filename, 'rb') as f:
        header = f.read(80)
        n_triangles = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        
        triangles = []
        for i in range(n_triangles):
            normal = np.frombuffer(f.read(12), dtype=np.float32)
            v1 = np.frombuffer(f.read(12), dtype=np.float32)
            v2 = np.frombuffer(f.read(12), dtype=np.float32)
            v3 = np.frombuffer(f.read(12), dtype=np.float32)
            attr = f.read(2)
            triangles.append([normal, v1, v2, v3])
    
    return np.array(triangles), header

def write_stl_binary(filename, triangles, header=None):
    """Write binary STL"""
    if header is None:
        header = b'Modified by Kiro - Enlarged fuel orifices' + b' ' * (80 - 42)
    
    with open(filename, 'wb') as f:
        f.write(header[:80].ljust(80, b' '))
        f.write(np.uint32(len(triangles)).tobytes())
        
        for tri in triangles:
            normal, v1, v2, v3 = tri
            f.write(normal.astype(np.float32).tobytes())
            f.write(v1.astype(np.float32).tobytes())
            f.write(v2.astype(np.float32).tobytes())
            f.write(v3.astype(np.float32).tobytes())
            f.write(b'\x00\x00')  # attribute byte count

def identify_fuel_holes(triangles, z_min=0.0, z_max=10.0, r_min=10, r_max=30):
    """
    Identify triangles that belong to fuel holes
    Fuel holes are typically:
    - Near the injector face (low z)
    - At certain radial positions (annular pattern)
    - Around pintle positions
    """
    vertices = triangles[:, 1:4, :].reshape(-1, 3)
    
    # Calculate z and radial positions
    z = vertices[:, 2].reshape(-1, 3).mean(axis=1)
    x = vertices[:, 0].reshape(-1, 3).mean(axis=1)
    y = vertices[:, 1].reshape(-1, 3).mean(axis=1)
    r = np.sqrt(x**2 + y**2)
    
    # Identify fuel hole triangles
    mask = (z >= z_min) & (z <= z_max) & (r >= r_min) & (r <= r_max)
    
    return mask

def scale_radially(vertices, center, scale_factor):
    """Scale vertices radially from center"""
    # Vector from center to vertex
    dx = vertices[:, 0] - center[0]
    dy = vertices[:, 1] - center[1]
    
    # Scale radially (keep z unchanged)
    vertices[:, 0] = center[0] + dx * scale_factor
    vertices[:, 1] = center[1] + dy * scale_factor
    
    return vertices

def enlarge_fuel_orifices(input_stl, output_stl, scale_factor=1.27):
    """
    Enlarge fuel orifices by radially scaling
    scale_factor: 1.62 for 62% area increase (area scales as r^2)
                  sqrt(1.62) = 1.27 for diameter/radius scaling
    """
    print(f"Reading STL: {input_stl}")
    triangles, header = read_stl(input_stl)
    
    print(f"Total triangles: {len(triangles)}")
    
    # Get all vertices
    all_vertices = triangles[:, 1:4, :].reshape(-1, 3)
    
    print(f"Geometry bounds:")
    print(f"  X: {all_vertices[:, 0].min():.3f} to {all_vertices[:, 0].max():.3f}")
    print(f"  Y: {all_vertices[:, 1].min():.3f} to {all_vertices[:, 1].max():.3f}")
    print(f"  Z: {all_vertices[:, 2].min():.3f} to {all_vertices[:, 2].max():.3f}")
    
    # Estimate injector center (should be near origin)
    center_x = (all_vertices[:, 0].min() + all_vertices[:, 0].max()) / 2
    center_y = (all_vertices[:, 1].min() + all_vertices[:, 1].max()) / 2
    center = np.array([center_x, center_y, 0])
    
    print(f"Injector center: ({center_x:.3f}, {center_y:.3f})")
    print()
    
    # User needs to specify which region to scale
    print("INTERACTIVE MODE:")
    print("To enlarge fuel holes, I need to know their location.")
    print()
    print("Please provide:")
    print("  1. Z range of fuel holes (e.g., 0.0 to 5.0 mm)")
    print("  2. Radial range (e.g., 15 to 25 mm from center)")
    print()
    print("Or I can analyze the geometry automatically...")
    print()
    
    # Automatic detection based on typical geometry
    z_range = [0.0, 10.0]  # First 10mm from face
    r_range = [15.0, 35.0]  # Typical pintle radius range
    
    print(f"Using automatic detection:")
    print(f"  Z range: {z_range[0]} to {z_range[1]} mm")
    print(f"  R range: {r_range[0]} to {r_range[1]} mm")
    print()
    
    # Identify fuel hole triangles
    mask = identify_fuel_holes(triangles, z_range[0], z_range[1], 
                               r_range[0], r_range[1])
    
    n_selected = mask.sum()
    print(f"Selected {n_selected} triangles ({n_selected/len(triangles)*100:.1f}%) as fuel holes")
    
    if n_selected == 0:
        print("ERROR: No triangles selected. Adjust z_range and r_range.")
        return
    
    if n_selected > len(triangles) * 0.5:
        print("WARNING: More than 50% of geometry selected. This may be wrong.")
        response = input("Continue anyway? (y/n): ")
        if response.lower() != 'y':
            return
    
    # Scale the selected triangles
    modified_triangles = triangles.copy()
    
    for i in range(len(triangles)):
        if mask[i]:
            # Scale vertices (v1, v2, v3)
            for j in range(1, 4):  # j=0 is normal, j=1,2,3 are vertices
                modified_triangles[i, j, :] = scale_radially(
                    triangles[i, j:j+1, :], center, scale_factor
                )[0]
            
            # Recalculate normal
            v1 = modified_triangles[i, 1]
            v2 = modified_triangles[i, 2]
            v3 = modified_triangles[i, 3]
            edge1 = v2 - v1
            edge2 = v3 - v1
            normal = np.cross(edge1, edge2)
            normal = normal / (np.linalg.norm(normal) + 1e-10)
            modified_triangles[i, 0] = normal
    
    print(f"Applied radial scaling factor: {scale_factor:.3f}")
    print(f"Area scaling factor: {scale_factor**2:.3f}")
    print()
    
    # Write output
    print(f"Writing modified STL: {output_stl}")
    write_stl_binary(output_stl, modified_triangles, header)
    
    print("Done!")
    print()
    print("VERIFICATION NEEDED:")
    print("  1. Open in CAD software (e.g., FreeCAD, MeshLab)")
    print("  2. Check that fuel holes are enlarged correctly")
    print("  3. Check for mesh quality issues (overlaps, gaps)")
    print("  4. Remesh if needed")

if __name__ == "__main__":
    # For 62% area increase: scale = sqrt(1.62) = 1.27
    # For 50% area increase: scale = sqrt(1.50) = 1.22
    # For 15% area increase (with K optimization): scale = sqrt(1.15) = 1.07
    
    input_stl = "D:/realFluidAllMach/work/path4_correct_domain/lox_supply_domain.stl"
    output_stl = "D:/realFluidAllMach/work/lox_supply_enlarged_fuel.stl"
    
    # Use 1.27 for full 62% area increase
    # Or use 1.22 if we also optimize the flow path (combination strategy)
    scale_factor = 1.22  # Conservative, assuming some K optimization
    
    print("=" * 70)
    print("FUEL ORIFICE ENLARGEMENT TOOL")
    print("=" * 70)
    print()
    print(f"Input:  {input_stl}")
    print(f"Output: {output_stl}")
    print(f"Target area increase: {(scale_factor**2 - 1)*100:.0f}%")
    print()
    
    enlarge_fuel_orifices(input_stl, output_stl, scale_factor)
