#!/usr/bin/env python3
"""
Visual geometry inspection of Vortex injector STL
Verify pintle dimensions match design document specifications
"""

import numpy as np
import struct
import json

def read_stl_binary(filename):
    """Read binary STL file"""
    with open(filename, 'rb') as f:
        header = f.read(80)
        n_triangles = struct.unpack('I', f.read(4))[0]
        
        print(f"Reading {filename}")
        print(f"Header: {header[:50].decode('ascii', errors='ignore')}")
        print(f"Triangles: {n_triangles:,}")
        print()
        
        # Read all triangles
        vertices = []
        for i in range(n_triangles):
            normal = struct.unpack('fff', f.read(12))
            v1 = struct.unpack('fff', f.read(12))
            v2 = struct.unpack('fff', f.read(12))
            v3 = struct.unpack('fff', f.read(12))
            attr = struct.unpack('H', f.read(2))
            
            vertices.extend([v1, v2, v3])
        
        vertices = np.array(vertices)
        print(f"Total vertices: {len(vertices):,}")
        return vertices

def analyze_geometry(vertices, name="Geometry"):
    """Analyze STL geometry statistics"""
    
    print("=" * 80)
    print(f"{name} - Geometry Analysis")
    print("=" * 80)
    print()
    
    # Bounding box
    x_min, y_min, z_min = vertices.min(axis=0)
    x_max, y_max, z_max = vertices.max(axis=0)
    
    print("Bounding Box:")
    print(f"  X: {x_min:.3f} to {x_max:.3f} mm (span: {x_max-x_min:.3f} mm)")
    print(f"  Y: {y_min:.3f} to {y_max:.3f} mm (span: {y_max-y_min:.3f} mm)")
    print(f"  Z: {z_min:.3f} to {z_max:.3f} mm (span: {z_max-z_min:.3f} mm)")
    print()
    
    # Center
    center_x = (x_min + x_max) / 2
    center_y = (y_min + y_max) / 2
    center_z = (z_min + z_max) / 2
    
    print(f"Geometric Center: ({center_x:.3f}, {center_y:.3f}, {center_z:.3f}) mm")
    print()
    
    return {
        "bbox": {
            "x_min": float(x_min), "x_max": float(x_max),
            "y_min": float(y_min), "y_max": float(y_max),
            "z_min": float(z_min), "z_max": float(z_max)
        },
        "center": {
            "x": float(center_x),
            "y": float(center_y),
            "z": float(center_z)
        },
        "span": {
            "x": float(x_max - x_min),
            "y": float(y_max - y_min),
            "z": float(z_max - z_min)
        }
    }

def find_pintle_features(vertices):
    """Identify pintle tip features by analyzing Z-plane cross-sections"""
    
    print("=" * 80)
    print("Pintle Feature Detection")
    print("=" * 80)
    print()
    
    # From injector_stl_analysis.md:
    # - Pintle tips are at z ~ 281 mm
    # - 3 pintles at 120° spacing
    # - Bolt circle radius R ~ 12.5 mm
    # - Center around (49, 51) mm
    
    z_tip = 281  # mm, expected pintle tip location
    z_tol = 5    # mm tolerance
    
    # Find vertices near pintle tip
    tip_vertices = vertices[(vertices[:, 2] > z_tip - z_tol) & 
                            (vertices[:, 2] < z_tip + z_tol)]
    
    print(f"Vertices near z={z_tip} mm (±{z_tol} mm): {len(tip_vertices):,}")
    
    if len(tip_vertices) == 0:
        print("WARNING: No vertices found at expected pintle tip location!")
        print(f"Expected z ~ {z_tip} mm")
        print(f"Actual z range: {vertices[:, 2].min():.1f} to {vertices[:, 2].max():.1f} mm")
        return None
    
    # Cluster analysis to find 3 pintle centers
    # Assuming pintles are at ~120° spacing around center
    
    center_x = 49.2  # From STL analysis doc
    center_y = 51.5
    
    # Convert to polar coordinates
    dx = tip_vertices[:, 0] - center_x
    dy = tip_vertices[:, 1] - center_y
    r = np.sqrt(dx**2 + dy**2)
    theta = np.arctan2(dy, dx) * 180 / np.pi
    
    # Expected bolt circle radius
    expected_r = 12.5  # mm
    
    # Find vertices on bolt circle
    on_circle = tip_vertices[(r > expected_r - 5) & (r < expected_r + 5)]
    
    print(f"Vertices on bolt circle (r ~ {expected_r} mm): {len(on_circle):,}")
    print()
    
    # Find 3 angular clusters (120° apart)
    angles = [-90, 30, 150]  # Expected pintle angles (approximate)
    
    pintles = []
    for i, target_angle in enumerate(angles):
        # Find vertices within ±30° of target
        dx_i = on_circle[:, 0] - center_x
        dy_i = on_circle[:, 1] - center_y
        theta_i = np.arctan2(dy_i, dx_i) * 180 / np.pi
        
        # Wrap angles
        angle_diff = (theta_i - target_angle + 180) % 360 - 180
        in_sector = on_circle[np.abs(angle_diff) < 30]
        
        if len(in_sector) > 0:
            pintle_center = in_sector.mean(axis=0)
            pintles.append({
                "id": i,
                "center": pintle_center.tolist(),
                "n_vertices": len(in_sector),
                "angle": float(np.arctan2(pintle_center[1]-center_y, 
                                         pintle_center[0]-center_x) * 180 / np.pi)
            })
            print(f"Pintle {i}: center=({pintle_center[0]:.2f}, {pintle_center[1]:.2f}, {pintle_center[2]:.2f}), "
                  f"angle={pintles[-1]['angle']:.1f}°, vertices={len(in_sector):,}")
    
    print()
    return pintles

def measure_pintle_dimensions(vertices, pintle_info):
    """Measure pintle slot dimensions"""
    
    print("=" * 80)
    print("Pintle Dimension Measurements")
    print("=" * 80)
    print()
    
    # From design document:
    # - Pintle diameter D_p = 9.0 mm
    # - LOX slot height h_LOX = 0.742 mm
    # - Fuel annulus gap t_C3H8 = 2.119 mm
    
    print("Target dimensions (from design document):")
    print(f"  Pintle diameter: 9.0 mm")
    print(f"  LOX slot height: 0.742 mm")
    print(f"  Fuel annulus gap: 2.119 mm")
    print()
    
    print("Measurement from STL:")
    print("  Note: STL is triangulated surface mesh, precise dimension")
    print("  measurement requires CAD software or detailed point cloud analysis.")
    print()
    
    if pintle_info is None or len(pintle_info) == 0:
        print("  Cannot measure - pintle features not detected")
        return None
    
    # Rough estimation from vertex spread
    measurements = []
    for pintle in pintle_info:
        # This is a very rough estimate - proper measurement needs CAD
        print(f"  Pintle {pintle['id']}: Detected at angle {pintle['angle']:.1f}°")
        measurements.append(pintle)
    
    print()
    print("RECOMMENDATION:")
    print("  For precise measurements:")
    print("  1. Open injector_v4.SLDPRT in SolidWorks")
    print("  2. Use 'Measure' tool on pintle cross-section")
    print("  3. Or export STEP file and use FreeCAD measurement tools")
    print()
    
    return measurements

def main():
    print("=" * 80)
    print("VORTEX INJECTOR - STL GEOMETRY INSPECTION")
    print("=" * 80)
    print()
    
    # Load injector STL
    stl_file = "D:/Vortex/cad/injector/injector_v4_clean.STL"
    
    try:
        vertices = read_stl_binary(stl_file)
    except FileNotFoundError:
        print(f"ERROR: File not found: {stl_file}")
        print()
        print("Available STL files:")
        import glob
        files = glob.glob("D:/Vortex/cad/injector/*.STL") + glob.glob("D:/Vortex/cad/injector/*.stl")
        for f in files:
            print(f"  {f}")
        return
    
    # Analyze overall geometry
    geo_stats = analyze_geometry(vertices, "Injector V4")
    
    # Find pintle features
    pintles = find_pintle_features(vertices)
    
    # Measure dimensions
    measurements = measure_pintle_dimensions(vertices, pintles)
    
    # Export results
    results = {
        "file": stl_file,
        "geometry": geo_stats,
        "pintles_detected": len(pintles) if pintles else 0,
        "pintles": pintles,
        "measurements": measurements,
        "design_reference": {
            "pintle_diameter_mm": 9.0,
            "lox_slot_height_mm": 0.742,
            "fuel_annulus_gap_mm": 2.119,
            "lox_throat_area_mm2": 20.99,
            "fuel_throat_area_mm2": 59.92
        }
    }
    
    output_file = "D:/realFluidAllMach/work/stl_geometry_inspection.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("=" * 80)
    print("Results saved to:", output_file)
    print("=" * 80)
    print()
    
    print("CONCLUSION:")
    print("  STL geometry loaded successfully")
    print(f"  Detected {len(pintles) if pintles else 0} pintle features")
    print("  For precise dimension verification, use CAD software")
    print()
    print("NEXT STEP: Set up cold flow CFD with current mesh")
    print()

if __name__ == "__main__":
    main()
