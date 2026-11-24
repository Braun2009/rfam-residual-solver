#!/usr/bin/env python3
"""
Analyze RFAM mesh to understand zone configuration and compare to design document
"""

import struct
import numpy as np

def read_rfamtet_header(filename):
    """Read RFAMTET mesh header to extract zone information"""
    
    print("=" * 80)
    print(f"ANALYZING MESH: {filename}")
    print("=" * 80)
    print()
    
    with open(filename, 'rb') as f:
        # Read magic number
        magic = f.read(10).decode('ascii', errors='ignore')
        print(f"Magic: {magic}")
        
        # Read counts (assuming binary format)
        f.seek(80)  # Skip header
        
        try:
            n_points = struct.unpack('I', f.read(4))[0]
            n_cells = struct.unpack('I', f.read(4))[0]
            n_faces = struct.unpack('I', f.read(4))[0]
            n_zones = struct.unpack('I', f.read(4))[0]
            
            print(f"Mesh Statistics:")
            print(f"  Points: {n_points:,}")
            print(f"  Cells:  {n_cells:,}")
            print(f"  Faces:  {n_faces:,}")
            print(f"  Zones:  {n_zones}")
            print()
            
            # Try to read zone info
            print(f"Zone Information:")
            print("-" * 80)
            
            # This is a simplified parser - actual format may differ
            # We need to check the RFAM documentation for exact binary format
            
        except Exception as e:
            print(f"Error parsing binary format: {e}")
            print()
            print("Note: RFAMTET is a binary format. Need to:")
            print("  1. Check RFAM documentation for exact binary structure")
            print("  2. Or run RFAM tool to query mesh info")
            print("  3. Or check if there's a companion .json or .txt file")

def check_companion_files(mesh_path):
    """Check for companion metadata files"""
    import os
    
    base_dir = os.path.dirname(mesh_path)
    mesh_name = os.path.basename(mesh_path)
    
    print()
    print("Checking for companion files...")
    print("-" * 80)
    
    # Look for common metadata files
    patterns = [
        "*.json",
        "*.txt", 
        "*audit*",
        "*zone*",
        "*boundary*"
    ]
    
    import glob
    
    for pattern in patterns:
        files = glob.glob(os.path.join(base_dir, pattern))
        if files:
            print(f"\nFound {pattern}:")
            for f in files[:5]:
                print(f"  {os.path.basename(f)}")

if __name__ == "__main__":
    mesh_file = "D:/realFluidAllMach/work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet"
    
    read_rfamtet_header(mesh_file)
    check_companion_files(mesh_file)
    
    print()
    print("=" * 80)
    print("RECOMMENDATION:")
    print("=" * 80)
    print()
    print("To properly analyze the mesh zones, we should:")
    print("  1. Check if there's mesh metadata in JSON format")
    print("  2. Or run a diagnostic RFAM command to query zones")
    print("  3. Or check the CFD setup scripts for zone definitions")
    print()
