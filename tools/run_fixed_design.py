#!/usr/bin/env python3
"""
Run RFAM cold flow with FIXED design parameters
Uses Option A: Increased fuel supply pressure to 5.26 MPa
"""

import json
import subprocess
import sys

# Load the fixed design
with open("D:/realFluidAllMach/work/injector_design_fixed.json", 'r') as f:
    data = json.load(f)

fixed_params = data['fixed_designs']['A']['parameters']

print("=" * 70)
print("VORTEX INJECTOR - FIXED DESIGN COLD FLOW TEST")
print("=" * 70)
print()
print("Fixed Design Parameters (Option A):")
print(f"  LOX supply:  {fixed_params['lox_supply_p_MPa']:.2f} MPa, {fixed_params['lox_supply_T_K']:.1f} K, {fixed_params['lox_total_mdot_kg_s']:.2f} kg/s")
print(f"  Fuel supply: {fixed_params['fuel_supply_p_MPa']:.2f} MPa, {fixed_params['fuel_supply_T_K']:.1f} K, {fixed_params['fuel_total_mdot_kg_s']:.2f} kg/s")
print(f"  Chamber:     {fixed_params['chamber_p_MPa']:.2f} MPa")
print()
print("Pressure feasibility:")
print(f"  LOX:  {data['fixed_designs']['A']['results']['lox']['required_dp_MPa']:.3f} MPa required, {data['fixed_designs']['A']['results']['lox']['available_dp_MPa']:.3f} MPa available (+{data['fixed_designs']['A']['results']['lox']['margin_percent']:.0f}%)")
print(f"  Fuel: {data['fixed_designs']['A']['results']['fuel']['required_dp_MPa']:.3f} MPa required, {data['fixed_designs']['A']['results']['fuel']['available_dp_MPa']:.3f} MPa available (+{data['fixed_designs']['A']['results']['fuel']['margin_percent']:.0f}%)")
print()

# Generate RFAM test script for A406
rfam_script = f"""#!/bin/bash
# RFAM Cold Flow Test with FIXED Design
cd /root/fine
export OMP_NUM_THREADS=8
export RFAM_TETRA_TABLE=/root/rfam-tetra/work/a406_rootcause/property_table_a.bin

./rfam_v20_tetra \\
  --mesh /root/rfam-tetra/work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet \\
  --steps 100 --monitor-every 10 \\
  --p0 {fixed_params['chamber_p_MPa']*1e6:.0f} --T0 400 \\
  --cfl 0.8 \\
  --lusgs --lusgs-relax 0.7 \\
  --outlet-p {fixed_params['chamber_p_MPa']*1e6:.0f} \\
  --backflow-T 400 --backflow-yO2 .784 \\
  --lox-pinj {fixed_params['lox_supply_p_MPa']*1e6:.0f} \\
  --fuel-pinj {fixed_params['fuel_supply_p_MPa']*1e6:.0f} \\
  --lox-T {fixed_params['lox_supply_T_K']:.1f} \\
  --fuel-T {fixed_params['fuel_supply_T_K']:.1f} \\
  --lox-mdot {fixed_params['lox_total_mdot_kg_s']:.2f} \\
  --fuel-mdot {fixed_params['fuel_total_mdot_kg_s']:.2f} \\
  --dump /tmp/rfam_fixed_design_s100.rfamtetf \\
  --dump-every 20 \\
  2>&1 | tee /tmp/rfam_fixed_design.log

echo "Exit code: $?"
"""

# Save script
with open("D:/realFluidAllMach/work/run_fixed_design_a406.sh", 'w') as f:
    f.write(rfam_script)

print("Generated A406 run script: D:/realFluidAllMach/work/run_fixed_design_a406.sh")
print()
print("To run on A406:")
print("  1. Transfer script: scp work/run_fixed_design_a406.sh root@yun.1617k.com:/tmp/")
print("  2. Execute: ssh root@yun.1617k.com 'bash /tmp/run_fixed_design_a406.sh'")
print()
print("Expected behavior with fixed design:")
print("  - Fuel inlet should no longer show entropy deficit warnings")
print("  - Both inlet pressures are now sufficient for design flow")
print("  - Convergence should be limited only by acoustic timestep (not pressure)")
print("  - Mass balance should close (inlet = outlet)")
print()
print("Note: Acoustic timestep limitation still exists (46 ps CFL limit)")
print("      This run will test if pressure fix resolves boundary compatibility")
print()
