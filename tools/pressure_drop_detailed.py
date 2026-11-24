#!/usr/bin/env python3
"""
Corrected pressure drop calculation with detailed step-by-step process
"""

import numpy as np

print("=" * 70)
print("DETAILED PRESSURE DROP CALCULATION - STEP BY STEP")
print("=" * 70)
print()

# ============================================================================
# GIVEN DATA
# ============================================================================
print("GIVEN CONDITIONS:")
print("-" * 70)

# LOX side
lox_supply_p = 6.0e6  # Pa
lox_density = 909.3  # kg/m3
lox_flow_total = 2.8  # kg/s
n_pintles = 3
lox_flow_per_pintle = lox_flow_total / n_pintles

print(f"LOX supply pressure:  {lox_supply_p/1e6:.2f} MPa")
print(f"LOX density:          {lox_density:.1f} kg/m3")
print(f"LOX total flow:       {lox_flow_total:.2f} kg/s")
print(f"LOX per pintle:       {lox_flow_per_pintle:.3f} kg/s")

# Fuel side
fuel_supply_p = 4.5e6  # Pa
fuel_density = 65.5  # kg/m3
fuel_flow_total = 1.05  # kg/s
fuel_flow_per_pintle = fuel_flow_total / n_pintles

print(f"Fuel supply pressure: {fuel_supply_p/1e6:.2f} MPa")
print(f"Fuel density:         {fuel_density:.1f} kg/m3")
print(f"Fuel total flow:      {fuel_flow_total:.2f} kg/s")
print(f"Fuel per pintle:      {fuel_flow_per_pintle:.3f} kg/s")

# Chamber
chamber_p = 4.0e6  # Pa
print(f"Chamber pressure:     {chamber_p/1e6:.2f} MPa")

print()

# ============================================================================
# AVAILABLE PRESSURE DROP
# ============================================================================
print("AVAILABLE PRESSURE DROP:")
print("-" * 70)

lox_dp_available = lox_supply_p - chamber_p
fuel_dp_available = fuel_supply_p - chamber_p

print(f"LOX:  {lox_supply_p/1e6:.2f} - {chamber_p/1e6:.2f} = {lox_dp_available/1e6:.2f} MPa available")
print(f"Fuel: {fuel_supply_p/1e6:.2f} - {chamber_p/1e6:.2f} = {fuel_dp_available/1e6:.2f} MPa available")
print()

# ============================================================================
# CURRENT DESIGN (from mesh analysis)
# ============================================================================
print("CURRENT INJECTOR GEOMETRY (from existing design):")
print("-" * 70)

# From your mesh: work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet
# LOX zones: 10, 11, 12 with areas ~49.5 mm2 each
lox_area_per_pintle_current = 49.5  # mm2
lox_area_per_pintle_m2 = lox_area_per_pintle_current * 1e-6

# Fuel zones: 20, 21, 22 with areas ~62.6 mm2 each  
fuel_area_per_pintle_current = 62.6  # mm2
fuel_area_per_pintle_m2 = fuel_area_per_pintle_current * 1e-6

print(f"LOX area per pintle:  {lox_area_per_pintle_current:.1f} mm2")
print(f"Fuel area per pintle: {fuel_area_per_pintle_current:.1f} mm2")
print()

# ============================================================================
# CALCULATE VELOCITIES
# ============================================================================
print("STEP 1: CALCULATE FLOW VELOCITIES")
print("-" * 70)
print("Formula: v = mdot / (rho * A)")
print()

# LOX velocity
lox_velocity = lox_flow_per_pintle / (lox_density * lox_area_per_pintle_m2)
print(f"LOX: v = {lox_flow_per_pintle:.4f} / ({lox_density:.1f} * {lox_area_per_pintle_m2:.6f})")
print(f"LOX velocity = {lox_velocity:.2f} m/s")
print()

# Fuel velocity
fuel_velocity = fuel_flow_per_pintle / (fuel_density * fuel_area_per_pintle_m2)
print(f"Fuel: v = {fuel_flow_per_pintle:.4f} / ({fuel_density:.1f} * {fuel_area_per_pintle_m2:.6f})")
print(f"Fuel velocity = {fuel_velocity:.2f} m/s")
print()

# ============================================================================
# CALCULATE DYNAMIC PRESSURE
# ============================================================================
print("STEP 2: CALCULATE DYNAMIC PRESSURE")
print("-" * 70)
print("Formula: q = 0.5 * rho * v^2")
print()

lox_dynamic_p = 0.5 * lox_density * lox_velocity**2
print(f"LOX: q = 0.5 * {lox_density:.1f} * {lox_velocity:.2f}^2")
print(f"LOX dynamic pressure = {lox_dynamic_p/1e6:.4f} MPa = {lox_dynamic_p/1e3:.2f} kPa")
print()

fuel_dynamic_p = 0.5 * fuel_density * fuel_velocity**2
print(f"Fuel: q = 0.5 * {fuel_density:.1f} * {fuel_velocity:.2f}^2")
print(f"Fuel dynamic pressure = {fuel_dynamic_p/1e6:.4f} MPa = {fuel_dynamic_p/1e3:.2f} kPa")
print()

# ============================================================================
# K FACTORS (CORRECTED)
# ============================================================================
print("STEP 3: ESTIMATE LOSS COEFFICIENTS (K factors)")
print("-" * 70)
print("Formula: dP = K * q = K * (0.5 * rho * v^2)")
print()

print("K factor components:")
print("  Manifold distribution: K = 0.1-0.2 (large plenum)")
print("  Feed tubes:            K = 0.2-0.4 (short, smooth)")
print("  Orifice:               K = 1.5-2.5 (sharp edge to rounded)")
print("  Total typical:         K = 2.0-3.0")
print()

# Conservative estimates for current geometry
lox_K_total = 2.5  # Manifold + tubes + orifice
fuel_K_total = 4.4  # Higher due to more complex path

print(f"Assumed K_lox:  {lox_K_total}")
print(f"Assumed K_fuel: {fuel_K_total}")
print()

# ============================================================================
# CALCULATE PRESSURE DROPS
# ============================================================================
print("STEP 4: CALCULATE PRESSURE DROP")
print("-" * 70)

lox_dp_calculated = lox_K_total * lox_dynamic_p
print(f"LOX: dP = {lox_K_total} * {lox_dynamic_p/1e6:.4f} MPa")
print(f"LOX pressure drop = {lox_dp_calculated/1e6:.3f} MPa")
print()

fuel_dp_calculated = fuel_K_total * fuel_dynamic_p
print(f"Fuel: dP = {fuel_K_total} * {fuel_dynamic_p/1e6:.4f} MPa")
print(f"Fuel pressure drop = {fuel_dp_calculated/1e6:.3f} MPa")
print()

# ============================================================================
# COMPARE TO AVAILABLE
# ============================================================================
print("=" * 70)
print("COMPARISON: REQUIRED vs AVAILABLE")
print("=" * 70)
print()

print(f"LOX:")
print(f"  Available: {lox_dp_available/1e6:.3f} MPa")
print(f"  Required:  {lox_dp_calculated/1e6:.3f} MPa")
if lox_dp_calculated <= lox_dp_available:
    margin = (lox_dp_available - lox_dp_calculated) / lox_dp_available * 100
    print(f"  Status: OK (+{margin:.0f}% margin)")
else:
    deficit = (lox_dp_calculated - lox_dp_available) / lox_dp_available * 100
    print(f"  Status: INSUFFICIENT (-{deficit:.0f}% deficit)")
print()

print(f"Fuel:")
print(f"  Available: {fuel_dp_available/1e6:.3f} MPa")
print(f"  Required:  {fuel_dp_calculated/1e6:.3f} MPa")
if fuel_dp_calculated <= fuel_dp_available:
    margin = (fuel_dp_available - fuel_dp_calculated) / fuel_dp_available * 100
    print(f"  Status: OK (+{margin:.0f}% margin)")
else:
    deficit = (fuel_dp_calculated - fuel_dp_available) / fuel_dp_available * 100
    print(f"  Status: INSUFFICIENT (-{deficit:.0f}% deficit)")
print()

# ============================================================================
# WHAT WENT WRONG IN MY PREVIOUS CALCULATION
# ============================================================================
print("=" * 70)
print("ERROR ANALYSIS: What went wrong in previous calculation")
print("=" * 70)
print()

print("MISTAKE in parametric_injector_design.py:")
print("  I used LOX area = 19.6 mm2 (4 holes x 2.5mm diameter)")
print("  This is 2.5× SMALLER than actual 49.5 mm2")
print("  Smaller area → higher velocity → much higher pressure drop")
print()

lox_area_wrong = 19.6e-6  # m2
lox_velocity_wrong = lox_flow_per_pintle / (lox_density * lox_area_wrong)
lox_dp_wrong = lox_K_total * 0.5 * lox_density * lox_velocity_wrong**2

print(f"With WRONG area (19.6 mm2):")
print(f"  Velocity:      {lox_velocity_wrong:.1f} m/s (vs actual {lox_velocity:.1f} m/s)")
print(f"  Pressure drop: {lox_dp_wrong/1e6:.2f} MPa (vs actual {lox_dp_calculated/1e6:.2f} MPa)")
print()

print("CONCLUSION:")
print("  My parametric design had TOO SMALL orifice areas")
print("  I was designing from scratch without checking actual geometry")
print("  Should have used the existing mesh dimensions as reference")
print()

# ============================================================================
# CORRECTED APPROACH
# ============================================================================
print("=" * 70)
print("CORRECTED APPROACH")
print("=" * 70)
print()

print("What we actually need to do:")
print()
print("1. Current fuel area (62.6 mm2) gives:")
print(f"   Velocity: {fuel_velocity:.1f} m/s")
print(f"   Pressure drop: {fuel_dp_calculated/1e6:.3f} MPa")
print(f"   But only {fuel_dp_available/1e6:.3f} MPa available")
print(f"   Deficit: {(fuel_dp_calculated - fuel_dp_available)/1e6:.3f} MPa")
print()

print("2. Required fuel area for 0.5 MPa drop:")
# dP = K * 0.5 * rho * (mdot/(rho*A))^2
# dP = K * 0.5 * mdot^2 / (rho * A^2)
# A^2 = K * 0.5 * mdot^2 / (rho * dP)
# A = mdot * sqrt(K * 0.5 / (rho * dP))

fuel_area_required_m2 = fuel_flow_per_pintle * np.sqrt(fuel_K_total * 0.5 / (fuel_density * fuel_dp_available))
fuel_area_required_mm2 = fuel_area_required_m2 * 1e6
area_increase = (fuel_area_required_mm2 / fuel_area_per_pintle_current - 1) * 100

print(f"   Required area: {fuel_area_required_mm2:.1f} mm2")
print(f"   Current area:  {fuel_area_per_pintle_current:.1f} mm2")
print(f"   Need to increase by: {area_increase:.0f}%")
print()

print("3. Verify LOX is OK:")
print(f"   LOX available: {lox_dp_available/1e6:.3f} MPa")
print(f"   LOX required:  {lox_dp_calculated/1e6:.3f} MPa")
lox_margin = (lox_dp_available - lox_dp_calculated) / lox_dp_available * 100
print(f"   Margin: +{lox_margin:.0f}%")
print(f"   Status: LOX is fine, no change needed")
print()

# Save corrected calculation
results = {
    "current_design": {
        "lox_area_mm2": lox_area_per_pintle_current,
        "fuel_area_mm2": fuel_area_per_pintle_current,
        "lox_velocity_m_s": float(lox_velocity),
        "fuel_velocity_m_s": float(fuel_velocity),
        "lox_dp_MPa": float(lox_dp_calculated / 1e6),
        "fuel_dp_MPa": float(fuel_dp_calculated / 1e6)
    },
    "available_dp": {
        "lox_MPa": float(lox_dp_available / 1e6),
        "fuel_MPa": float(fuel_dp_available / 1e6)
    },
    "required_modification": {
        "fuel_area_required_mm2": float(fuel_area_required_mm2),
        "area_increase_percent": float(area_increase),
        "lox_modification": "None needed"
    },
    "K_factors_used": {
        "lox": lox_K_total,
        "fuel": fuel_K_total,
        "note": "These are estimates; CFD will give actual values"
    }
}

import json
with open("D:/realFluidAllMach/work/pressure_drop_corrected.json", 'w') as f:
    json.dump(results, f, indent=2)

print("Corrected results saved to: D:/realFluidAllMach/work/pressure_drop_corrected.json")
print()
