#!/usr/bin/env python3
"""
Vortex Injector Redesign for Fixed 0.5 MPa Fuel Pressure Drop
Constraint: Regenerative cooling system limits fuel dP to 0.5 MPa
Solution: Enlarge fuel orifices to reduce pressure drop
"""

import numpy as np
import json

def calculate_required_area(mdot_total, dp_available, rho, K_total, n_pintle):
    """
    Calculate required orifice area for given pressure drop
    From: dP = K * ρ * (v)^2 / 2, where v = mdot/(ρ*A)
    Solve for A: A = mdot / sqrt(2*dP*ρ/K)
    """
    mdot_per_pintle = mdot_total / n_pintle
    # Add 20% safety margin by designing for 80% of available dP
    dp_design = dp_available * 0.8
    
    # A = mdot / (rho * sqrt(2*dp/(K*rho)))
    required_area = mdot_per_pintle / (rho * np.sqrt(2 * dp_design / (K_total * rho)))
    return required_area

def redesign_injector():
    """Redesign injector to work with fixed 0.5 MPa fuel pressure drop"""
    
    # System constraints (FIXED by regen cooling)
    fuel_dp_available = 0.5e6  # Pa (CANNOT CHANGE)
    fuel_supply_p = 4.5e6  # Pa
    chamber_p = 4.0e6  # Pa
    
    # Design requirements
    fuel_mdot_total = 1.05  # kg/s
    n_pintle = 3
    
    # Fluid properties
    fuel_rho = 65.5  # kg/m3 at 4.5 MPa, 453.5 K
    
    # Loss coefficients (conservative estimates)
    K_manifold_tee = 0.9
    K_orifice = 2.5  # Sharp-edged orifice
    K_sudden_expansion = 1.0
    K_total = K_manifold_tee + K_orifice + K_sudden_expansion
    
    # Current design
    current_area_mm2 = 62.6  # mm2 per pintle
    current_area_m2 = current_area_mm2 * 1e-6
    
    print("=" * 70)
    print("VORTEX INJECTOR REDESIGN FOR FIXED FUEL PRESSURE DROP")
    print("=" * 70)
    print()
    
    print("SYSTEM CONSTRAINTS (Fixed by Regenerative Cooling):")
    print(f"  Fuel supply pressure: {fuel_supply_p/1e6:.2f} MPa")
    print(f"  Chamber pressure:     {chamber_p/1e6:.2f} MPa")
    print(f"  Available fuel dP:    {fuel_dp_available/1e6:.2f} MPa (CANNOT CHANGE)")
    print()
    
    print("DESIGN REQUIREMENTS:")
    print(f"  Fuel flow rate:       {fuel_mdot_total:.2f} kg/s total")
    print(f"  Per pintle:           {fuel_mdot_total/n_pintle:.3f} kg/s")
    print(f"  Number of pintles:    {n_pintle}")
    print()
    
    # Calculate required area
    required_area_m2 = calculate_required_area(
        fuel_mdot_total, fuel_dp_available, fuel_rho, K_total, n_pintle
    )
    required_area_mm2 = required_area_m2 * 1e6
    
    # Calculate area increase
    area_ratio = required_area_mm2 / current_area_mm2
    area_increase_pct = (area_ratio - 1) * 100
    
    print("CURRENT DESIGN:")
    print(f"  Fuel annulus area:    {current_area_mm2:.1f} mm2 per pintle")
    print(f"  Equivalent diameter:  {np.sqrt(4*current_area_mm2/np.pi):.1f} mm")
    print()
    
    print("REQUIRED DESIGN (for 0.5 MPa dP with 20% safety margin):")
    print(f"  Fuel annulus area:    {required_area_mm2:.1f} mm2 per pintle")
    print(f"  Equivalent diameter:  {np.sqrt(4*required_area_mm2/np.pi):.1f} mm")
    print(f"  Area increase:        +{area_increase_pct:.0f}%")
    print()
    
    # Verify the design
    mdot_per_pintle = fuel_mdot_total / n_pintle
    v = mdot_per_pintle / (fuel_rho * required_area_m2)
    dp_actual = K_total * fuel_rho * v**2 / 2.0
    margin = (fuel_dp_available - dp_actual) / fuel_dp_available * 100
    
    print("VERIFICATION:")
    print(f"  Flow velocity:        {v:.1f} m/s")
    print(f"  Pressure drop:        {dp_actual/1e6:.3f} MPa")
    print(f"  Safety margin:        {margin:.0f}%")
    print()
    
    # Practical implementation suggestions
    print("=" * 70)
    print("IMPLEMENTATION OPTIONS:")
    print("=" * 70)
    print()
    
    # Option 1: Annular gap enlargement
    print("Option 1: Enlarge Fuel Annulus Gap")
    # Assuming annulus between inner and outer diameters
    # Current approximate dimensions (estimate from area)
    d_equiv_current = np.sqrt(4 * current_area_mm2 / np.pi)
    d_equiv_required = np.sqrt(4 * required_area_mm2 / np.pi)
    
    print(f"  Current equivalent diameter: {d_equiv_current:.1f} mm")
    print(f"  Required equivalent diameter: {d_equiv_required:.1f} mm")
    print(f"  Increase diameter by: {d_equiv_required - d_equiv_current:.1f} mm")
    print()
    print("  Practical implementation:")
    print("  - If annular gap h = 1.0 mm, increase to ~1.6 mm")
    print("  - Or increase annular diameter proportionally")
    print("  - Maintain pintle structural integrity")
    print()
    
    # Option 2: Multiple orifices
    print("Option 2: Multiple Smaller Orifices per Pintle")
    # Each orifice area
    n_orifices_per_pintle = 4
    area_per_orifice = required_area_mm2 / n_orifices_per_pintle
    d_orifice = np.sqrt(4 * area_per_orifice / np.pi)
    print(f"  Number of orifices per pintle: {n_orifices_per_pintle}")
    print(f"  Area per orifice: {area_per_orifice:.1f} mm2")
    print(f"  Diameter per orifice: {d_orifice:.1f} mm")
    print()
    print("  Practical implementation:")
    print("  - Drill 4 holes around pintle circumference")
    print("  - Better atomization than single annulus")
    print("  - May improve mixing")
    print()
    
    # Option 3: Reduce loss coefficients
    print("Option 3: Optimize Flow Path (Reduce K factors)")
    print("  Current total K: 4.4 (tee 0.9 + orifice 2.5 + expansion 1.0)")
    print()
    print("  Potential improvements:")
    print("  - Streamline manifold junction: K_tee 0.9 → 0.4 (save 0.5)")
    print("  - Rounded orifice entrance: K_orifice 2.5 → 1.5 (save 1.0)")
    print("  - Gradual expansion: K_expansion 1.0 → 0.3 (save 0.7)")
    print("  - New total K: ~2.2 (50% reduction)")
    print()
    
    # Calculate area needed with optimized K
    K_optimized = 2.2
    required_area_opt_m2 = calculate_required_area(
        fuel_mdot_total, fuel_dp_available, fuel_rho, K_optimized, n_pintle
    )
    required_area_opt_mm2 = required_area_opt_m2 * 1e6
    area_increase_opt_pct = (required_area_opt_mm2 / current_area_mm2 - 1) * 100
    
    print(f"  With optimized K=2.2:")
    print(f"  - Required area: {required_area_opt_mm2:.1f} mm2 (+{area_increase_opt_pct:.0f}%)")
    print(f"  - Combined with Option 1: less enlargement needed")
    print()
    
    # Save results
    results = {
        "constraint": {
            "fuel_dp_available_MPa": fuel_dp_available / 1e6,
            "reason": "Fixed by regenerative cooling system",
            "fuel_supply_p_MPa": fuel_supply_p / 1e6,
            "chamber_p_MPa": chamber_p / 1e6
        },
        "current_design": {
            "fuel_area_mm2_per_pintle": current_area_mm2,
            "equivalent_diameter_mm": float(d_equiv_current),
            "feasible": False,
            "reason": "Requires 1.05 MPa, only 0.5 MPa available"
        },
        "redesign_option_1": {
            "name": "Enlarge annular area",
            "required_area_mm2_per_pintle": float(required_area_mm2),
            "equivalent_diameter_mm": float(d_equiv_required),
            "area_increase_percent": float(area_increase_pct),
            "pressure_drop_MPa": dp_actual / 1e6,
            "safety_margin_percent": float(margin),
            "feasible": True
        },
        "redesign_option_2": {
            "name": "Multiple orifices per pintle",
            "orifices_per_pintle": n_orifices_per_pintle,
            "area_per_orifice_mm2": float(area_per_orifice),
            "diameter_per_orifice_mm": float(d_orifice),
            "feasible": True
        },
        "redesign_option_3": {
            "name": "Optimize flow path (reduce K)",
            "K_current": K_total,
            "K_optimized": K_optimized,
            "required_area_mm2_with_optimization": float(required_area_opt_mm2),
            "area_increase_percent": float(area_increase_opt_pct),
            "feasible": True
        },
        "recommendation": "Combine Option 1 (enlarge area) + Option 3 (optimize flow path) for minimum geometry change"
    }
    
    with open("D:/realFluidAllMach/work/injector_redesign_fixed_dp.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    print("=" * 70)
    print("RECOMMENDATION:")
    print("=" * 70)
    print()
    print("BEST APPROACH: Combination Strategy")
    print("  1. Enlarge fuel annulus area by ~40-50% (not full 59%)")
    print("  2. Streamline flow path to reduce K by ~30%")
    print("  3. Result: Achieves 0.5 MPa target with moderate geometry changes")
    print()
    print("CRITICAL: We CANNOT change the 0.5 MPa pressure drop")
    print("          This is constrained by the regenerative cooling system")
    print("          The injector geometry must adapt to this constraint")
    print()
    print("Results saved to: D:/realFluidAllMach/work/injector_redesign_fixed_dp.json")
    print()
    
    return results

if __name__ == "__main__":
    results = redesign_injector()
    
    print()
    print("NEXT STEP: Run CFD with redesigned geometry")
    print("  - Update CAD model with larger fuel orifices")
    print("  - Remesh the LOX domain")
    print("  - Run cold flow with 0.5 MPa fuel dP")
    print()
