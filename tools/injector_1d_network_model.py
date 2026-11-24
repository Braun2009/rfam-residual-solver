#!/usr/bin/env python3
"""
1-D Network Model for Vortex Triple-Pintle Injector
Estimates 3-way flow split with manifold junction losses

Purpose: Interim design guidance while detailed 3-D CFD is blocked by:
- RFAM acoustic timestep limitation
- Fluent boundary zone grouping issue

Model: Shared LOX manifold → 3 parallel pintle branches
Each branch: orifice → pintle gap → chamber
Assumes: incompressible, isothermal, fully-developed losses
"""

import numpy as np
import json
from scipy.optimize import fsolve
import sys

# Design parameters from Vortex configuration
DESIGN = {
    "lox_supply_p_MPa": 6.0,
    "lox_supply_T_K": 133.5,
    "lox_total_mdot_kg_s": 2.8,
    "chamber_p_MPa": 4.0,
    
    "fuel_supply_p_MPa": 4.5,
    "fuel_supply_T_K": 453.5,
    "fuel_total_mdot_kg_s": 1.05,
    
    # Pintle geometry (from CAD)
    "pintle_count": 3,
    "lox_orifice_diameter_mm": 5.0,  # Per V6 design
    "lox_slot_area_mm2": 49.5,  # Per zone from mesh analysis
    
    "fuel_orifice_diameter_mm": 4.5,
    "fuel_annulus_area_mm2": 62.6,
    
    # Loss coefficients (typical for sharp-edged orifices and gaps)
    "K_orifice": 2.5,  # Sharp-edged orifice
    "K_sudden_expansion": 1.0,  # Pintle gap → chamber
    "K_manifold_tee": 0.9,  # Manifold branch junction
    
    # Fluid properties (from CoolProp @ supply conditions)
    "lox_rho_kg_m3": 909.3,
    "fuel_rho_kg_m3": 65.5,
}

def pressure_drop_incompressible(mdot, area_m2, rho_kg_m3, K_loss):
    """
    Incompressible pressure drop: ΔP = K * (ρ * v²/2)
    where v = mdot / (rho * A)
    """
    if mdot <= 0 or area_m2 <= 0:
        return 0.0
    v = mdot / (rho_kg_m3 * area_m2)
    return K_loss * rho_kg_m3 * v**2 / 2.0

def pintle_pressure_drop(mdot, params, propellant="lox"):
    """
    Total pressure drop through one pintle branch:
    1. Manifold tee junction loss
    2. Orifice/slot loss
    3. Sudden expansion to chamber loss
    """
    if propellant == "lox":
        area = params["lox_slot_area_mm2"] * 1e-6
        rho = params["lox_rho_kg_m3"]
        K_total = params["K_manifold_tee"] + params["K_orifice"] + params["K_sudden_expansion"]
    else:  # fuel
        area = params["fuel_annulus_area_mm2"] * 1e-6
        rho = params["fuel_rho_kg_m3"]
        K_total = params["K_manifold_tee"] + params["K_orifice"] + params["K_sudden_expansion"]
    
    return pressure_drop_incompressible(mdot, area, rho, K_total)

def solve_flow_split(params):
    """
    Solve for flow split assuming equal pressure drop through all branches.
    
    Unknowns: mdot1, mdot2, mdot3 (kg/s per pintle)
    Constraints:
    1. mdot1 + mdot2 + mdot3 = total_mdot
    2. ΔP1 = ΔP2 = ΔP3 = supply_p - chamber_p
    
    For symmetric geometry → mdot1 = mdot2 = mdot3
    """
    
    propellants = ["lox", "fuel"]
    results = {}
    
    for prop in propellants:
        if prop == "lox":
            supply_p = params["lox_supply_p_MPa"] * 1e6
            chamber_p = params["chamber_p_MPa"] * 1e6
            total_mdot = params["lox_total_mdot_kg_s"]
        else:
            supply_p = params["fuel_supply_p_MPa"] * 1e6
            chamber_p = params["chamber_p_MPa"] * 1e6
            total_mdot = params["fuel_total_mdot_kg_s"]
        
        available_dp = supply_p - chamber_p
        n_pintle = params["pintle_count"]
        
        # For symmetric geometry, each pintle gets equal flow
        mdot_per_pintle = total_mdot / n_pintle
        
        # Check if pressure drop is consistent
        dp_required = pintle_pressure_drop(mdot_per_pintle, params, prop)
        
        # Calculate velocity and pressure loss breakdown
        if prop == "lox":
            area = params["lox_slot_area_mm2"] * 1e-6
            rho = params["lox_rho_kg_m3"]
        else:
            area = params["fuel_annulus_area_mm2"] * 1e-6
            rho = params["fuel_rho_kg_m3"]
        
        velocity = mdot_per_pintle / (rho * area)
        dynamic_pressure = rho * velocity**2 / 2.0
        
        # Individual loss components
        dp_manifold = params["K_manifold_tee"] * dynamic_pressure
        dp_orifice = params["K_orifice"] * dynamic_pressure
        dp_expansion = params["K_sudden_expansion"] * dynamic_pressure
        
        results[prop] = {
            "total_mdot_kg_s": total_mdot,
            "mdot_per_pintle_kg_s": mdot_per_pintle,
            "split_percentage": [100.0/n_pintle] * n_pintle,
            "available_dp_Pa": available_dp,
            "required_dp_Pa": dp_required,
            "dp_margin_percent": (available_dp - dp_required) / available_dp * 100,
            "velocity_m_s": velocity,
            "dynamic_pressure_Pa": dynamic_pressure,
            "breakdown": {
                "manifold_tee_Pa": dp_manifold,
                "orifice_Pa": dp_orifice,
                "sudden_expansion_Pa": dp_expansion,
            }
        }
        
        # Check feasibility
        if dp_required > available_dp:
            results[prop]["feasible"] = False
            results[prop]["message"] = f"Required ΔP ({dp_required/1e6:.2f} MPa) exceeds available ({available_dp/1e6:.2f} MPa)"
        else:
            results[prop]["feasible"] = True
            results[prop]["message"] = "Flow is feasible"
    
    return results

def asymmetric_split_sensitivity(params, pintle_variations):
    """
    Sensitivity analysis: what if one pintle has different area/loss?
    """
    # TODO: Implement non-uniform pintle analysis
    # This would solve: ΔP1 = ΔP2 = ΔP3 with different K or A values
    pass

def main():
    print("=" * 70)
    print("Vortex Triple-Pintle Injector - 1-D Network Model")
    print("=" * 70)
    print()
    
    print("Design Parameters:")
    print(f"  LOX: {DESIGN['lox_supply_p_MPa']:.1f} MPa, {DESIGN['lox_supply_T_K']:.1f} K, {DESIGN['lox_total_mdot_kg_s']:.2f} kg/s")
    print(f"  Fuel: {DESIGN['fuel_supply_p_MPa']:.1f} MPa, {DESIGN['fuel_supply_T_K']:.1f} K, {DESIGN['fuel_total_mdot_kg_s']:.2f} kg/s")
    print(f"  Chamber: {DESIGN['chamber_p_MPa']:.1f} MPa")
    print(f"  Pintles: {DESIGN['pintle_count']}")
    print()
    
    # Solve symmetric case
    results = solve_flow_split(DESIGN)
    
    print("Results (Symmetric Geometry Assumption):")
    print()
    
    for prop, res in results.items():
        print(f"{prop.upper()} Flow Split:")
        print(f"  Total flow: {res['total_mdot_kg_s']:.3f} kg/s")
        print(f"  Per pintle: {res['mdot_per_pintle_kg_s']:.3f} kg/s ({res['split_percentage'][0]:.1f}%)")
        print(f"  Available ΔP: {res['available_dp_Pa']/1e6:.3f} MPa")
        print(f"  Required ΔP: {res['required_dp_Pa']/1e6:.3f} MPa")
        print(f"  Margin: {res['dp_margin_percent']:.1f}%")
        print(f"  Velocity: {res['velocity_m_s']:.1f} m/s")
        print(f"  Feasible: {res['feasible']} - {res['message']}")
        print(f"  Loss breakdown:")
        for component, dp in res['breakdown'].items():
            print(f"    {component}: {dp/1e6:.3f} MPa ({dp/res['required_dp_Pa']*100:.1f}%)")
        print()
    
    # Export results
    output = {
        "model": "1D_network_symmetric",
        "design_parameters": DESIGN,
        "results": results,
        "assumptions": [
            "Incompressible flow",
            "Isothermal (no heat transfer)",
            "Fully-developed losses (K factors)",
            "Symmetric geometry (equal split)",
            "No cross-flow or interaction effects",
            "Steady state"
        ],
        "limitations": [
            "Does not capture 3-D recirculation zones",
            "Does not model pintle-to-pintle interaction",
            "Ignores unsteady effects and acoustic coupling",
            "K factors are empirical estimates",
            "No validation against CFD or test data"
        ],
        "validation_status": "INTERIM_ESTIMATE_ONLY",
        "recommendation": "Use for design guidance only. Validate with 3-D CFD when solver issues are resolved."
    }
    
    output_file = "D:/realFluidAllMach/work/injector_1d_flow_split.json"
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"Results saved to: {output_file}")
    print()
    print("IMPORTANT NOTES:")
    print("  - This is a 1-D estimate based on loss coefficient correlations")
    print("  - For symmetric geometry, flow split is exactly 33.33% per pintle")
    print("  - Actual split may vary due to 3-D effects, manufacturing tolerances")
    print("  - Validate with detailed CFD when available")
    print()

if __name__ == "__main__":
    main()
