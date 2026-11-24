#!/usr/bin/env python3
"""
Vortex Injector Design Fixer
Automatically adjusts operating conditions to achieve 2.8 kg/s cold flow feasibility
"""

import numpy as np
import json
from scipy.optimize import fsolve

def pressure_drop(mdot, area_m2, rho_kg_m3, K_total):
    """ΔP = K * (ρ * v²/2) where v = mdot / (rho * A)"""
    if mdot <= 0 or area_m2 <= 0:
        return 0.0
    v = mdot / (rho_kg_m3 * area_m2)
    return K_total * rho_kg_m3 * v**2 / 2.0

def evaluate_design(params):
    """Evaluate if design is feasible"""
    results = {}
    
    for prop in ["lox", "fuel"]:
        if prop == "lox":
            supply_p = params["lox_supply_p_MPa"] * 1e6
            chamber_p = params["chamber_p_MPa"] * 1e6
            total_mdot = params["lox_total_mdot_kg_s"]
            area = params["lox_slot_area_mm2"] * 1e-6
            rho = params["lox_rho_kg_m3"]
        else:
            supply_p = params["fuel_supply_p_MPa"] * 1e6
            chamber_p = params["chamber_p_MPa"] * 1e6
            total_mdot = params["fuel_total_mdot_kg_s"]
            area = params["fuel_annulus_area_mm2"] * 1e-6
            rho = params["fuel_rho_kg_m3"]
        
        n_pintle = params["pintle_count"]
        K_total = params["K_manifold_tee"] + params["K_orifice"] + params["K_sudden_expansion"]
        
        mdot_per_pintle = total_mdot / n_pintle
        available_dp = supply_p - chamber_p
        required_dp = pressure_drop(mdot_per_pintle, area, rho, K_total)
        
        velocity = mdot_per_pintle / (rho * area)
        margin_percent = (available_dp - required_dp) / available_dp * 100
        
        results[prop] = {
            "available_dp_MPa": available_dp / 1e6,
            "required_dp_MPa": required_dp / 1e6,
            "margin_percent": margin_percent,
            "velocity_m_s": velocity,
            "feasible": required_dp <= available_dp
        }
    
    return results

def fix_option_a_increase_fuel_pressure(base_params):
    """Option A: Increase fuel supply pressure"""
    params = base_params.copy()
    
    # Calculate required fuel supply pressure
    chamber_p = params["chamber_p_MPa"] * 1e6
    total_mdot = params["fuel_total_mdot_kg_s"]
    area = params["fuel_annulus_area_mm2"] * 1e-6
    rho = params["fuel_rho_kg_m3"]
    n_pintle = params["pintle_count"]
    K_total = params["K_manifold_tee"] + params["K_orifice"] + params["K_sudden_expansion"]
    
    mdot_per_pintle = total_mdot / n_pintle
    required_dp = pressure_drop(mdot_per_pintle, area, rho, K_total)
    
    # Add 20% safety margin
    required_dp_with_margin = required_dp * 1.2
    required_supply_p = (chamber_p + required_dp_with_margin) / 1e6
    
    params["fuel_supply_p_MPa"] = required_supply_p
    
    return params, f"Increase fuel supply pressure to {required_supply_p:.2f} MPa (was 4.5 MPa)"

def fix_option_c_enlarge_fuel_orifice(base_params):
    """Option C: Enlarge fuel orifice area"""
    params = base_params.copy()
    
    # Calculate required area
    supply_p = params["fuel_supply_p_MPa"] * 1e6
    chamber_p = params["chamber_p_MPa"] * 1e6
    available_dp = supply_p - chamber_p
    
    total_mdot = params["fuel_total_mdot_kg_s"]
    rho = params["fuel_rho_kg_m3"]
    n_pintle = params["pintle_count"]
    K_total = params["K_manifold_tee"] + params["K_orifice"] + params["K_sudden_expansion"]
    
    mdot_per_pintle = total_mdot / n_pintle
    
    # From ΔP = K * ρ * (mdot/(rho*A))² / 2
    # Solve for A: A = mdot / sqrt(2 * ΔP / K)
    # Add 20% margin by reducing ΔP used
    design_dp = available_dp / 1.2
    
    required_area = mdot_per_pintle / (rho * np.sqrt(2 * design_dp / (K_total * rho)))
    required_area_mm2 = required_area * 1e6
    
    params["fuel_annulus_area_mm2"] = required_area_mm2
    
    original_area = base_params["fuel_annulus_area_mm2"]
    area_increase = (required_area_mm2 / original_area - 1) * 100
    
    # Estimate new diameter (assuming annular area)
    # For annulus: A = π/4 * (D_outer² - D_inner²)
    # Rough estimate: scale diameter by sqrt(area_ratio)
    original_d_equiv = np.sqrt(4 * original_area / np.pi)
    new_d_equiv = np.sqrt(4 * required_area_mm2 / np.pi)
    
    return params, f"Increase fuel annulus area from {original_area:.1f} mm2 to {required_area_mm2:.1f} mm2 (+{area_increase:.0f}%), equivalent diameter {original_d_equiv:.1f}->{new_d_equiv:.1f} mm"

def fix_option_d_reduce_chamber_pressure(base_params):
    """Option D: Reduce chamber pressure"""
    params = base_params.copy()
    
    # Calculate required chamber pressure
    fuel_supply_p = params["fuel_supply_p_MPa"] * 1e6
    total_mdot = params["fuel_total_mdot_kg_s"]
    area = params["fuel_annulus_area_mm2"] * 1e-6
    rho = params["fuel_rho_kg_m3"]
    n_pintle = params["pintle_count"]
    K_total = params["K_manifold_tee"] + params["K_orifice"] + params["K_sudden_expansion"]
    
    mdot_per_pintle = total_mdot / n_pintle
    required_dp = pressure_drop(mdot_per_pintle, area, rho, K_total)
    
    # Add 20% margin
    required_dp_with_margin = required_dp * 1.2
    max_chamber_p = (fuel_supply_p - required_dp_with_margin) / 1e6
    
    params["chamber_p_MPa"] = max_chamber_p
    
    pressure_reduction = (base_params["chamber_p_MPa"] - max_chamber_p) / base_params["chamber_p_MPa"] * 100
    
    return params, f"Reduce chamber pressure to {max_chamber_p:.2f} MPa (was 4.0 MPa, -{pressure_reduction:.0f}%)"

def main():
    # Original design
    base_design = {
        "lox_supply_p_MPa": 6.0,
        "lox_supply_T_K": 133.5,
        "lox_total_mdot_kg_s": 2.8,
        "chamber_p_MPa": 4.0,
        
        "fuel_supply_p_MPa": 4.5,
        "fuel_supply_T_K": 453.5,
        "fuel_total_mdot_kg_s": 1.05,
        
        "pintle_count": 3,
        "lox_slot_area_mm2": 49.5,
        "fuel_annulus_area_mm2": 62.6,
        
        "K_orifice": 2.5,
        "K_sudden_expansion": 1.0,
        "K_manifold_tee": 0.9,
        
        "lox_rho_kg_m3": 909.3,
        "fuel_rho_kg_m3": 65.5,
    }
    
    print("=" * 70)
    print("VORTEX INJECTOR DESIGN FIXER")
    print("=" * 70)
    print()
    
    print("Original Design:")
    original_results = evaluate_design(base_design)
    print(f"  LOX:  {original_results['lox']['required_dp_MPa']:.3f} MPa required, {original_results['lox']['available_dp_MPa']:.3f} MPa available -> {'OK' if original_results['lox']['feasible'] else 'FAIL'} ({original_results['lox']['margin_percent']:+.0f}%)")
    print(f"  Fuel: {original_results['fuel']['required_dp_MPa']:.3f} MPa required, {original_results['fuel']['available_dp_MPa']:.3f} MPa available -> {'OK' if original_results['fuel']['feasible'] else 'FAIL'} ({original_results['fuel']['margin_percent']:+.0f}%)")
    print()
    
    # Generate fix options
    fixes = [
        ("A", "Increase Fuel Supply Pressure", fix_option_a_increase_fuel_pressure),
        ("C", "Enlarge Fuel Orifice Area", fix_option_c_enlarge_fuel_orifice),
        ("D", "Reduce Chamber Pressure", fix_option_d_reduce_chamber_pressure),
    ]
    
    fixed_designs = {}
    
    for option, name, fix_func in fixes:
        print(f"─── Option {option}: {name} ───")
        fixed_params, description = fix_func(base_design)
        print(f"  Change: {description}")
        
        results = evaluate_design(fixed_params)
        print(f"  LOX:  {results['lox']['required_dp_MPa']:.3f} MPa required, {results['lox']['available_dp_MPa']:.3f} MPa available -> {'OK' if results['lox']['feasible'] else 'FAIL'} ({results['lox']['margin_percent']:+.0f}%)")
        print(f"  Fuel: {results['fuel']['required_dp_MPa']:.3f} MPa required, {results['fuel']['available_dp_MPa']:.3f} MPa available -> {'OK' if results['fuel']['feasible'] else 'FAIL'} ({results['fuel']['margin_percent']:+.0f}%)")
        print()
        
        fixed_designs[option] = {
            "name": name,
            "description": description,
            "parameters": fixed_params,
            "results": {
                "lox": {
                    "available_dp_MPa": results['lox']['available_dp_MPa'],
                    "required_dp_MPa": results['lox']['required_dp_MPa'],
                    "margin_percent": results['lox']['margin_percent'],
                    "velocity_m_s": results['lox']['velocity_m_s'],
                    "feasible": bool(results['lox']['feasible'])
                },
                "fuel": {
                    "available_dp_MPa": results['fuel']['available_dp_MPa'],
                    "required_dp_MPa": results['fuel']['required_dp_MPa'],
                    "margin_percent": results['fuel']['margin_percent'],
                    "velocity_m_s": results['fuel']['velocity_m_s'],
                    "feasible": bool(results['fuel']['feasible'])
                }
            },
            "both_feasible": bool(results['lox']['feasible'] and results['fuel']['feasible'])
        }
    
    # Export all options
    output = {
        "original_design": {
            "parameters": base_design,
            "results": {
                "lox": {
                    "available_dp_MPa": original_results['lox']['available_dp_MPa'],
                    "required_dp_MPa": original_results['lox']['required_dp_MPa'],
                    "margin_percent": original_results['lox']['margin_percent'],
                    "velocity_m_s": original_results['lox']['velocity_m_s'],
                    "feasible": bool(original_results['lox']['feasible'])
                },
                "fuel": {
                    "available_dp_MPa": original_results['fuel']['available_dp_MPa'],
                    "required_dp_MPa": original_results['fuel']['required_dp_MPa'],
                    "margin_percent": original_results['fuel']['margin_percent'],
                    "velocity_m_s": original_results['fuel']['velocity_m_s'],
                    "feasible": bool(original_results['fuel']['feasible'])
                }
            },
            "feasible": bool(original_results['lox']['feasible'] and original_results['fuel']['feasible'])
        },
        "fixed_designs": fixed_designs,
        "recommendation": "Option A (increase fuel supply pressure) is recommended - minimal system changes"
    }
    
    with open("D:/realFluidAllMach/work/injector_design_fixed.json", 'w') as f:
        json.dump(output, f, indent=2)
    
    print("===" + "=" * 67)
    print("RECOMMENDATION: Option A (Increase Fuel Supply Pressure)")
    print("===" + "=" * 67)
    print(f"  New fuel supply: {fixed_designs['A']['parameters']['fuel_supply_p_MPa']:.2f} MPa (was 4.5 MPa)")
    print(f"  All other parameters unchanged")
    print(f"  Both propellants now feasible with safety margins")
    print()
    print("Results saved to: D:/realFluidAllMach/work/injector_design_fixed.json")
    print()
    
    return fixed_designs['A']['parameters']

if __name__ == "__main__":
    fixed_params = main()
