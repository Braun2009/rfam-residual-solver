#!/usr/bin/env python3
"""
Injector Design Space Exploration
Generate and evaluate multiple design variants
Optimize for: pressure drop + TMR + atomization quality
"""

import numpy as np
import pandas as pd
import json

class InjectorDesignSpace:
    """Systematic exploration of injector parameter space"""
    
    def __init__(self):
        # Fixed parameters
        self.lox_flow = 2.8  # kg/s total
        self.fuel_flow = 1.05  # kg/s total
        self.n_pintles = 3
        self.lox_density = 909.3  # kg/m3
        self.fuel_density = 65.5  # kg/m3
        self.chamber_p = 4.0e6  # Pa
        
        # Supply pressures
        self.lox_supply_p = 6.0e6  # Pa
        self.fuel_supply_p = 4.5e6  # Pa
        
        # Available pressure drops
        self.lox_dp_max = self.lox_supply_p - self.chamber_p
        self.fuel_dp_max = self.fuel_supply_p - self.chamber_p
        
        # Current baseline (from mesh)
        self.baseline_lox_area = 49.5  # mm2 per pintle
        self.baseline_fuel_area = 62.6  # mm2 per pintle
        
    def calculate_performance(self, lox_area_mm2, fuel_area_mm2, 
                            fuel_swirl_deg=0, K_lox=2.5, K_fuel=4.4):
        """
        Calculate all performance metrics for a design variant
        
        Returns: dict with pressure drops, TMR, velocities, etc.
        """
        
        # Per-pintle flow rates
        lox_mdot = self.lox_flow / self.n_pintles
        fuel_mdot = self.fuel_flow / self.n_pintles
        
        # Convert areas to m2
        lox_area = lox_area_mm2 * 1e-6
        fuel_area = fuel_area_mm2 * 1e-6
        
        # Calculate velocities
        lox_vel = lox_mdot / (self.lox_density * lox_area)
        fuel_vel = fuel_mdot / (self.fuel_density * fuel_area)
        
        # Dynamic pressures
        lox_q = 0.5 * self.lox_density * lox_vel**2
        fuel_q = 0.5 * self.fuel_density * fuel_vel**2
        
        # Pressure drops
        lox_dp = K_lox * lox_q
        fuel_dp = K_fuel * fuel_q
        
        # Adjust fuel velocity for swirl
        fuel_swirl_rad = np.deg2rad(fuel_swirl_deg)
        fuel_vel_axial = fuel_vel * np.cos(fuel_swirl_rad)
        
        # TMR (momentum flux ratio)
        lox_momentum_flux = self.lox_density * lox_vel**2
        fuel_momentum_flux = self.fuel_density * fuel_vel_axial**2
        TMR = lox_momentum_flux / fuel_momentum_flux if fuel_momentum_flux > 0 else 0
        
        # Relative velocity (for atomization)
        v_rel = abs(lox_vel - fuel_vel_axial)
        
        # Pressure margins
        lox_margin = (self.lox_dp_max - lox_dp) / self.lox_dp_max * 100
        fuel_margin = (self.fuel_dp_max - fuel_dp) / self.fuel_dp_max * 100
        
        # Feasibility
        lox_feasible = lox_dp <= self.lox_dp_max
        fuel_feasible = fuel_dp <= self.fuel_dp_max
        
        # TMR quality score (optimal = 3-8)
        if 3 <= TMR <= 8:
            TMR_score = 1.0
        elif 2 <= TMR < 3:
            TMR_score = 0.8
        elif 8 < TMR <= 12:
            TMR_score = 0.8
        elif 1.5 <= TMR < 2:
            TMR_score = 0.5
        else:
            TMR_score = 0.2
        
        # Atomization quality (based on relative velocity)
        # Higher v_rel → better atomization
        # Normalize to 0-1 scale (optimal around 40-60 m/s)
        if 40 <= v_rel <= 80:
            atom_score = 1.0
        elif 20 <= v_rel < 40:
            atom_score = 0.7
        elif v_rel >= 80:
            atom_score = 0.9
        else:
            atom_score = 0.4
        
        # Overall score (weighted)
        # Weights: pressure margin (0.4), TMR (0.4), atomization (0.2)
        if lox_feasible and fuel_feasible:
            score = 0.3 * (lox_margin/100) + 0.3 * (fuel_margin/100) + 0.4 * TMR_score + 0.2 * atom_score
        else:
            score = 0  # Infeasible
        
        return {
            'lox_area_mm2': lox_area_mm2,
            'fuel_area_mm2': fuel_area_mm2,
            'fuel_swirl_deg': fuel_swirl_deg,
            'lox_velocity_m_s': lox_vel,
            'fuel_velocity_m_s': fuel_vel,
            'fuel_velocity_axial_m_s': fuel_vel_axial,
            'lox_dp_MPa': lox_dp / 1e6,
            'fuel_dp_MPa': fuel_dp / 1e6,
            'lox_margin_pct': lox_margin,
            'fuel_margin_pct': fuel_margin,
            'TMR': TMR,
            'v_relative_m_s': v_rel,
            'lox_feasible': lox_feasible,
            'fuel_feasible': fuel_feasible,
            'TMR_score': TMR_score,
            'atomization_score': atom_score,
            'overall_score': score
        }
    
    def explore_design_space(self):
        """Generate multiple design variants and rank them"""
        
        designs = []
        
        # Parameter ranges to explore
        lox_areas = [35, 40, 45, 49.5, 50, 55]  # mm2
        fuel_areas = [70, 80, 85, 90, 95, 100, 110]  # mm2
        swirl_angles = [0, 15, 30, 45]  # degrees
        
        # Generate all combinations
        for lox_a in lox_areas:
            for fuel_a in fuel_areas:
                for swirl in swirl_angles:
                    result = self.calculate_performance(lox_a, fuel_a, swirl)
                    designs.append(result)
        
        # Convert to DataFrame
        df = pd.DataFrame(designs)
        
        # Sort by overall score
        df = df.sort_values('overall_score', ascending=False)
        
        return df
    
    def summarize_top_designs(self, df, n=10):
        """Print summary of top N designs"""
        
        print("=" * 100)
        print(f"TOP {n} INJECTOR DESIGNS (Ranked by Overall Score)")
        print("=" * 100)
        print()
        
        for i, row in df.head(n).iterrows():
            rank = list(df.index).index(i) + 1
            
            print(f"RANK #{rank} | Score: {row['overall_score']:.3f}")
            print("-" * 100)
            print(f"  LOX area:       {row['lox_area_mm2']:.1f} mm2 (baseline: {self.baseline_lox_area})")
            print(f"  Fuel area:      {row['fuel_area_mm2']:.1f} mm2 (baseline: {self.baseline_fuel_area})")
            print(f"  Fuel swirl:     {row['fuel_swirl_deg']:.0f}°")
            print()
            print(f"  LOX velocity:   {row['lox_velocity_m_s']:.1f} m/s")
            print(f"  Fuel velocity:  {row['fuel_velocity_m_s']:.1f} m/s (axial: {row['fuel_velocity_axial_m_s']:.1f} m/s)")
            print(f"  Relative vel:   {row['v_relative_m_s']:.1f} m/s")
            print()
            print(f"  LOX dP:         {row['lox_dp_MPa']:.3f} MPa (margin: {row['lox_margin_pct']:+.0f}%)")
            print(f"  Fuel dP:        {row['fuel_dp_MPa']:.3f} MPa (margin: {row['fuel_margin_pct']:+.0f}%)")
            print()
            print(f"  TMR:            {row['TMR']:.2f} {'[OPTIMAL]' if 3<=row['TMR']<=8 else '[SUBOPTIMAL]'}")
            print(f"  TMR score:      {row['TMR_score']:.2f}")
            print(f"  Atomization:    {row['atomization_score']:.2f}")
            print()

def main():
    print("=" * 100)
    print("VORTEX INJECTOR DESIGN SPACE EXPLORATION")
    print("=" * 100)
    print()
    
    # Create design space explorer
    explorer = InjectorDesignSpace()
    
    print("CONSTRAINTS:")
    print(f"  LOX max dP:  {explorer.lox_dp_max/1e6:.2f} MPa")
    print(f"  Fuel max dP: {explorer.fuel_dp_max/1e6:.2f} MPa")
    print()
    print("OBJECTIVES:")
    print("  1. Feasible pressure drops (LOX and Fuel)")
    print("  2. TMR in range 3-8 (optimal impingement)")
    print("  3. Good relative velocity (atomization)")
    print("  4. Maximize pressure margins (safety)")
    print()
    print("Exploring design space...")
    print()
    
    # Generate all design variants
    df = explorer.explore_design_space()
    
    print(f"Generated {len(df)} design variants")
    print(f"Feasible designs: {len(df[df['overall_score'] > 0])}")
    print()
    
    # Show top 10 designs
    explorer.summarize_top_designs(df, n=10)
    
    # Export full results
    df.to_csv("D:/realFluidAllMach/work/injector_design_space.csv", index=False)
    print("=" * 100)
    print("Full results exported to: D:/realFluidAllMach/work/injector_design_space.csv")
    print()
    
    # Extract top 3 for detailed analysis
    top3 = df.head(3).to_dict('records')
    
    output = {
        "design_space_size": len(df),
        "feasible_designs": int(len(df[df['overall_score'] > 0])),
        "baseline": {
            "lox_area_mm2": explorer.baseline_lox_area,
            "fuel_area_mm2": explorer.baseline_fuel_area,
            "performance": explorer.calculate_performance(
                explorer.baseline_lox_area, 
                explorer.baseline_fuel_area
            )
        },
        "top_3_designs": top3
    }
    
    with open("D:/realFluidAllMach/work/injector_top_designs.json", 'w') as f:
        json.dump(output, f, indent=2)
    
    print("Top 3 designs exported to: D:/realFluidAllMach/work/injector_top_designs.json")
    print()
    
    print("=" * 100)
    print("RECOMMENDATION")
    print("=" * 100)
    print()
    
    best = df.iloc[0]
    print(f"Best design: LOX={best['lox_area_mm2']:.0f} mm2, Fuel={best['fuel_area_mm2']:.0f} mm2, Swirl={best['fuel_swirl_deg']:.0f}°")
    print(f"  TMR: {best['TMR']:.2f} (optimal range)")
    print(f"  Fuel dP: {best['fuel_dp_MPa']:.3f} MPa (margin: {best['fuel_margin_pct']:+.0f}%)")
    print(f"  LOX dP: {best['lox_dp_MPa']:.3f} MPa (margin: {best['lox_margin_pct']:+.0f}%)")
    print()
    print("NEXT STEPS:")
    print("  1. Generate STL geometry for top 3 designs")
    print("  2. Run CFD cold flow validation")
    print("  3. Compare performance metrics")
    print("  4. Select final design for hot-fire testing")
    print()

if __name__ == "__main__":
    main()
