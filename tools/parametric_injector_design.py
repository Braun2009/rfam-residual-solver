#!/usr/bin/env python3
"""
Parametric Vortex Injector Design Generator
Creates STL geometry from design parameters
Clean-sheet design with single unified manifolds
"""

import numpy as np
import json

class InjectorDesign:
    """Parametric injector with integrated manifolds"""
    
    def __init__(self):
        # Design parameters (all in mm)
        
        # Chamber geometry (YOU NEED TO PROVIDE)
        self.chamber_id = 80  # mm - PLACEHOLDER
        self.chamber_length = 300  # mm - PLACEHOLDER
        self.throat_diameter = 40  # mm - PLACEHOLDER
        
        # Manifold geometry
        self.lox_manifold_id = 40
        self.lox_manifold_od = 80
        self.lox_manifold_height = 20
        
        self.fuel_manifold_id = 85
        self.fuel_manifold_od = 120
        self.fuel_manifold_height = 15
        
        # Pintle configuration
        self.n_pintles = 3
        self.pintle_diameter = 10
        self.pintle_length = 50
        self.pintle_spacing_radius = 30  # from center
        
        # LOX orifices (at pintle tip)
        self.lox_holes_per_pintle = 4
        self.lox_hole_diameter = 2.5
        self.lox_recess = 2  # mm, LOX recessed into pintle
        
        # Fuel swirl slots (around pintle)
        self.fuel_slots_per_pintle = 8
        self.fuel_slot_width = 3.0
        self.fuel_slot_length = 8.0
        self.fuel_swirl_angle = 30  # degrees
        
        # Torch igniter
        self.torch_enable = True
        self.torch_lox_diameter = 3.0
        self.torch_fuel_od = 5.0
        
        # Flow requirements
        self.lox_flow_kg_s = 2.8
        self.fuel_flow_kg_s = 1.05
        self.lox_density = 909.3  # kg/m3
        self.fuel_density = 65.5  # kg/m3
        
        # Pressure constraints
        self.fuel_dp_max = 0.5  # MPa - CRITICAL
        
    def calculate_areas(self):
        """Calculate required flow areas"""
        
        # LOX area per pintle
        lox_area_per_hole = np.pi * (self.lox_hole_diameter/2)**2
        self.lox_area_per_pintle = lox_area_per_hole * self.lox_holes_per_pintle
        self.lox_area_total = self.lox_area_per_pintle * self.n_pintles
        
        # Fuel area per pintle
        fuel_area_per_slot = self.fuel_slot_width * self.fuel_slot_length
        self.fuel_area_per_pintle = fuel_area_per_slot * self.fuel_slots_per_pintle
        self.fuel_area_total = self.fuel_area_per_pintle * self.n_pintles
        
        print(f"LOX area per pintle:  {self.lox_area_per_pintle:.1f} mm2 ({self.lox_holes_per_pintle} holes)")
        print(f"LOX total area:       {self.lox_area_total:.1f} mm2")
        print(f"Fuel area per pintle: {self.fuel_area_per_pintle:.1f} mm2 ({self.fuel_slots_per_pintle} slots)")
        print(f"Fuel total area:      {self.fuel_area_total:.1f} mm2")
        
    def calculate_pressure_drops(self):
        """Estimate pressure drops using K factors"""
        
        # LOX side
        lox_mdot_per_pintle = self.lox_flow_kg_s / self.n_pintles
        lox_area_m2 = self.lox_area_per_pintle * 1e-6
        lox_velocity = lox_mdot_per_pintle / (self.lox_density * lox_area_m2)
        lox_dynamic_p = 0.5 * self.lox_density * lox_velocity**2
        
        # K factors: manifold (0.1) + feed tube (0.3) + orifice rounded (1.5)
        lox_K_total = 0.1 + 0.3 + 1.5
        lox_dp = lox_K_total * lox_dynamic_p
        
        # Fuel side
        fuel_mdot_per_pintle = self.fuel_flow_kg_s / self.n_pintles
        fuel_area_m2 = self.fuel_area_per_pintle * 1e-6
        fuel_velocity = fuel_mdot_per_pintle / (self.fuel_density * fuel_area_m2)
        fuel_dynamic_p = 0.5 * self.fuel_density * fuel_velocity**2
        
        # K factors: manifold (0.05) + feed tube (0.15) + swirl slots (2.0)
        fuel_K_total = 0.05 + 0.15 + 2.0
        fuel_dp = fuel_K_total * fuel_dynamic_p
        
        print()
        print("PRESSURE DROP ESTIMATES:")
        print(f"LOX velocity:  {lox_velocity:.1f} m/s")
        print(f"LOX dP:        {lox_dp/1e6:.3f} MPa (K={lox_K_total})")
        print(f"Fuel velocity: {fuel_velocity:.1f} m/s")
        print(f"Fuel dP:       {fuel_dp/1e6:.3f} MPa (K={fuel_K_total})")
        print()
        
        if fuel_dp > self.fuel_dp_max * 1e6:
            print(f"WARNING:  WARNING: Fuel dP exceeds constraint ({self.fuel_dp_max} MPa)")
            print(f"   Need to increase fuel area or reduce K factors")
            
            # Calculate required area
            required_area = fuel_area_m2 * (fuel_dp / (self.fuel_dp_max * 1e6))
            print(f"   Required fuel area: {required_area*1e6:.0f} mm2 (current: {self.fuel_area_per_pintle:.0f} mm2)")
        else:
            margin = (self.fuel_dp_max * 1e6 - fuel_dp) / (self.fuel_dp_max * 1e6) * 100
            print(f"OK: Fuel dP within constraint (+{margin:.0f}% margin)")
        
        return lox_dp, fuel_dp
    
    def generate_geometry(self):
        """Generate 3D geometry (simplified - actual CAD would be more detailed)"""
        
        print()
        print("GEOMETRY GENERATION:")
        print("This is a parametric design - actual STL generation requires CAD software")
        print("Recommended workflow:")
        print("  1. Use OpenSCAD with these parameters")
        print("  2. Or FreeCAD Python scripting")
        print("  3. Or write custom mesh generation code")
        print()
        print("Key dimensions for CAD:")
        print(f"  LOX manifold: ID={self.lox_manifold_id} OD={self.lox_manifold_od} H={self.lox_manifold_height}")
        print(f"  Fuel manifold: ID={self.fuel_manifold_id} OD={self.fuel_manifold_od} H={self.fuel_manifold_height}")
        print(f"  Pintles: {self.n_pintles} × D={self.pintle_diameter} at R={self.pintle_spacing_radius}")
        print(f"  LOX holes: {self.lox_holes_per_pintle} × D={self.lox_hole_diameter}")
        print(f"  Fuel slots: {self.fuel_slots_per_pintle} × {self.fuel_slot_width}×{self.fuel_slot_length}")
        
    def export_parameters(self, filename):
        """Export design parameters to JSON"""
        
        params = {
            "design_type": "parametric_pintle_injector",
            "version": "1.0",
            "chamber": {
                "id_mm": self.chamber_id,
                "length_mm": self.chamber_length,
                "throat_diameter_mm": self.throat_diameter
            },
            "manifolds": {
                "lox": {
                    "id_mm": self.lox_manifold_id,
                    "od_mm": self.lox_manifold_od,
                    "height_mm": self.lox_manifold_height
                },
                "fuel": {
                    "id_mm": self.fuel_manifold_id,
                    "od_mm": self.fuel_manifold_od,
                    "height_mm": self.fuel_manifold_height
                }
            },
            "pintles": {
                "count": self.n_pintles,
                "diameter_mm": self.pintle_diameter,
                "length_mm": self.pintle_length,
                "spacing_radius_mm": self.pintle_spacing_radius
            },
            "lox_injection": {
                "holes_per_pintle": self.lox_holes_per_pintle,
                "hole_diameter_mm": self.lox_hole_diameter,
                "recess_mm": self.lox_recess,
                "area_per_pintle_mm2": float(self.lox_area_per_pintle),
                "total_area_mm2": float(self.lox_area_total)
            },
            "fuel_injection": {
                "slots_per_pintle": self.fuel_slots_per_pintle,
                "slot_width_mm": self.fuel_slot_width,
                "slot_length_mm": self.fuel_slot_length,
                "swirl_angle_deg": self.fuel_swirl_angle,
                "area_per_pintle_mm2": float(self.fuel_area_per_pintle),
                "total_area_mm2": float(self.fuel_area_total)
            },
            "torch": {
                "enabled": self.torch_enable,
                "lox_diameter_mm": self.torch_lox_diameter,
                "fuel_od_mm": self.torch_fuel_od
            },
            "flow_requirements": {
                "lox_kg_s": self.lox_flow_kg_s,
                "fuel_kg_s": self.fuel_flow_kg_s
            },
            "constraints": {
                "fuel_dp_max_MPa": self.fuel_dp_max
            }
        }
        
        with open(filename, 'w') as f:
            json.dump(params, f, indent=2)
        
        print(f"Parameters exported to: {filename}")

def main():
    print("=" * 70)
    print("VORTEX INJECTOR - PARAMETRIC DESIGN TOOL")
    print("=" * 70)
    print()
    
    # Create design
    design = InjectorDesign()
    
    print("DESIGN PARAMETERS:")
    print(f"Flow rates: {design.lox_flow_kg_s} kg/s LOX, {design.fuel_flow_kg_s} kg/s fuel")
    print(f"Pintles: {design.n_pintles}")
    print(f"Fuel dP constraint: {design.fuel_dp_max} MPa (from regen cooling)")
    print()
    
    # Calculate areas
    design.calculate_areas()
    
    # Estimate pressure drops
    lox_dp, fuel_dp = design.calculate_pressure_drops()
    
    # Generate geometry (conceptual)
    design.generate_geometry()
    
    # Export parameters
    design.export_parameters("D:/realFluidAllMach/work/injector_parametric_design.json")
    
    print()
    print("=" * 70)
    print("NEXT STEPS:")
    print("=" * 70)
    print()
    print("1. Review parameters and adjust if needed")
    print("2. Provide actual chamber dimensions (ID, length, throat)")
    print("3. Choose CAD tool:")
    print("   - OpenSCAD (text-based, parametric)")
    print("   - FreeCAD (Python API)")
    print("   - Custom mesh generation")
    print("4. Generate STL geometry")
    print("5. Run CFD validation")
    print()

if __name__ == "__main__":
    main()
