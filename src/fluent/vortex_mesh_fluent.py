#!/usr/bin/env python3
"""vortex_mesh_fluent.py - Fluent Meshing watertight workflow for
vortex_full360_fluid.step (full-360 chamber + plume + injector).

Steps: import -> local sizing (curvature+proximity, wall/cone boxes) ->
surface mesh -> describe -> update boundaries (physical names from gmsh
export are preserved as labels) -> poly-hexcore volume mesh.
"""
import sys
import ansys.fluent.core as pyfluent

STEP = "/root/vortex_fluent/vortex_full360_fluid.step"

def main():
    meshing = pyfluent.launch_fluent(
        precision="double", processor_count=16, mode="meshing",
        start_timeout=600, ui_mode="no_gui")
    w = meshing.workflow
    w.InitializeWorkflow(WorkflowType="Watertight Geometry")

    imp = w.TaskObject["Import Geometry"]
    imp.Arguments = {
        "FileName": STEP,
        "LengthUnit": "mm",
        "CadImportOptions": {"OneObjectPer": "File"},
    }
    imp.Execute()
    print("IMPORT_OK")

    # ---- local sizing: wall first layer + cone + throat + shock ----
    sizing = w.TaskObject["Add Local Sizing"]
    sizing.Arguments = {
        "AddChild": "yes",
        "BOIControlName": "wall_prox",
        "BOIFaceLabelList": ["wall_chamber", "wall_plate", "wall_plate_back"],
        "BOIGrowthRate": 1.2,
        "BOISize": 0.4,
    }
    sizing.Execute()
    sizing2 = w.TaskObject["Add Local Sizing"]
    sizing2.Arguments = {
        "AddChild": "yes",
        "BOIControlName": "cone_prox",
        "BOIFaceLabelList": ["inlet_lox", "inlet_fuel", "inlet_torch"],
        "BOIGrowthRate": 1.2,
        "BOISize": 0.8,
    }
    sizing2.Execute()
    print("SIZING_OK")

    # ---- surface mesh ----
    sm = w.TaskObject["Generate the Surface Mesh"]
    sm.Arguments = {
        "CFDSurfaceMeshControls": {
            "MaxSize": 8.0,
            "MinSize": 0.4,
            "CurvatureNormalAngle": 18,
            "CellsPerGap": 1,
        },
    }
    sm.Execute()
    print("SURFACE_MESH_OK")

    # ---- describe + update boundaries ----
    w.TaskObject["Describe Geometry"].Arguments = {
        "SetupType": "The geometry consists of only fluid regions with no voids"}
    w.TaskObject["Describe Geometry"].Execute()
    ub = w.TaskObject["Update Boundaries"]
    ub.Arguments = {
        "LabelType": "zone",
    }
    ub.Execute()
    print("BOUNDARIES_OK")

    # ---- volume mesh ----
    vm = w.TaskObject["Generate the Volume Mesh"]
    vm.Arguments = {
        "VolumeFill": "poly-hexcore",
        "VolumeMeshControls": {
            "HexMaxCellLength": 12.0,
            "BufferLayers": 2,
            "PeelLayers": 2,
        },
    }
    vm.Execute()
    print("VOLUME_MESH_OK")

    meshing.tui.mesh.check_mesh()
    meshing.tui.file.write_case("/root/vortex_fluent/vortex360_mesh.cas.h5")
    print("MESH_CASE_WRITTEN")
    meshing.exit()

if __name__ == "__main__":
    main()
