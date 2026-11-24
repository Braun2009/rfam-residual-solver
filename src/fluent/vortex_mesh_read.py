#!/usr/bin/env python3
"""vortex_mesh_read.py - Fluent Meshing Fault-tolerant workflow: import
gmsh surface mesh (.msh), describe, update boundaries, volume fill."""
import sys
import ansys.fluent.core as pyfluent

MSH = "/root/vortex_fluent/vortex_full360_surface.msh"

def main():
    meshing = pyfluent.launch_fluent(
        precision="double", processor_count=16, mode="meshing",
        start_timeout=600, ui_mode="no_gui")
    w = meshing.workflow
    w.InitializeWorkflow(WorkflowType="Fault-tolerant Meshing")
    imp = w.TaskObject["Import CAD and Part Management"]
    print("IMP_SCHEMA:", sorted(imp.Arguments.get_state()))
    print("IMP_OBJ:", imp.Arguments.ObjectSetting.get_state())
    imp.Arguments = {"FileName": MSH, "LengthUnit": "mm"}
    imp.Execute()
    print("IMPORT_OK")
    dg = w.TaskObject["Describe Geometry and Flow"]
    dg.Arguments = {"SetupType": "The geometry consists of only fluid regions with no voids"}
    dg.Execute()
    print("DESCRIBE_OK")
    ub = w.TaskObject["Update Boundaries"]
    ub.Arguments = {"LabelType": "label"}
    ub.Execute()
    print("BOUNDARIES_OK")
    sm = w.TaskObject["Generate the Surface Mesh"]
    sm.Execute()
    print("SURFACE_MESH_OK")
    cr = w.TaskObject["Compute Regions"]
    cr.Execute()
    print("REGIONS_OK")
    ub2 = w.TaskObject["Update Boundaries"]
    ub2.Execute()
    print("BOUNDARIES2_OK")
    vm = w.TaskObject["Generate the Volume Mesh"]
    vm.Arguments = {"VolumeFill": "tet", "VolumeMeshControls": {"MaxCellLength": 8.0}}
    vm.Execute()
    print("VOLUME_FILL_OK")
    tui = meshing.tui
    tui = meshing.tui
    solver = meshing.switch_to_solver()
    solver.mesh.check()
    solver.file.write(file_name="/root/vortex_fluent/vortex360_mesh.cas.h5",
                      file_type="case")
    print("MESH_CASE_WRITTEN")
    meshing.exit()

if __name__ == "__main__":
    main()
