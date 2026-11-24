#!/usr/bin/env python3
"""Batch Fluent Fault-tolerant Meshing driver for the split Vortex segments.

Reuses the validated object-TUI route that produced strict_tui_prism_final:
  import zoned ASCII STL (mm) -> describe -> surface -> regions -> boundary
  types -> TUI tet with size-field -> scoped prisms (inner only) -> auto-fill
  -> check -> write .msh.h5

Segments:
  inner : chamber+nozzle x[0,406.5]; boundaries mass-flow-inlet x2, wall x2,
          pressure-outlet (nozzle exit); 12 prism layers, 0.05 mm first,
          rate 1.25 on walls.
  plume : x[406.5,1406.5]; boundaries mass-flow-inlet (nozzle exit profile
          inlet), wall (aft stand-in), pressure-outlet (farfield); pure tet.
"""
import argparse
import json
import os
import sys
import time
import traceback

os.environ.update(AWP_ROOT261='/ansys_inc/v261',
                  ANSYSLMD_LICENSE_FILE='1055@localhost',
                  ANSYSCL_FNP_PATH='1055@localhost')

import ansys.fluent.core as pyfluent
from ansys.fluent.core import FluentMode, UIMode

SEGMENTS = {
    "inner": {
        "stl": "/root/vortex_fluent/gmsh_split_segments/inner_surf/vortex_inner_406p5_zoned.stl",
        "boundary_types": {
            "inlet_torch": "mass-flow-inlet",
            "inlet_spray": "mass-flow-inlet",
            "wall_plate": "wall",
            "wall_chamber": "wall",
            "outlet_nozzle": "pressure-outlet",
        },
        "prism": {"name": "wall_prism_005", "first": 0.05, "layers": 12, "rate": 1.25},
        "out": "vortex_inner_406p5_tui_prism.msh.h5",
    },
    "plume": {
        "stl": "/root/vortex_fluent/gmsh_split_segments/plume_surf2/vortex_plume_406p5_zoned.stl",
        "boundary_types": {
            "inlet_nozzle_exit": "mass-flow-inlet",
            "wall_nozzle_aft": "wall",
            "farfield": "pressure-outlet",
        },
        "prism": None,
        "out": "vortex_plume_406p5_tui.msh.h5",
    },
}


def mark(label):
    print("RFAM_STAGE", label, time.strftime("%H:%M:%S"), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", choices=tuple(SEGMENTS), required=True)
    parser.add_argument("--processors", type=int, default=8)
    args = parser.parse_args()
    cfg = SEGMENTS[args.segment]
    base = "/root/vortex_fluent/gmsh_split_segments"
    out_dir = base + "/" + args.segment + "_mesh"
    os.makedirs(out_dir, exist_ok=True)
    out_mesh = out_dir + "/" + cfg["out"]
    meta = {
        "segment": args.segment,
        "stl": cfg["stl"],
        "boundary_types": cfg["boundary_types"],
        "prism": cfg["prism"],
        "out_mesh": out_mesh,
    }
    with open(out_dir + "/mesh_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
        f.write("\n")
    s = None
    try:
        s = pyfluent.launch_fluent(precision="double", processor_count=args.processors,
                                   mode=FluentMode.MESHING, ui_mode=UIMode.NO_GUI,
                                   start_timeout=1200)
        w = s.workflow
        w.InitializeWorkflow(WorkflowType="Fault-tolerant Meshing")
        mark("WF_INIT_OK")
        imp = w.TaskObject["Import CAD and Part Management"]
        imp.Arguments = {
            "FMDFileName": cfg["stl"], "FileLengthUnit": "mm", "LengthUnit": "mm",
            "CreateObjectPer": "One per part", "Route": "Native", "OneZonePer": "body",
            "PartPerBody": False, "IgnoreSolidNames": False, "EdgeExtraction": "no",
        }
        imp.Execute()
        mark("IMPORT_OK")
        d = w.TaskObject["Describe Geometry and Flow"]
        d.Arguments = {
            "FlowType": "Internal flow through object", "GeometryOptions": True,
            "AddEnclosure": "No", "CloseCaps": "No", "LocalRefinementRegions": "No",
            "DescribeGeometryAndFlowOptions": {
                "IdentifyRegions": "Yes", "CloseLeakages": "No", "PorousRegions": "No",
                "ZeroThickness": "No", "ExtractEdgeFeatures": "No", "MovingObjects": "No",
                "EnableOverset": "No", "EnablePrimeWrapper": "No",
            },
        }
        d.Execute()
        mark("DESCRIBE_OK")
        sm = w.TaskObject["Generate the Surface Mesh"]
        sm.Arguments = {
            "SurfaceQuality": 0.8, "SaveSurfaceMesh": True, "AdvancedOptions": True,
            "SeparateSurface": "No", "ProjectOnGeometry": "No",
            "AutoAssignZoneTypes": "No", "GlobalMin": 0.2,
        }
        sm.Execute()
        mark("SURFACE_OK")
        w.TaskObject["Compute Regions"].Execute()
        mark("REGIONS_OK")
        ub = w.TaskObject["Update Boundaries"]
        zone_list = sorted(cfg["boundary_types"])
        ub.Arguments = {
            "SelectionType": "zone",
            "BoundaryZoneList": zone_list,
            "BoundaryZoneTypeList": [cfg["boundary_types"][z] for z in zone_list],
        }
        ub.Execute()
        mark("BOUNDARIES_OK")
        obj_name = os.path.basename(cfg["stl"]).rsplit(".", 1)[0]
        ur = w.TaskObject["Update Region Settings"]
        ur.Arguments = {
            "MainFluidRegion": obj_name,
            "AllRegionNameList": [obj_name],
            "AllRegionMeshMethodList": ["tet mesh"],
            "AllRegionTypeList": ["fluid"],
            "AllRegionVolumeFillList": ["tet"],
            "AllRegionLeakageSizeList": ["none"],
            "AllRegionOversetComponenList": ["no"],
        }
        ur.Execute()
        mark("RSET_OK")
        # Boundary types must be set AFTER Update Region Settings: the region
        # task otherwise resets the zones back to wall (validated by failure).
        # The TUI loops over zone/type pairs until a blank line exits the loop.
        type_cmds = "/boundary/manage/type\n"
        for zone_name, bc_type in cfg["boundary_types"].items():
            type_cmds += f"{zone_name}\n{bc_type}\n"
        tui = s.execute_tui
        tui(type_cmds + "\n")
        mark("BC_TYPES_TUI_OK")

        tui = s.execute_tui
        tui("/mesh/tet/controls/cell_sizing\nsize-field\n/mesh/tet/controls/use-max-cell-size\nno")
        mark("TET_CONTROLS_OK")
        tui("/mesh/tet/init")
        mark("TET_INIT_OK")
        tui(f"/mesh/auto-mesh\n{obj_name}")
        mark("AUTOMESH_OK")
        # manage/type keeps prompting; the trailing blank line exits the loop.
        tui(f"/mesh/manage/type\n{obj_name}-tet-cells\nfluid\n\n")
        mark("FLUID_TYPE_OK")
        # Explicitly set the volumetric region type so the scoped-prism
        # control can resolve its scope (validated on the final-grid run).
        tui(f"/objects/volumetric-regions/change-type\n{obj_name}\n*\n\nfluid\n")
        mark("REGION_TYPE_OK")
        if cfg["prism"]:
            p = cfg["prism"]
            tui(f"/mesh/scoped-prisms/create\n{p['name']}\nuniform\n{p['first']}\n"
                f"{p['layers']}\n{p['rate']}\n{obj_name}\nfluid-regions\nonly-walls\n")
            mark("PRISM_CONTROL_OK")
        tui(f"/objects/volumetric-regions/delete-cells\n{obj_name}\n*\n\n")
        if cfg["prism"]:
            # 'Grow Prisms' answer 'scoped' makes auto-fill grow the layers
            # defined by the scoped-prism control (validated on final-grid run).
            tui(f"/objects/volumetric-regions/auto-fill-volume\n{obj_name}\n*\n\nscoped\n\n\n")
        else:
            tui(f"/objects/volumetric-regions/auto-fill-volume\n{obj_name}\n*\n\n")
        mark("VOLUME_OK")
        tui("/mesh/check-mesh")
        mark("MESH_CHECK_OK")
        tui("/mesh/check-quality")
        mark("QUALITY_OK")
        if cfg["prism"]:
            tui("/mesh/prism/improve/smooth-improve-prism-cells")
            mark("PRISM_IMPROVE_OK")
        tui(f"/file/write-mesh\n{out_mesh}")
        mark("MESH_WRITTEN")
        print("VORTEX_SPLIT_MESH_OK", args.segment, out_mesh, flush=True)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
    finally:
        if s is not None:
            try:
                s.exit()
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    main()
