#!/usr/bin/env python3
"""Assert separated near-exit channels remain background in A44 diagnostic field.

This is a regression for the *observed failure to fill*, not acceptance of the
injector design or ignition. It depends on the audited mesh and A44 field.
"""
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from tetra_field_to_vtu import read_mesh,read_field
from audit_tetra_preexit_components import components


def main():
    mesh=ROOT/'work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet'
    field=ROOT/'work/a406_rootcause/aroute_a44_final.rfamtetf'
    points,cells,_=read_mesh(mesh)
    _,f,y,_=read_field(field)
    for z,expected in [(276.8,71),(277.68,35),(280.9,5),(283.3,2)]:
        origin,labels,area,region_area,_=components(points,cells,z/1000)
        assert len(region_area)==expected,(z,len(region_area))
        outside=region_area.argmax()
        independent=labels!=outside
        assert independent.any(),z
        oxy=y[origin[independent],0].astype(float)
        assert np.max(np.abs(oxy-.784))<1e-6,(z,oxy.min(),oxy.max())
        outside=labels==region_area.argmax()
        external_y=y[origin[outside],0].astype(float)
        assert np.max(np.abs(external_y-.784))<1e-6,(z,external_y.min(),external_y.max())
        print('PREEXIT_BACKGROUND z_mm',z,'separate_regions',expected-1,
              'internal_Y_range',(oxy.min(),oxy.max()),
              'external_Y_range',(external_y.min(),external_y.max()))
    print('PREEXIT_FIELD_REGRESSION_PASS (unconverged diagnostic)')


if __name__=='__main__':main()
