#!/usr/bin/env python3
"""The delivery mesh's 64 inlet branches remain disjoint before merging."""
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from audit_injector_passage_connectivity import selection
from tetra_field_to_vtu import read_mesh


def main():
    mesh=ROOT/'work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet'
    tags,details=selection(mesh)
    expected={'lox':(3,91289),'fuel':(60,204832),'torch':(1,8710)}
    assert len(tags)==1761002
    assert int(np.count_nonzero(tags))==304831
    for family,(branches,cells) in expected.items():
        info=details[family]
        assert info['components']==branches and info['cells']==cells,(family,info)
        assert info['highest_selected_vertex_z_mm']<info['cutoff_z_mm']
    points,cells,_=read_mesh(mesh)
    board=.2709527893066406
    exterior=np.flatnonzero((tags==0)&(points[cells][:,:,2].min(axis=1)>.2858028))
    assert len(exterior)>900000,len(exterior)
    print('INJECTOR_PASSAGE_MASK_PASS separated_branches=64 selected_cells=304831')


if __name__=='__main__':main()
