#!/usr/bin/env python3
"""Regression: the audited A-route mesh has 64 matchable upstream passages."""
import json
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from tetra_field_to_vtu import read_mesh
from audit_tetra_passage_flux import passage_cut


def main():
    mesh=ROOT/'work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet'
    audit=ROOT/'work/mesh_rebuild/delivery/reports/mesh_audit.json'
    points,cells,_=read_mesh(mesh)
    patches=json.loads(audit.read_text(encoding='utf-8'))['inlets']
    assert len(patches)==64
    _,labels,areas,region_areas,xy=passage_cut(points,cells,.2705)
    inlet=np.asarray([p['centroid_cad_mm'][:2] for p in patches])/1000
    diff=np.linalg.norm(xy[:,None,:]-inlet[None,:,:],axis=2)
    match=diff.argmin(axis=1)
    assert len(set(match))==64
    assert diff[np.arange(64),match].max()<3e-6
    assert np.bincount(labels,minlength=64).min()>0
    assert abs(areas.sum()*1e6-335.4650942318865)<1e-4
    count={f:sum(patches[j]['design_family']==f for j in match)
           for f in ('fuel','lox','torch')}
    assert count=={'fuel':60,'lox':3,'torch':1},count
    assert np.all(region_areas>0)
    # Above the board the chamber joins the passages. Such a slice cannot
    # silently be reused as the same 64 distinct channels.
    try:passage_cut(points,cells,.2712)
    except ValueError:pass
    else:raise AssertionError('joined chamber was accepted as 64 passages')
    print('TETRA_PASSAGE_CUT_PASS')


if __name__=='__main__':main()
