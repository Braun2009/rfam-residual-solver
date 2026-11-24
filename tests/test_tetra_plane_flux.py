#!/usr/bin/env python3
"""Analytic tetra/plane cut and uniform axial-flux regression."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from audit_tetra_plane_flux import cut_area_by_cell


def main():
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    cells=np.array([[0,1,2,3]],dtype=np.int32)
    ids,areas=cut_area_by_cell(points,cells,.5)
    assert len(ids)==1 and ids[0]==0
    assert abs(areas[0]-.125)<1e-12,areas
    rho,uz,y=2.,3.,.8
    assert abs(float(np.sum(rho*uz*areas))-.75)<1e-12
    assert abs(float(np.sum(rho*uz*areas*y))-.6)<1e-12
    try:cut_area_by_cell(points,cells,2.)
    except ValueError:pass
    else:raise AssertionError('cut outside mesh accepted')
    print('TETRA_PLANE_FLUX_PASS')


if __name__=='__main__':main()
