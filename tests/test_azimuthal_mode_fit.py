#!/usr/bin/env python3
"""A nonuniform angular mesh must not create a mode in a constant field."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from coldflow_characteristics import fit_azimuthal_modes


def main():
    rng=np.random.default_rng(14)
    theta=np.r_[rng.uniform(-.8,.6,700),rng.uniform(1.5,2.9,120)]
    weight=rng.uniform(.1,4.,len(theta))
    baseline=fit_azimuthal_modes(theta,weight,np.full(len(theta),.784))
    assert baseline['m1']<1e-12 and baseline['m3']<1e-12,baseline
    field=.784+.03*np.cos(theta+.4)+.12*np.cos(3*theta-1.2)
    modes=fit_azimuthal_modes(theta,weight,field)
    assert abs(modes['mean']-.784)<1e-12,modes
    assert abs(modes['m1']-.03)<1e-12,modes
    assert abs(modes['m3']-.12)<1e-12,modes
    assert abs(modes['phase_deg']-np.degrees(1.2)/3)<1e-12,modes
    try:fit_azimuthal_modes(np.zeros(20),np.ones(20),np.ones(20))
    except ValueError:pass
    else:raise AssertionError('rank-deficient angular sampling was accepted')
    print('AZIMUTHAL_MODE_FIT_PASS')


if __name__=='__main__':main()
