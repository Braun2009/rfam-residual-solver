#!/usr/bin/env python3
"""0D chemistry regression only; not a rocket injector ignition test."""
import sys
from pathlib import Path
import warnings

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from audit_conditional_ignition import STOICH_MASS_RATIO,ignition_delay
import cantera as ct


def main():
    assert abs(STOICH_MASS_RATIO-3.6283)<.001
    assert 1.35<STOICH_MASS_RATIO/2.667<1.37
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',UserWarning)
        gas=ct.Solution(str(ROOT/'tests/reference/sandiego.yaml'))
    assert gas.n_species==59 and gas.n_reactions==271
    for T,reference in [(1000.,4.515e-4),(1200.,3.5e-5)]:
        delay,_=ignition_delay(gas,T,4e6,STOICH_MASS_RATIO,.005,400.)
        assert delay is not None and abs(delay/reference-1.)<.05,(T,delay)
    print('CONDITIONAL_0D_CHEMISTRY_PASS (NOT hardware ignition)')


if __name__=='__main__':main()
