#!/usr/bin/env python3
"""Conditional homogeneous O2/C3H8 ignition delay (NOT injector ignition).

Ideal-gas constant-volume Cantera chemistry at a specified uniform T,p,O/F.
Cryogenic LOX, mixing, torch spark deposition, heat losses and flame transport
are outside this model; never infer hardware ignition from this table.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import cantera as ct

MW_O2=31.9988
MW_C3H8=44.0956
STOICH_MASS_RATIO=5.*MW_O2/MW_C3H8


def ignition_delay(gas,temperature,pressure,of_ratio,max_time,rise):
    gas.TPX=temperature,pressure,{'O2':of_ratio/MW_O2,'C3H8':1./MW_C3H8}
    reactor=ct.IdealGasReactor(gas,energy='on')
    network=ct.ReactorNet([reactor])
    threshold=temperature+rise
    prev_t=0.
    prev_T=temperature
    while network.time<max_time:
        t=network.step()
        new_T=reactor.T
        if new_T>=threshold:
            delay=prev_t+(threshold-prev_T)*(t-prev_t)/(new_T-prev_T)
            return float(delay),float(new_T)
        if t>=max_time:break
        prev_t,prev_T=t,new_T
    return None,float(reactor.T)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mechanism',required=True)
    parser.add_argument('--out',required=True)
    parser.add_argument('--pressure-MPa',type=float,default=4.)
    parser.add_argument('--temperatures',type=float,nargs='+',default=[800.,900.,1000.,1100.,1200.])
    parser.add_argument('--of-ratios',type=float,nargs='+',default=[STOICH_MASS_RATIO,2.667])
    parser.add_argument('--max-time-ms',type=float,default=20.)
    parser.add_argument('--rise-K',type=float,default=400.)
    args=parser.parse_args()
    path=Path(args.mechanism)
    source=hashlib.sha256(path.read_bytes()).hexdigest()
    result={'scope':'conditional 0D uniform ideal-gas chemistry; NOT injector or torch ignition',
            'model':'adiabatic constant-volume Cantera IdealGasReactor',
            'mechanism':str(path),'mechanism_sha256':source,'cantera_version':ct.__version__,
            'pressure_Pa':args.pressure_MPa*1e6,'temperature_rise_K':args.rise_K,
            'observation_window_s':args.max_time_ms*1e-3,
            'stoichiometric_of_mass_ratio':STOICH_MASS_RATIO,
            'missing_physics':['LOX liquid and phase transfer','spatial mixture fraction and turbulent mixing',
                               'spark energy and hot torch boundary','wall heat losses and unsteady flow'],
            'cases':[]}
    for of in args.of_ratios:
        if of<=0:raise ValueError('positive mass O/F required')
        for temperature in args.temperatures:
            if temperature<600.:raise ValueError('this gas-only screen is not for cryogenic states')
            gas=ct.Solution(str(path))
            delay,last=ignition_delay(gas,temperature,args.pressure_MPa*1e6,of,
                                      args.max_time_ms*1e-3,args.rise_K)
            case={'of_mass_ratio':of,'equivalence_ratio':STOICH_MASS_RATIO/of,
                  'initial_T_K':temperature,'ignition_delay_s':delay,
                  'final_observed_T_K':last,'criterion':'T reaches initial_T + rise_K'}
            result['cases'].append(case)
            print(f'O/F={of:.4f} phi={STOICH_MASS_RATIO/of:.3f} T0={temperature:.0f} K '
                  f'IDT={delay if delay is not None else "not reached"} s',flush=True)
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    Path(args.out).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('CONDITIONAL_0D_IGNITION_WRITTEN',args.out)


if __name__=='__main__':main()
