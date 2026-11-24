#!/usr/bin/env python3
"""Independent HEOS binary equilibrium isentropes; no forced mixture phase."""
import json,math
from pathlib import Path
import numpy as np
from scipy.optimize import root
from CoolProp import AbstractState,CoolProp as CP
import CoolProp


def main():
    rows=[]
    for y,p,t in [(.9,4e6,170),(.99,4e6,150),(.73,4e6,250),(.1,4e6,250)]:
        s=AbstractState('HEOS','Oxygen&Propane');s.set_mass_fractions([y,1-y])
        s.update(CP.PT_INPUTS,p,t);rho=s.rhomass();entropy=s.smass()
        molar_mass=[CP.PropsSI('M','Oxygen'),CP.PropsSI('M','Propane')]
        xliq=s.mole_fractions_liquid();xvap=s.mole_fractions_vapor();beta_mole=s.Q()
        phase=[]
        for x,phase_id,density in [(xliq,CP.iphase_liquid,s.saturated_liquid_keyed_output(CP.iDmass)),
                                   (xvap,CP.iphase_gas,s.saturated_vapor_keyed_output(CP.iDmass))]:
            endpoint=AbstractState('HEOS','Oxygen&Propane');endpoint.set_mole_fractions(x)
            # Fixed-density endpoint evaluation verifies the flashed endpoints;
            # this does not select or impose the bulk equilibrium phase.
            endpoint.specify_phase(phase_id);endpoint.update(CP.DmassT_INPUTS,density,t);phase.append(endpoint)
        fugacity_residual=max(abs(math.log(phase[0].fugacity(i)/phase[1].fugacity(i))) for i in range(2))
        if fugacity_residual>1e-7:raise RuntimeError('phase fugacities do not close')
        results=[]
        for eps in [1e-3,5e-4,2.5e-4,1.25e-4]:
            pp=[];states=[]
            for sign in [-1,1]:
                target=rho*(1+sign*eps)
                def residual(z):
                    s.update(CP.PT_INPUTS,p*math.exp(z[0]),t*math.exp(z[1]));return [(s.rhomass()-target)/rho,(s.smass()-entropy)/max(abs(entropy),1)]
                sol=root(residual,[0,0],tol=1e-10)
                err=residual(sol.x)
                if max(abs(x) for x in err)>1e-8:raise RuntimeError(str(sol))
                pp.append(s.p());states.append({'rho':s.rhomass(),'s_error_J_kgK':s.smass()-entropy,'p':s.p(),'T':s.T()})
            a2=(pp[1]-pp[0])/(2*rho*eps)
            results.append({'relative_density_perturbation':eps,'a_m_s':math.sqrt(a2),'states':states})
        rows.append({'Y_O2':y,'p':p,'T':t,'liquid_mole_fractions':xliq,'vapor_mole_fractions':xvap,
                     'molar_vapor_fraction':beta_mole,'max_log_fugacity_residual':fugacity_residual,'results':results})
        print('isentrope',y,p,t,[r['a_m_s'] for r in results],flush=True)
    out=Path(__file__).resolve().parents[1]/'work/coldflow_B_verified/thermo_audit/isentrope_reference.json'
    out.write_text(json.dumps({'coolprop_version':CoolProp.__version__,'coolprop_git':CoolProp.__gitrevision__,'rows':rows},indent=2))

if __name__=='__main__':main()
