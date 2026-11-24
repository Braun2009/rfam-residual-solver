#!/usr/bin/env python3
"""Frozen CoolProp 7.2.0 HEOS values and comparison to a tetra backend.

Energy comparison uses each species' zero-density ideal enthalpy at 298.15 K,
matching RFAM's sensible-enthalpy convention. Raw CP h/u are retained as well.
"""
import argparse,csv,json,math
from pathlib import Path
import CoolProp
from CoolProp import AbstractState, CoolProp as CP


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True);parser.add_argument('--rfam')
    args=parser.parse_args();out=Path(args.out);out.mkdir(exist_ok=True,parents=True)
    if CoolProp.__version__!='7.2.0':raise RuntimeError('Requires CoolProp 7.2.0')
    refs=[]
    for fluid in ['Oxygen','Propane']:
        s=AbstractState('HEOS',fluid);s.update(CP.DmassT_INPUTS,1e-10,298.15)
        refs.append(s.hmass_idealgas())
    cases=[('lox_supply',6e6,133.5,1.),('lox_chamber_4MPa',4e6,133.5,1.),
           ('fuel_supply',4.5e6,453.5,0.),('fuel_chamber_4MPa',4e6,453.5,0.)]
    cases += [(f'{fluid}_{p:g}_{t:g}',p,t,y) for fluid,y in [('o2',1.),('fuel',0.)]
              for p in [1e5,5e5,4e6] for t in [150.,300.,400.,453.5,600.]]
    cases += [(f'mix_{p:g}_{t:g}_{y:g}',p,t,y) for p in [5e5,4e6]
              for t in [150.,250.,400.,600.] for y in [.1,.5,2.8/3.85,.95]]
    rows=[]
    for cid,p,t,y in cases:
        row={'id':cid,'p_Pa':p,'T_K':t,'Y_O2':y}
        try:
            s=AbstractState('HEOS','Oxygen' if y==1 else 'Propane' if y==0 else 'Oxygen&Propane')
            if 0<y<1:s.set_mass_fractions([y,1-y])
            s.update(CP.PT_INPUTS,p,t)
            h0=y*refs[0]+(1-y)*refs[1]
            row.update(status='SUCCESS',rho=s.rhomass(),h=s.hmass()-h0,e=s.umass()-h0,
                       h_raw=s.hmass(),e_raw=s.umass(),h_reference=h0,
                       phase_index=int(s.phase()),quality=s.Q())
            if s.phase()==CP.iphase_twophase:
                row['status']='TWO_PHASE_ACOUSTICS_SEPARATE'
                row['liquid_mole_fractions']=s.mole_fractions_liquid()
                row['vapor_mole_fractions']=s.mole_fractions_vapor()
            else:row.update(a=s.speed_sound(),cp=s.cpmass(),cv=s.cvmass())
        except Exception as exc:row.update(status='UNSUPPORTED',detail=str(exc))
        rows.append(row)
    with (out/'probe_inputs.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['id','p_Pa','T_K','Y_O2']);w.writerows(cases)
    report={'coolprop_version':CoolProp.__version__,'coolprop_git':CoolProp.__gitrevision__,
            'backend':'HEOS','h_ideal_298K_Jkg':refs,'rows':rows,
            'limits':['Binary flash failures are recorded, never replaced by pure-fluid averaging.',
                      'Two-phase sound speed requires a separate equilibrium acoustic closure.']}
    if args.rfam:
        legacy={r['id']:r for r in csv.DictReader(Path(args.rfam).open())}
        for r in rows:
            if r['id'] not in legacy:raise RuntimeError('Missing RFAM probe row')
            old=legacy[r['id']];r['rfam']=old
            if r['status'] in ['SUCCESS','TWO_PHASE_ACOUSTICS_SEPARATE'] and old['status']=='SUCCESS':
                r['rfam_error_percent']={k:100*(float(old[k])-r[k])/max(abs(r[k]),1.) for k in ['rho','a','cp','cv','h','e'] if k in r and math.isfinite(float(old[k]))}
        with (out/'comparison.csv').open('w',newline='') as f:
            columns=['id','p_Pa','T_K','Y_O2','cp_status','rfam_status']+[f'{k}_error_pct' for k in ['rho','a','cp','cv','h','e']]
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
            for r in rows:w.writerow({**{k:r[k] for k in ['id','p_Pa','T_K','Y_O2']},'cp_status':r['status'],
                 'rfam_status':r['rfam']['status'],**{k+'_error_pct':v for k,v in r.get('rfam_error_percent',{}).items()}})
    (out/'coolprop_audit.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    for r in rows[:4]:print(r['id'],{k:r.get(k) for k in ['rho','h','e','a','cp','phase_index','rfam_error_percent']})
    print('COOLPROP_AUDIT_WRITTEN',out,'cases',len(rows),'supported',sum(r['status']=='SUCCESS' for r in rows))

if __name__=='__main__':main()
