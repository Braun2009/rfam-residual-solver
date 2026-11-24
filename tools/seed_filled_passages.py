#!/usr/bin/env python3
"""Create a development seed with only audited internal passages filled.

Not a transient startup prediction: external fluid cells are copied verbatim
from the unseeded step-0 field, and the cold torch remains off. Results seeded
this way cannot validate self-filling of the injector channels; they need an
independent steady-state and grid/initialization comparison.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import struct
import numpy as np
from CoolProp import CoolProp as CP
from CoolProp.CoolProp import AbstractState


def checksum(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as file:
        for part in iter(lambda:file.read(1<<20),b''):h.update(part)
    return h.hexdigest()


def seed(mesh_path,field_path,mask_path,out_path):
    identity=checksum(mesh_path)
    with np.load(mask_path) as data:
        tags=np.array(data['tag'],copy=True)
        stored=str(data['mesh_sha256'])
    if stored!=identity:raise ValueError('mask mesh SHA does not match selected geometry')
    raw=Path(field_path).read_bytes()
    if raw[:9]!=b'RFAMTETF1':raise ValueError('field magic mismatch')
    cells,species=struct.unpack_from('<QQ',raw,9)
    if len(tags)!=cells or species!=2:raise ValueError('mesh/mask/field shape mismatch')
    with Path(mesh_path).open('rb') as file:
        header=file.read(32)
    np_,nc,nf=struct.unpack_from('<QQQ',header,8)
    if nc!=cells:raise ValueError('mesh/field cell counts differ')
    fields=np.frombuffer(raw,'<f4',count=8*cells,offset=25+12*cells).reshape(cells,8).copy()
    y=np.frombuffer(raw,'<f4',count=2*cells,offset=25+44*cells).reshape(cells,2).copy()
    # Identical sensible-enthalpy origin as CoolPropThermoModel.C:
    # ideal-gas h_i(298.15 K) weighted by mass fraction.
    offset={}
    for fluid in ('Oxygen','Propane'):
        sample=AbstractState('HEOS',fluid)
        sample.update(CP.DmassT_INPUTS,1e-10,298.15)
        offset[fluid]=sample.hmass_idealgas()
    families={
        'lox':dict(tag=1,fluid='Oxygen',p=4.2045e6,T=129.48,mdot_per_zone=2.8/3,
                   inlet_area_m2=(49.545+49.5514+49.3312)/3*1e-6,Y=(1.,0.)),
        'fuel':dict(tag=2,fluid='Propane',p=4.1295e6,T=448.2,mdot_per_zone=1.05/3,
                    inlet_area_m2=(62.6356+62.635+62.6042)/3*1e-6,Y=(0.,1.))}
    summary={}
    for name,d in families.items():
        chosen=tags==d['tag']
        sample=AbstractState('HEOS',d['fluid'])
        sample.update(CP.PT_INPUTS,d['p'],d['T'])
        rho=sample.rhomass()
        velocity=d['mdot_per_zone']/(rho*d['inlet_area_m2'])
        fields[chosen,0]=d['p'];fields[chosen,1]=d['T']
        fields[chosen,2]=0.;fields[chosen,3]=0.;fields[chosen,4]=velocity
        fields[chosen,5]=rho;fields[chosen,6]=sample.speed_sound()
        fields[chosen,7]=sample.hmass()-offset[d['fluid']]
        y[chosen,0],y[chosen,1]=d['Y']
        summary[name]={'cells':int(chosen.sum()),'p_Pa':d['p'],'T_K':d['T'],
            'rho_kg_m3':rho,'u_z_m_s':velocity,'sound_m_s':sample.speed_sound(),
            'h_sensible_J_kg':sample.hmass()-offset[d['fluid']],
            'species':list(d['Y']),'p_inlet_design_kg_s':d['mdot_per_zone']}
    untouched=(tags==0)|(tags==3)
    if not np.array_equal(fields[untouched],
        np.frombuffer(raw,'<f4',count=8*cells,offset=25+12*cells).reshape(cells,8)[untouched]):
        raise AssertionError('exterior state changed')
    if not np.array_equal(y[untouched],
        np.frombuffer(raw,'<f4',count=2*cells,offset=25+44*cells).reshape(cells,2)[untouched]):
        raise AssertionError('exterior composition changed')
    if not np.any(tags==3):raise AssertionError('torch mask missing')
    out=Path(out_path);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('wb') as file:
        file.write(raw[:25+12*cells])
        file.write(fields.astype('<f4',copy=False).tobytes())
        file.write(y.astype('<f4',copy=False).tobytes())
    meta={'purpose':'topology-restricted filled-channel STEADY development seed; NOT startup or ignition result',
          'mesh_sha256':identity,'mask_sha256':checksum(mask_path),
          'source_field_sha256':checksum(field_path),'output_field_sha256':checksum(out),
          'cells_unchanged_exterior_and_torch':int(np.count_nonzero((tags==0)|(tags==3))),
          'torch_on':False,'families':summary,
          'limitations':['all filled-channel velocity initialized axially; actual LOX outlet is oblique',
                         'sharp interface at channel end may cause EOS projection/stiffness',
                         'not independently validated against a second initialization or mesh']}
    Path(str(out)+'.seed.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return meta


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--field',required=True)
    ap.add_argument('--mask',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    meta=seed(args.mesh,args.field,args.mask,args.out)
    print('PASSAGE_SEED_WRITTEN',args.out,'lox',meta['families']['lox']['cells'],
          'fuel',meta['families']['fuel']['cells'],'torch_on',meta['torch_on'])


if __name__=='__main__':main()
