#!/usr/bin/env python3
"""Inventory propellant filling in the audited, disjoint internal passages.

A pure-stream Y threshold measures a tagged concentration front, not a
physical-time propagation speed. No claim is made about post-exit mixing.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from tetra_field_to_vtu import read_mesh,read_field

FAMILIES={'lox':(1,280.9),'fuel':(2,276.8),'torch':(3,283.3)}


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as file:
        for chunk in iter(lambda:file.read(1<<20),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--mask',required=True,help='NPZ from audit_injector_passage_connectivity.py')
    ap.add_argument('--field',action='append',required=True,help='LABEL=FIELD.rfamtetf')
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    _,cells,volume=read_mesh(args.mesh)
    with np.load(args.mask) as data:
        tags=data['tag'].copy()
        identity=str(data['mesh_sha256'])
    mesh_sha=sha256(args.mesh)
    if identity!=mesh_sha or len(tags)!=len(cells):
        raise ValueError('passage mask does not match the authoritative mesh')
    result={'scope':'unconverged internal-passage inventory only, not outlet or ignition',
            'mesh_sha256':mesh_sha,'mask_sha256':sha256(args.mask),'fields':{}}
    for entry in args.field:
        label,path=entry.split('=',1)
        cent,f,y,ns=read_field(path)
        if len(f)!=len(tags) or ns!=2:raise ValueError('field and mask mismatch')
        z=cent[:,2].astype(float)
        yo=y[:,0].astype(float)
        data={}
        for family,(tag,limit) in FAMILIES.items():
            selected=tags==tag
            pure=selected & (yo>.95 if family=='lox' else yo<.2 if family=='fuel' else np.zeros(len(yo),dtype=bool))
            last=selected & (z>limit*1e-3-5e-4)
            premerge_pure=last&pure
            last_bg=last&(np.abs(yo-.784)<1e-5)
            threshold_values=([.8,.85,.9,.95,.99] if family=='lox'
                              else [.7,.6,.4,.2,.1] if family=='fuel' else [])
            thresholds=[]
            for cutoff in threshold_values:
                matched=selected&((yo>cutoff) if family=='lox' else (yo<cutoff))
                thresholds.append({'Y_O2_condition':('>' if family=='lox' else '<')+str(cutoff),
                    'cells':int(np.count_nonzero(matched)),
                    'furthest_centroid_z_mm':float(z[matched].max()*1e3) if matched.any() else None})
            data[family]={'cells':int(np.count_nonzero(selected)),
                'threshold_fronts':thresholds,
                'volume_cm3':float(np.sum(volume[selected])*1e6),
                'pure_stream_cells':int(np.count_nonzero(pure)),
                'pure_stream_volume_fraction':float(np.sum(volume[pure])/np.sum(volume[selected])),
                'furthest_pure_centroid_z_mm':float(z[pure].max()*1e3) if pure.any() else None,
                'gap_to_premerge_cut_mm':float(limit-z[pure].max()*1e3) if pure.any() else None,
                'last_half_mm_cells':int(np.count_nonzero(last)),
                'last_half_mm_pure_cells':int(np.count_nonzero(premerge_pure)),
                'last_half_mm_background_cells':int(np.count_nonzero(last_bg)),
                'max_pressure_MPa':float(np.max(f[selected,0])/1e6)}
        result['fields'][label]={'source':path,'sha256':sha256(path),'family':data}
        print(label,[(k,round(v['pure_stream_volume_fraction']*100,2),
                      v['gap_to_premerge_cut_mm'],v['last_half_mm_pure_cells'])
                     for k,v in data.items()])
    Path(args.out).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('PASSAGE_FILL_INVENTORY_WRITTEN',args.out)


if __name__=='__main__':main()
