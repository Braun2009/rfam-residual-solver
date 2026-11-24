#!/usr/bin/env python3
"""Compare A-route development fields to the original unseeded prefill.

The centroid-partition flux is a diagnostic surface, not a CAD injector exit.
Pseudo-time iterations are not physical time; no transient mass integration is
inferred from this audit.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tetra_field_to_vtu import read_mesh,read_field


def parse_spec(spec):
    label,path=spec.split('=',1)
    if not label or not path:
        raise ValueError('expected LABEL=PATH')
    return label,Path(path)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as src:
        for chunk in iter(lambda:src.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--baseline',required=True)
    ap.add_argument('--characteristics',action='append',required=True,
                    help='LABEL=corrected-characteristics.json')
    ap.add_argument('--field',action='append',required=True,help='LABEL=field.rfamtetf')
    ap.add_argument('--crossing',action='append',default=[],help='LABEL=centroid-partition CSV')
    ap.add_argument('--board-z',type=float,default=.2709527893066406)
    ap.add_argument('--injector-tip-z',type=float,default=.2858028)
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    characteristics={label:json.loads(path.read_text(encoding='utf-8'))
                     for label,path in map(parse_spec,args.characteristics)}
    fields=dict(map(parse_spec,args.field))
    crossings=dict(map(parse_spec,args.crossing))
    if args.baseline not in characteristics or args.baseline not in fields:
        raise ValueError('baseline must appear in characteristics and fields')
    if set(characteristics)!=set(fields):
        raise ValueError('labels in --field and --characteristics must match')
    _,cells,volume=read_mesh(args.mesh)
    samples={}
    for label,path in fields.items():
        cent,state,y,ns=read_field(path)
        if len(cells)!=len(state) or ns!=2:
            raise ValueError('field/mesh mismatch for '+label)
        downstream=cent[:,2]>=args.board_z
        post_tip=cent[:,2]>=args.injector_tip_z
        oxygen=y[:,0].astype(float)
        rho=state[:,5].astype(float)
        pure=(oxygen>.95)&downstream
        post_tip_pure=(oxygen>.95)&post_tip
        item=characteristics[label]
        if abs(item['injector_tip_z_m']-args.injector_tip_z)>1e-9:
            raise ValueError('field and characteristics use different injector tip')
        info={
            'field':str(path),'corrected_characteristics':item['label'],
            'board_downstream_excess_O2_kg':float(np.sum(volume[downstream]*rho[downstream]*(oxygen[downstream]-.784))),
            'board_downstream_pure_O2_mass_kg':float(np.sum(volume[pure]*rho[pure])),
            'board_downstream_pure_O2_cells':int(np.count_nonzero(pure)),
            'post_tip_pure_O2_mass_kg':float(np.sum(volume[post_tip_pure]*rho[post_tip_pure])),
            'post_tip_pure_O2_cells':int(np.count_nonzero(post_tip_pure)),
            'board_downstream_pure_O2_max_z_mm':float((cent[pure,2].max()-args.board_z)*1e3) if pure.any() else None,
            'eos_two_phase_cells':item['two_phase']['cells'],
            'eos_two_phase_board_downstream_cells':item['two_phase']['board_downstream_cells'],
            'eos_two_phase_post_tip_cells':item['two_phase']['post_tip_cells'],
            'cryogenic_or_dense_proxy_cells':item['cryogenic_or_dense']['cells'],
            'board_core_Mrr_Mzz_motion_ratio':item.get('board_core_motion',{}).get('radial_to_axial_motion_ratio'),
            'board_annulus_YO2_m3_over_mean':item['board_annulus_pattern']['Y_O2']['m3_over_mean'],
            'board_annulus_YO2_m3_over_m1':item['board_annulus_pattern']['Y_O2']['m3_over_m1'],
            'board_annulus_YO2_phase_deg':item['board_annulus_pattern']['Y_O2']['phase_deg']}
        if label in crossings:
            with crossings[label].open(newline='') as f:
                rows=list(csv.DictReader(f))
            info['centroid_partition_flux']={
                'note':'centroid partition, not CAD injector exit; instantaneous from float32 field',
                'faces':len(rows),
                'area_mm2':sum(float(r['area_m2']) for r in rows)*1e6,
                'net_mass_kg_s':sum(float(r['mass_to_chamber_kg_s']) for r in rows),
                'net_O2_kg_s':sum(float(r['O2_to_chamber_kg_s']) for r in rows),
                'net_fuel_kg_s':sum(float(r['fuel_to_chamber_kg_s']) for r in rows)}
        samples[label]=info
    baseline=samples[args.baseline]
    for label,item in samples.items():
        item['delta_from_prefill']={key:item[key]-baseline[key]
            for key in ('board_downstream_excess_O2_kg','board_downstream_pure_O2_mass_kg',
                        'post_tip_pure_O2_mass_kg','board_annulus_YO2_m3_over_mean',
                        'board_annulus_YO2_m3_over_m1','board_core_Mrr_Mzz_motion_ratio')
            if item[key] is not None and baseline[key] is not None}
    report={'status':'unconverged diagnostic, not a CFD design or ignition verdict',
            'baseline':args.baseline,'board_z_m':args.board_z,
            'injector_tip_z_m':args.injector_tip_z,'fields':samples,
            'input_sha256':{'mesh':digest(args.mesh),
                'fields':{label:digest(path) for label,path in fields.items()},
                'phase':{label:digest(str(path)+'.phase.csv')
                         for label,path in fields.items() if Path(str(path)+'.phase.csv').exists()},
                'characteristics':{label:digest(path) for label,path in map(parse_spec,args.characteristics)},
                'crossing':{label:digest(path) for label,path in crossings.items()}},
            'remaining_evidence':['CAD outlet surface and per-stream mass/species/energy/momentum flux',
              'near-field species inventories and mixing/recirculation stable across iterations',
              'signed vector momentum flux through a closed control surface around the axis',
              'EOS phase plus projection residuals with an admitted two-phase transport closure',
              'steady local and boundary residual gates, mesh/model sensitivity']}
    Path(args.out).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for label,data in samples.items():
        print(label,'board_downstream_max_O2_mm',data['board_downstream_pure_O2_max_z_mm'],
              'board_downstream_excess_O2_mg',data['board_downstream_excess_O2_kg']*1e6,
              'post_tip_pureO2_cells',data['post_tip_pure_O2_cells'],
              'board_annulus_m3/mean',data['board_annulus_YO2_m3_over_mean'],
              'EOS_twophase',data['eos_two_phase_cells'],
              'cross_mass',data.get('centroid_partition_flux',{}).get('net_mass_kg_s'))
    print('COLD_FLOW_DEVELOPMENT_AUDIT_WRITTEN',args.out)


if __name__=='__main__':
    main()
