#!/usr/bin/env python3
"""Build a machine-readable external decision package for A cold-flow gates."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

from audit_injector_flow_allocation import audit as audit_flow_allocation


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def verify_flow_allocation(allocation):
    if (allocation.get('schema') != 'rfam.injector-flow-allocation-audit/v1'
            or allocation.get('status') != 'CORE_BRANCH_CONFIRMED_OPERATING_DETAILS_PENDING'
            or allocation.get('production_validated') is not False
            or allocation.get('ignition_authorized') is not False):
        raise ValueError('unverified Core Torch flow allocation audit')
    try:
        sources=allocation['sources']
        for name in ('designer','old_design_json','mesh_audit'):
            if sha(sources[name]['path']) != sources[name]['sha256']:
                raise ValueError('stale Core Torch flow allocation audit: '+name)
        design=allocation['design_kg_s']
        core=allocation['core_torch_kg_s']
        pintles=allocation['pintles_kg_s']
        fraction=allocation['assumed_core_fraction']
        close=lambda x,y: (isinstance(x,(int,float)) and not isinstance(x,bool)
                           and isinstance(y,(int,float)) and not isinstance(y,bool)
                           and math.isfinite(x) and math.isfinite(y)
                           and math.isclose(x,y,rel_tol=1e-7,abs_tol=1e-9))
        if (not 0 < fraction < 1
                or not close(core['lox']+pintles['lox_total'],design['lox_total'])
                or not close(core['main_fuel']+pintles['main_fuel_total'],design['main_fuel_total'])
                or not close(core['lox']+core['main_fuel'],core['total'])
                or not close(core['total'],fraction*design['main_feed_total_excluding_film'])
                or not close(design['lox_total']+design['main_fuel_total'],
                             design['main_feed_total_excluding_film'])
                or not close(pintles['lox_total'],3*pintles['lox_per_pintle'])
                or not close(pintles['main_fuel_total'],3*pintles['main_fuel_per_pintle'])
                or allocation['old_json']['double_counts_torch'] is not True
                or allocation['current_A_diagnostic']['full_branch_allocation_represented'] is not False):
            raise ValueError('Core Torch mass budget does not close')
        expected=audit_flow_allocation(*(Path(sources[name]['path']) for name in
                                          ('designer','old_design_json','mesh_audit')))
        if allocation != expected:
            raise ValueError('Core Torch allocation differs from source audit')
    except (KeyError, TypeError, OSError) as exc:
        raise ValueError('incomplete Core Torch flow allocation audit') from exc


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True)
    ap.add_argument('--state-definition',required=True)
    ap.add_argument('--exit-semantics',required=True)
    ap.add_argument('--fuel-envelope',required=True)
    ap.add_argument('--flow-allocation',required=True)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--design-json',required=True)
    ap.add_argument('--geometry-manifest',required=True)
    ap.add_argument('--chemistry-status', default=None)
    args=ap.parse_args()
    state=json.loads(Path(args.state_definition).read_text(encoding='utf-8'))
    exit_gate=json.loads(Path(args.exit_semantics).read_text(encoding='utf-8'))
    envelope=json.loads(Path(args.fuel_envelope).read_text(encoding='utf-8'))
    allocation=json.loads(Path(args.flow_allocation).read_text(encoding='utf-8'))
    verify_flow_allocation(allocation)
    chemistry = json.loads(Path(args.chemistry_status).read_text(encoding='utf-8')) if args.chemistry_status else {
        'status': 'NOT_A_RUNTIME_IGNITION_TEST',
        'sundials_available': False,
        'tetra_sandiego_59_species_coupling': False,
        'torch_or_spark_source': False,
        'source': 'not-run-local'}
    requested={
        'boundary_definition': {
            'status': state.get('boundary_definition_status'),
            'question': 'Are listed injector p/T static or stagnation values?',
            'required': ['static_or_stagnation_for_each_supply',
                         'fuel_manifold_pressure_temperature_history',
                         'fixed_mdot_is_design_or_supply_outcome',
                         'pressure_drop_definition_and_loss_model'],
        },
        'exit_surface_semantics': {
            'status': 'MISSING' if not exit_gate.get('true_exit_surfaces_available') else 'AVAILABLE',
            'required_labels': exit_gate.get('required_exit_labels'),
            'question': 'Provide named/geometry-mapped curved LOX/fuel/torch exit surfaces.',
            'not_accepted_as_exit': ['axial z cuts', 'pre-exit channel components',
                                     'z thresholds', 'largest exterior connected component'],
        },
        'operating_point': {
            'fuel_envelope_scope': envelope.get('scope'),
            'fuel_envelope_case': 'historical A with torch=0 and 0.35 kg/s per fuel zone; not an audit of torch-allocated 0.33775 kg/s per zone',
            'current_rows': envelope.get('rows'),
            'question': 'Approve the supply-state contract or provide authoritative manifold histories.',
        },
        'core_torch_flow_allocation': {
            'status': allocation['status'],
            'authority': allocation['authority'],
            'design_kg_s': allocation['design_kg_s'],
            'assumed_core_fraction': allocation['assumed_core_fraction'],
            'core_torch_kg_s': allocation['core_torch_kg_s'],
            'pintles_kg_s': allocation['pintles_kg_s'],
            'old_json': allocation['old_json'],
            'current_A_diagnostic': allocation['current_A_diagnostic'],
            'unresolved': allocation['unresolved'],
        },
        'chemistry_readiness': chemistry,
        'production_gates_before_ignition': [
            'accepted physical cold flow beyond true exits',
            'per-outlet mass/species/energy/vector-momentum budgets',
            'external mixing and backflow metrics',
            'phase/transport closure and mesh sensitivity',
            'five stable residual/budget samples',
            'chemistry state mapping with actual torch/spark inputs',
        ],
    }
    files={}
    for label,path in [('state_definition',args.state_definition),('exit_semantics',args.exit_semantics),
                        ('fuel_envelope',args.fuel_envelope),('flow_allocation',args.flow_allocation),
                        ('mesh',args.mesh),
                        ('design_json',args.design_json),('geometry_manifest',args.geometry_manifest)]:
        files[label]={'path':str(path),'sha256':sha(path)}
    report={'schema':'rfam.injector-boundary-decision/v1','status':'BLOCKED_PENDING_EXTERNAL_AUTHORITY',
            'files':files,'decision_requests':requested,
            'current_claim':'No ignition test or reacting injector run is authorized from current fields.'}
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('INJECTOR_BOUNDARY_DECISION_PACKAGE',report['status'],out)


if __name__=='__main__':main()
