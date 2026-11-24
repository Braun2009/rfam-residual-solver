#!/usr/bin/env python3
"""Audit Core Torch mass allocation without modifying CFD inputs or geometry."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def designer_defaults(path: Path) -> dict[str, float]:
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    config = next((node for node in tree.body
                   if isinstance(node, ast.ClassDef) and node.name == 'VortexConfig'), None)
    if config is None:
        raise ValueError('VortexConfig not found')
    values = {}
    for node in config.body:
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            try:
                values[node.target.id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
    names = ('mdot_LOX', 'mdot_C3H8_main', 'mdot_C3H8_film',
             'core_torch_frac', 'n_pintles', 'OF')
    if not all(name in values for name in names):
        raise ValueError('missing designer flow allocation inputs')
    result = {name: values[name] for name in names}
    if (not all(isinstance(result[name], (int, float)) and not isinstance(result[name], bool)
                and math.isfinite(result[name]) for name in names)
        or result['mdot_LOX'] <= 0 or result['mdot_C3H8_main'] <= 0
        or result['mdot_C3H8_film'] < 0 or result['OF'] <= 0
        or not 0 < result['core_torch_frac'] < 1
        or type(result['n_pintles']) is not int or result['n_pintles'] <= 0):
        raise ValueError('invalid designer flow allocation inputs')
    return result


def audit(designer: Path, old_design: Path, mesh_audit: Path) -> dict:
    values = designer_defaults(designer)
    old = json.loads(old_design.read_text(encoding='utf-8'))
    mesh = json.loads(mesh_audit.read_text(encoding='utf-8'))
    lox = values['mdot_LOX']
    fuel = values['mdot_C3H8_main']
    fraction = values['core_torch_frac']
    count = values['n_pintles']
    torch_total = (lox + fuel) * fraction
    torch_lox = torch_total * values['OF'] / (1 + values['OF'])
    torch_fuel = torch_total / (1 + values['OF'])
    pintle_lox, pintle_fuel = lox - torch_lox, fuel - torch_fuel
    # The original designer's O/F is rounded to seven significant figures.
    if not math.isclose(lox / fuel, values['OF'], rel_tol=1e-6):
        raise ValueError('designer O/F and main-feed totals disagree')
    old_pintle_lox = old['per_pintle']['LOX_core']['mdot_kg_s'] * count
    old_pintle_fuel = old['per_pintle']['C3H8_ring']['mdot_kg_s'] * count
    old_torch = old['torch_center']['mdot_kg_s']
    old_lox, old_fuel = old['totals_kg_s']['LOX'], old['totals_kg_s']['C3H8_main']
    if not all(math.isfinite(x) and x >= 0 for x in
               (old_pintle_lox, old_pintle_fuel, old_torch, old_lox, old_fuel)):
        raise ValueError('invalid old design mass flows')
    zones = mesh['boundary_zones']
    expected = {10: 'lox_0', 11: 'lox_1', 12: 'lox_2',
                20: 'fuel_0', 21: 'fuel_1', 22: 'fuel_2', 30: 'torch'}
    if any(zones.get(str(zone)) != label for zone, label in expected.items()):
        raise ValueError('mesh inlet zone mapping differs from allocation audit')
    patches = mesh['inlets']
    counts = {family: sum(p['design_family'] == family for p in patches)
              for family in ('lox', 'fuel', 'torch')}
    if counts != {'lox': 3, 'fuel': 60, 'torch': 1}:
        raise ValueError('mesh inlet patch counts differ from allocation audit')
    old_excess = (old_pintle_lox + old_pintle_fuel + old_torch) - (old_lox + old_fuel)
    if not math.isclose(old_lox, lox, abs_tol=1e-10) or not math.isclose(old_fuel, fuel, abs_tol=1e-10):
        raise ValueError('old-design main totals differ from designer')
    return {
        'schema': 'rfam.injector-flow-allocation-audit/v1',
        'status': 'CORE_BRANCH_CONFIRMED_OPERATING_DETAILS_PENDING',
        'scope': 'design mass budget only; no CFD boundary, geometry, or ignition authorization',
        'sources': {name: {'path': str(path), 'sha256': sha256(path)}
                    for name, path in [('designer', designer), ('old_design_json', old_design),
                                       ('mesh_audit', mesh_audit)]},
        'authority': 'User confirms Core Torch must be included in flow allocation. Subtraction from main LOX and C3H8 totals and the 3.5% fraction follow the existing designer, not a newly measured or approved operating schedule.',
        'assumed_core_fraction': fraction,
        'design_kg_s': {'lox_total': lox, 'main_fuel_total': fuel,
                        'film_fuel_separate': values['mdot_C3H8_film'],
                        'main_feed_total_excluding_film': lox + fuel},
        'core_torch_kg_s': {'lox': torch_lox, 'main_fuel': torch_fuel,
                            'total': torch_total},
        'pintles_kg_s': {'lox_total': pintle_lox, 'main_fuel_total': pintle_fuel,
                         'lox_per_pintle': pintle_lox / count,
                         'main_fuel_per_pintle': pintle_fuel / count},
        'mesh_inlet_zones': {label: zone for zone, label in expected.items()},
        'mesh_inlet_patch_counts': counts,
        'old_json': {
            'claimed_main_feed_total_kg_s': old_lox + old_fuel,
            'pintles_plus_torch_kg_s': old_pintle_lox + old_pintle_fuel + old_torch,
            'overallocation_kg_s': old_excess,
            'double_counts_torch': old_excess > 0.1 * torch_total,
            'precision_note': 'Original JSON flows are rounded to four decimal places.'},
        'current_A_diagnostic': {
            'pintle_lox_kg_s': lox, 'pintle_fuel_kg_s': fuel,
            'torch_kg_s': 0, 'full_branch_allocation_represented': False},
        'unresolved': ['confirm 3.5% torch fraction and valve timing',
                       'independent LOX and propane torch inlet geometry and supply states',
                       'cold-flow model with both physical central feeds',
                       'named true exits and verified outlet fluxes'],
        'production_validated': False,
        'ignition_authorized': False,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--designer', required=True, type=Path)
    ap.add_argument('--old-design-json', required=True, type=Path)
    ap.add_argument('--mesh-audit', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    result = audit(args.designer, args.old_design_json, args.mesh_audit)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    print('INJECTOR_FLOW_ALLOCATION', result['status'], args.out)


if __name__ == '__main__':
    main()
