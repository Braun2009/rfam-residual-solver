#!/usr/bin/env python3
"""Check the Core Torch design split and detect the old extra-torch budget."""
from __future__ import annotations

import copy
import json
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from audit_injector_flow_allocation import audit
from build_injector_boundary_decision_package import verify_flow_allocation


def close(actual: float, expected: float) -> None:
    assert math.isclose(actual, expected, abs_tol=1e-7, rel_tol=0), (actual, expected)


def main() -> None:
    with tempfile.TemporaryDirectory() as directory:
        base = Path(directory)
        designer = base / 'designer.py'
        designer.write_text('''class VortexConfig:
    mdot_LOX: float = 2.8
    mdot_C3H8_main: float = 1.05
    mdot_C3H8_film: float = 0.45
    core_torch_frac: float = 0.035
    n_pintles: int = 3
    OF: float = 2.666667
''', encoding='utf-8')
        old = base / 'old.json'
        old.write_text(json.dumps({
            'totals_kg_s': {'LOX': 2.8, 'C3H8_main': 1.05},
            'per_pintle': {'LOX_core': {'mdot_kg_s': 0.9333},
                           'C3H8_ring': {'mdot_kg_s': 0.35}},
            'torch_center': {'mdot_kg_s': 0.1348},
        }), encoding='utf-8')
        mesh = base / 'mesh.json'
        zones = {str(zone): label for zone, label in {
            10: 'lox_0', 11: 'lox_1', 12: 'lox_2',
            20: 'fuel_0', 21: 'fuel_1', 22: 'fuel_2', 30: 'torch',
        }.items()}
        mesh.write_text(json.dumps({
            'boundary_zones': zones,
            'inlets': ([{'design_family': 'lox'}] * 3 +
                       [{'design_family': 'fuel'}] * 60 +
                       [{'design_family': 'torch'}]),
        }), encoding='utf-8')
        result = audit(designer, old, mesh)
        close(result['core_torch_kg_s']['lox'], 0.098)
        close(result['core_torch_kg_s']['main_fuel'], 0.03675)
        close(result['core_torch_kg_s']['total'], 0.13475)
        close(result['pintles_kg_s']['lox_total'] + result['core_torch_kg_s']['lox'], 2.8)
        close(result['pintles_kg_s']['main_fuel_total'] + result['core_torch_kg_s']['main_fuel'], 1.05)
        close(result['pintles_kg_s']['main_fuel_per_pintle'], 0.33775)
        close(result['design_kg_s']['film_fuel_separate'], 0.45)
        assert result['old_json']['double_counts_torch']
        close(result['old_json']['overallocation_kg_s'], 0.1347)
        assert not result['current_A_diagnostic']['full_branch_allocation_represented']
        assert not result['ignition_authorized'] and not result['production_validated']
        verify_flow_allocation(result)
        broken = copy.deepcopy(result)
        broken['core_torch_kg_s']['lox'] += 0.1
        try:
            verify_flow_allocation(broken)
        except ValueError as error:
            assert 'mass budget does not close' in str(error)
        else:
            raise AssertionError('tampered Core Torch mass budget was accepted')
        broken = copy.deepcopy(result)
        broken['core_torch_kg_s']['lox'] += 0.01
        broken['pintles_kg_s']['lox_total'] -= 0.01
        broken['pintles_kg_s']['lox_per_pintle'] -= 0.01 / 3
        broken['core_torch_kg_s']['main_fuel'] -= 0.01
        broken['pintles_kg_s']['main_fuel_total'] += 0.01
        broken['pintles_kg_s']['main_fuel_per_pintle'] += 0.01 / 3
        try:
            verify_flow_allocation(broken)
        except ValueError as error:
            assert 'differs from source audit' in str(error)
        else:
            raise AssertionError('source-inconsistent but balanced torch split was accepted')
        broken = copy.deepcopy(result)
        broken['sources']['designer']['sha256'] = '0' * 64
        try:
            verify_flow_allocation(broken)
        except ValueError as error:
            assert 'stale Core Torch' in str(error)
        else:
            raise AssertionError('stale designer identity was accepted')
        zones['30'] = 'fuel_3'
        mesh.write_text(json.dumps({'boundary_zones': zones,
                                    'inlets': [{'design_family': 'lox'}] * 3 +
                                              [{'design_family': 'fuel'}] * 60 +
                                              [{'design_family': 'torch'}]}), encoding='utf-8')
        try:
            audit(designer, old, mesh)
        except ValueError as error:
            assert 'zone mapping differs' in str(error)
        else:
            raise AssertionError('incorrect Core Torch zone mapping was accepted')
    print('INJECTOR_FLOW_ALLOCATION_AUDIT_PASS')


if __name__ == '__main__':
    main()
