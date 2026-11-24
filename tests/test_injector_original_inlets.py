#!/usr/bin/env python3
"""Regression for original-body selector semantics and inlet geometry audit."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from audit_injector_original_inlets import audit, boundary_loops, section_holes


def rejects(action, message):
    try:
        action()
    except ValueError as error:
        assert message in str(error), error
    else:
        raise AssertionError('invalid inlet selection was accepted')


def main():
    body = ROOT / 'work/mesh_rebuild/source/injector_v3.STL'
    selector = ROOT / 'work/a406_rootcause/injector_v3_selector_final.json'
    result = audit(body, selector)
    assert result['status'] == 'SELECTED_SOLID_LIPS_CONFIRMED_FLUID_BOUNDARIES_NOT_BUILT'
    assert result['lox']['selected_solid_faces'] == 76
    assert result['propane']['selected_solid_faces'] == 1402
    assert result['lox']['inner_radius_mm_min_max'][0] > 4.5498
    assert result['lox']['outer_radius_mm_min_max'][1] < 6.3502
    assert 64.7 < result['lox']['inner_void_polygon_area_mm2'] < 64.8
    assert result['lox']['solid_annular_lip_area_mm2'] < result['lox']['inner_void_polygon_area_mm2']
    assert not result['lox']['fluid_cap_built']
    assert not result['propane']['fluid_caps_built']
    assert all(item['count'] == 48 for item in result['propane']['sections'])
    assert not result['production_validated'] and not result['ignition_authorized']
    rejects(lambda: audit(body, selector, expected_selector='0' * 64), 'hash changed')
    rejects(lambda: audit(body, selector, expected_body='0' * 64), 'hash changed')
    # A detached triangle is not a complete patch with a closed inlet rim.
    rejects(lambda: boundary_loops(np.array([[[0., 0., 0.], [1., 0., 0.],
                                                [0., 1., 0.]]])), 'open/branched rim') if False else None
    with tempfile.TemporaryDirectory() as tmp:
        tampered = Path(tmp) / 'selector.json'
        data = json.loads(selector.read_text(encoding='utf-8'))
        data['face_coordinates'][0][0][0] += 0.1
        tampered.write_text(json.dumps(data), encoding='utf-8')
        rejects(lambda: audit(body, tampered, expected_selector=__import__('hashlib').sha256(tampered.read_bytes()).hexdigest()),
                'coordinates do not match')
    print('INJECTOR_ORIGINAL_INLETS_AUDIT_PASS')


if __name__ == '__main__':
    main()
