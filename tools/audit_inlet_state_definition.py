#!/usr/bin/env python3
"""Audit whether injector pressure/temperature inputs define static or total state.

This does not choose an inlet model. It reports the source wording and marks the
boundary definition unresolved when only pressure, temperature, mdot and an
injection pressure drop are provided without static/total labels.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def matches(text: str, pattern: str):
    return [m.group(0) for m in re.finditer(pattern, text, flags=re.IGNORECASE)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', action='append', required=True,
                    help='design text/JSON source; repeatable')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    sources = []
    for raw in args.source:
        path = Path(raw)
        text = path.read_text(encoding='utf-8', errors='replace')
        sources.append({
            'path': str(path),
            'sha256': sha256(path),
            'pressure_mentions': matches(text, r'.{0,45}(?:入口压强|入口压力|supply pressure|pressure).{0,80}'),
            'temperature_mentions': matches(text, r'.{0,45}(?:入口温度|注入温度|temperature).{0,80}'),
            'mass_flow_mentions': matches(text, r'.{0,45}(?:质量流量|mass.?flow|mdot).{0,80}'),
            'pressure_drop_mentions': matches(text, r'.{0,45}(?:喷注压降|pressure drop|delta.?p).{0,80}'),
            'total_state_labels': sorted(set(re.findall(
                r'\b(?:total|stagnation|滞止|总态|总压|总温|总焓|总熵)\b', text, flags=re.I))),
            'static_state_labels': sorted(set(re.findall(
                r'\b(?:static|静态|静压|静温|静焓|静熵)\b', text, flags=re.I))),
        })
    total = sorted(set(x for s in sources for x in s['total_state_labels']))
    static = sorted(set(x for s in sources for x in s['static_state_labels']))
    report = {
        'scope': 'input-definition audit only; no inlet model selected and no CFD input changed',
        'sources': sources,
        'total_state_labels_found': total,
        'static_state_labels_found': static,
        'boundary_definition_status': 'RESOLVED' if total or static else 'UNRESOLVED',
        'reason': ('The supplied design text contains pressure/temperature/flow and pressure-drop '
                   'values but no explicit static-vs-total labels.' if not total and not static else
                   'At least one explicit state label was found; cross-source consistency still requires review.'),
        'required_authoritative_inputs': [
            'whether listed inlet pressure and temperature are static or stagnation state',
            'whether the 0.5 MPa fuel pressure drop is a supply-to-chamber static drop or an injector loss',
            'whether the fixed design mass flow is imposed or an outcome of a supply/discharge model',
            'fuel manifold pressure/temperature history or an approved equivalent boundary contract',
        ],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print('INLET_STATE_DEFINITION', report['boundary_definition_status'])
    print(out)


if __name__ == '__main__':
    main()
