#!/usr/bin/env python3
"""Mixed-state recovery audit.

Feeds the material-contact mixed conserved states printed by
`audit_tetra_linearization --contact` through the production rho-e recovery
(`probe_rhoe_recovery`, seeded from the 4 MPa/400 K ambient state) and records
which fractions the equilibrium flash rejects. Same inputs as the in-process
contact test; this driver separates the failure inventory from the pressure
error table.
"""
import argparse
import json
import re
import subprocess
from pathlib import Path

PROBE = '/mnt/d/realFluidAllMach/build-phase-foundation/probe_rhoe_recovery'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contact-log', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--seed-p', type=float, default=4e6)
    parser.add_argument('--seed-T', type=float, default=400.0)
    args = parser.parse_args()

    rows = []
    for line in Path(args.contact_log).read_text().splitlines():
        m = re.match(r'MATERIAL_CONTACT_INPUT kind=(\S+) volume_fraction=(\S+) '
                     r'rho=(\S+) e=(\S+) yO2=(\S+)', line)
        if not m:
            continue
        kind, frac, rho, energy, yO2 = (m.group(1), float(m.group(2)),
                                       float(m.group(3)), float(m.group(4)),
                                       float(m.group(5)))
        command = ['wsl', 'bash', '-lc',
                   f'{PROBE} {rho:.17g} {energy:.17g} {yO2:.17g} '
                   f'{args.seed_p:.17g} {args.seed_T:.17g} run']
        proc = subprocess.run(command, capture_output=True, text=True, timeout=120)
        text = proc.stdout
        success = re.search(r'MODEL_RESULT success=(\d)', text)
        residuals = re.search(r'MODEL_RESIDUALS density=(\S+) energy=(\S+)', text)
        native = re.search(r'NATIVE_RHOE success=(\d)(?: error=(.*))?', text)
        state = re.search(r'MODEL_STATE p=(\S+) T=(\S+) rho=(\S+) e=(\S+)', text)
        rows.append(dict(
            kind=kind, volume_fraction=frac, rho=rho, e=energy, yO2=yO2,
            recovery_success=bool(success and success.group(1) == '1'),
            recovery_residuals=({'density': float(residuals.group(1)),
                                 'energy': float(residuals.group(2))}
                                if residuals else None),
            recovered_state=({'p': float(state.group(1)), 'T': float(state.group(2)),
                              'rho': float(state.group(3)), 'e': float(state.group(4))}
                             if state else None),
            native_dmassumass_success=bool(native and native.group(1) == '1'),
            native_error=(native.group(2).strip() if native and native.group(2) else None),
        ))
    rejected = [row for row in rows if not row['recovery_success']]
    report = {'probe': PROBE, 'seed': {'p': args.seed_p, 'T': args.seed_T},
              'states': rows, 'rejected': len(rejected), 'total': len(rows),
              'note': ('native CoolProp mixture DmassUmass is not implemented '
                       '(PHSU_D_flash not ready for mixtures); the wrapper Newton '
                       'is the only mixture rho-e path')}
    Path(args.out).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('MIXED_STATE_RECOVERY', f'{len(rows) - len(rejected)}/{len(rows)} ok')
    for row in rejected:
        print('  rejected', row['kind'], 'f=%.3g' % row['volume_fraction'],
              'yO2=%.3f' % row['yO2'], 'rho=%.6g' % row['rho'])
    print(args.out)


if __name__ == '__main__':
    main()
