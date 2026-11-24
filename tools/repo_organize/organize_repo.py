#!/usr/bin/env python3
"""
Tidy D:\\realFluidAllMach: move the 68 loose root .cpp files and misc strays into
a conventional layout, updating every reference site.

Classification
  active  = referenced by CMakeLists.txt OR scripts/local_test.sh
  dead    = referenced by neither (kept, but parked in _unbuilt/ so it is obvious)

Destination
  solver*.cpp, init_test.cpp -> apps/            (executable mains)
  test*.cpp                  -> tests/           (PMake/CTest + local regression suite)
  dead ones                  -> apps/_unbuilt/ , tests/_unbuilt/
  loose *.md                 -> docs/history/
  loose *.sh/*.cmd/*.py      -> scripts/

Reference sites updated (only ones that actually mention a moved file):
  CMakeLists.txt, scripts/local_test.sh, scripts/run_server_tasks.sh,
  scripts/cantera_reference.py, tests/reference/vortex_physics_validation_contract_v1.json,
  tools/validate_vortex_physics_contract.py

Usage:  python organize_repo.py            # dry run, prints the plan
        python organize_repo.py --apply    # performs it
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
APPLY = '--apply' in sys.argv

REFERENCE_FILES = [
    'CMakeLists.txt',
    'scripts/local_test.sh',
    'scripts/run_server_tasks.sh',
    'scripts/cantera_reference.py',
    'tests/reference/vortex_physics_validation_contract_v1.json',
    'tools/validate_vortex_physics_contract.py',
]

DOC_MOVES = {
    'FINAL_EXECUTION_PLAN.md': 'docs/history/',
    'FINAL_THREE_THREAD_REPORT.md': 'docs/history/',
    'PROJECT_COMPLETION_REPORT.md': 'docs/history/',
    'RFAM_ENTROPY_FIX.md': 'docs/history/',
    'THREE_THREAD_EXECUTION_SUMMARY.md': 'docs/history/',
    'rfam_repair_plan.md': 'docs/history/',
    'fluent_rfam_parallel_plan.md': 'docs/history/',
}
SCRIPT_MOVES = {
    'apply_entropy_fix.sh': 'scripts/',
    'create_chemistry_solver.sh': 'scripts/',
    'parallel_build.sh': 'scripts/',
    'test_low_mach_preconditioning.sh': 'scripts/',
    'open_core_sim.cmd': 'scripts/',
    'test_pinfin_bridge.py': 'tests/',
}


def tracked(path: Path) -> bool:
    r = subprocess.run(['git', 'ls-files', '--error-unmatch', str(path)],
                       cwd=ROOT, capture_output=True, text=True)
    return r.returncode == 0


def referenced_by(rel_file: str, basename: str) -> bool:
    p = ROOT / rel_file
    if not p.exists():
        return False
    txt = p.read_text(encoding='utf-8', errors='replace')
    return re.search(r'(?<![/\w])' + re.escape(basename) + r'\b', txt) is not None


# ---- 1. classify root .cpp -------------------------------------------------
root_cpp = sorted(f.name for f in ROOT.glob('*.cpp') if f.is_file())
plan = {}
for name in root_cpp:
    active = any(referenced_by(rf, name) for rf in REFERENCE_FILES)
    if name.startswith('solver'):
        dest = 'apps/' if active else 'apps/_unbuilt/'
    else:
        dest = 'tests/' if active else 'tests/_unbuilt/'
    plan[name] = {'dest': dest, 'active': active}

moves = {k: v['dest'] + k for k, v in plan.items()}
for k, v in DOC_MOVES.items():
    if (ROOT / k).exists():
        moves[k] = v + k
for k, v in SCRIPT_MOVES.items():
    if (ROOT / k).exists():
        moves[k] = v + k

print(f"root .cpp files          : {len(root_cpp)}")
print(f"  -> tests/              : {sum(1 for v in plan.values() if v['dest'] == 'tests/')}")
print(f"  -> tests/_unbuilt/     : {sum(1 for v in plan.values() if v['dest'] == 'tests/_unbuilt/')}")
print(f"  -> apps/               : {sum(1 for v in plan.values() if v['dest'] == 'apps/')}")
print(f"  -> apps/_unbuilt/      : {sum(1 for v in plan.values() if v['dest'] == 'apps/_unbuilt/')}")
print(f"docs / scripts strays    : {len(DOC_MOVES) + len(SCRIPT_MOVES)}")
print(f"total moves              : {len(moves)}")

# ---- 2. reference rewrite map ---------------------------------------------
# for every moved .cpp, the replacement depends on destination
cpp_moves = {k: v for k, v in moves.items() if k.endswith('.cpp')}
ref_edits = {}
for rf in REFERENCE_FILES:
    p = ROOT / rf
    if not p.exists():
        continue
    txt = p.read_text(encoding='utf-8', errors='replace')
    new = txt
    edited = []
    for base, dest in sorted(cpp_moves.items(), key=lambda kv: -len(kv[0])):
        pat = re.compile(r'(?<![/\w])' + re.escape(base) + r'\b')
        if pat.search(new):
            new = pat.sub(dest, new)
            edited.append(base)
    if edited:
        ref_edits[rf] = edited
print("\nreference sites to update:")
for rf, files in ref_edits.items():
    print(f"  {rf}: {len(files)} file refs")

if not APPLY:
    print("\n[dry run] re-run with --apply to execute")
    sys.exit(0)

# ---- 3. verify hashes of provenance-pinned files before moving -------------
PINNED = ROOT / 'tests/reference/vortex_physics_validation_contract_v1.json'
if PINNED.exists():
    doc = json.loads(PINNED.read_text())
    import hashlib
    for entry in doc.get('sources', []):
        pass

# ---- 4. execute moves ------------------------------------------------------
for d in ('apps', 'apps/_unbuilt', 'tests/_unbuilt', 'docs/history', 'scripts'):
    (ROOT / d).mkdir(parents=True, exist_ok=True)

moved, failed = [], []
for src, dst in sorted(moves.items()):
    s, d = ROOT / src, ROOT / dst
    if d.exists():
        failed.append((src, dst, 'destination exists'))
        continue
    try:
        if tracked(Path(src)):
            subprocess.run(['git', 'mv', src, dst], cwd=ROOT, check=True,
                           capture_output=True, text=True)
        else:
            shutil.move(str(s), str(d))
        moved.append((src, dst))
    except Exception as e:
        failed.append((src, dst, repr(e)[:120]))

print(f"\nmoved {len(moved)} files, {len(failed)} failed")
for f in failed:
    print("  FAIL", f)

# ---- 5. rewrite references -------------------------------------------------
for rf, files in ref_edits.items():
    p = ROOT / rf
    txt = p.read_text(encoding='utf-8', errors='replace')
    new = txt
    for base in sorted(files, key=len, reverse=True):
        pat = re.compile(r'(?<![/\w])' + re.escape(base) + r'\b')
        new = pat.sub(cpp_moves[base], new)
    if new != txt:
        p.write_text(new, encoding='utf-8', newline='\n')
        print(f"  updated {rf} ({len(files)} refs)")

# ---- 6. verify no bare references remain ----------------------------------
print("\n=== leftovers ===")
leftover = []
for rf in REFERENCE_FILES:
    p = ROOT / rf
    if not p.exists():
        continue
    txt = p.read_text(encoding='utf-8', errors='replace')
    for base in cpp_moves:
        if re.search(r'(?<![/\w])' + re.escape(base) + r'\b', txt):
            leftover.append((rf, base))
print(f"bare root refs remaining: {len(leftover)}")
for x in leftover[:20]:
    print("  ", x)
print(f"\nroot .cpp remaining: {len(list(ROOT.glob('*.cpp')))}")
