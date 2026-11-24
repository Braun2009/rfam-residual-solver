#!/usr/bin/env bash
# Measure one solver step with the full San Diego mechanism on the A grid.
set -uo pipefail
BIN=${BIN:-./build-full/rfam_v20_tetra_react}
MESH=${MESH:-work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet}
STEPS=${STEPS:-2}
t0=$(date +%s)
$BIN --mesh "$MESH" --steps "$STEPS" --monitor-every 1 --blueprint \
  --thermo legacy --inlet-mode mass --p0 4e6 --T0 3400 --guard-T-max 4300 \
  --lox-pinj 6e6 --lox-T 300 --lox-mdot 2.70 \
  --fuel-pinj 4.5e6 --fuel-T 450 --fuel-mdot 1.01 --torch-mdot 0 \
  --reacting --mechanism "${MECH:-sandiego}" --chem-max-substeps 2 \
  --molecular-viscous --turbulence sst \
  --lusgs --cfl 0.3 --lusgs-dt-scale 0.05 > /tmp/sd_step.log 2>&1
t1=$(date +%s)
n=$(grep -c '^STEP' /tmp/sd_step.log || true)
echo "MECH=${MECH:-sandiego} STEP_lines=$n WALL=$((t1-t0))s"
grep -E 'chemistry_mechanism|molecular_viscous=wilke|SandiegoReduced|San Diego' /tmp/sd_step.log | head -3
grep '^STEP' /tmp/sd_step.log | tail -1 | cut -c1-150
grep -E 'FAILED|REJECTED' /tmp/sd_step.log | tail -2 | cut -c1-160
