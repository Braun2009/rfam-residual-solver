#!/usr/bin/env bash
# Per-step wall-clock comparison: legacy (6-species reacting, Wilke + SST)
# against CoolProp (2-species real fluid, HEOS transport + SST) on one grid.
# Both runs use the same mesh, CFL and turbulence settings; the backends differ.
set -uo pipefail
MESH=${MESH:-work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet}
STEPS=${STEPS:-2}
BIN_REACT=${BIN_REACT:-./build-full/rfam_v20_tetra_react}
BIN_FLOW=${BIN_FLOW:-./build-full/rfam_v20_tetra}

echo "mesh=$MESH steps=$STEPS"
t0=$(date +%s)
$BIN_REACT --mesh "$MESH" --steps "$STEPS" --monitor-every 1 --blueprint \
  --thermo legacy --inlet-mode mass --p0 4e6 --T0 3400 --guard-T-max 4300 \
  --lox-pinj 6e6 --lox-T 3400 --lox-mdot 2.70 \
  --fuel-pinj 4.5e6 --fuel-T 3400 --fuel-mdot 1.01 --torch-mdot 0 \
  --reacting --mechanism twostep --chem-max-substeps 2 \
  --molecular-viscous --turbulence sst \
  --lusgs --cfl 0.5 --lusgs-dt-scale 0.05 > /tmp/bench_legacy.log 2>&1
t1=$(date +%s)
leg_steps=$(grep -c '^STEP' /tmp/bench_legacy.log || true)
echo "LEGACY   steps=$leg_steps wall=$((t1 - t0))s  per_step=$(( (t1 - t0) / (leg_steps > 0 ? leg_steps : 1) ))s"

t0=$(date +%s)
$BIN_FLOW --mesh "$MESH" --steps "$STEPS" --monitor-every 1 --outlet-exact \
  --thermo coolprop --inlet-mode mass --p0 4e6 --T0 300 \
  --lox-pinj 6e6 --lox-T 130 --lox-mdot 2.70 \
  --fuel-pinj 4.5e6 --fuel-T 450 --fuel-mdot 1.01 --torch-mdot 0 \
  --molecular-viscous --turbulence sst \
  --lusgs --cfl 0.5 --lusgs-dt-scale 0.02 > /tmp/bench_coolprop.log 2>&1
t1=$(date +%s)
cp_steps=$(grep -c '^STEP' /tmp/bench_coolprop.log || true)
echo "COOLPROP steps=$cp_steps wall=$((t1 - t0))s  per_step=$(( (t1 - t0) / (cp_steps > 0 ? cp_steps : 1) ))s"

echo "--- legacy STEP tail"
grep '^STEP' /tmp/bench_legacy.log | tail -1 | cut -c1-160
echo "--- coolprop STEP tail"
grep '^STEP' /tmp/bench_coolprop.log | tail -1 | cut -c1-160
