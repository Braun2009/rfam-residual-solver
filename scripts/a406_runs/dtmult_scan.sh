#!/bin/bash
cd /root/realFluidAllMach
export RFAM_SOURCE_DIR=/root/realFluidAllMach RFAM_TEQP_DATA=/root/teqp/dev/fluiddata OMP_NUM_THREADS=36
: > dtmult_scan_summary.log
for mult in 5 10 20; do
  echo "=== dtMult=$mult start $(date +%H:%M:%S) ===" >> dtmult_scan_summary.log
  ./build-a406/rfam_v20_sector_smoke --sandiego --reacting     --allow-single-phase-reacting-diagnostic --les     --dual-time --dt-mult $mult --steps 50     --restart cardd_sandiego_iter5700.chk     --monitor-every 25     --field-csv dtmult_scan_${mult}_field.csv --wall-csv dtmult_scan_${mult}_wall.csv     > dtmult_scan_${mult}.log 2>&1
  echo "=== dtMult=$mult exit=$? end $(date +%H:%M:%S) ===" >> dtmult_scan_summary.log
done
echo ALL_SCAN_DONE >> dtmult_scan_summary.log
