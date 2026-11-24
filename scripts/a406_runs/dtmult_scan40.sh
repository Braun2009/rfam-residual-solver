#!/bin/bash
cd /root/realFluidAllMach
export RFAM_SOURCE_DIR=/root/realFluidAllMach RFAM_TEQP_DATA=/root/teqp/dev/fluiddata OMP_NUM_THREADS=36
./build-a406/rfam_v20_sector_smoke --sandiego --reacting   --allow-single-phase-reacting-diagnostic --les   --dual-time --dt-mult 40 --steps 50   --restart cardd_sandiego_iter5700.chk   --monitor-every 25   --field-csv dtmult_scan_40_field.csv --wall-csv dtmult_scan_40_wall.csv   > dtmult_scan_40.log 2>&1
echo SCAN40_DONE exit=$? >> dtmult_scan_summary.log
