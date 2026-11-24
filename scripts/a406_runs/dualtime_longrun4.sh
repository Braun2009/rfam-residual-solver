#!/bin/bash
cd /root/realFluidAllMach
export RFAM_SOURCE_DIR=/root/realFluidAllMach RFAM_TEQP_DATA=/root/teqp/dev/fluiddata OMP_NUM_THREADS=36
./build-a406/rfam_v20_sector_smoke --sandiego --reacting   --allow-single-phase-reacting-diagnostic --les   --dual-time --dt-mult 20 --steps 200   --restart cardd_dualtime3_iter.chk   --checkpoint cardd_dualtime4_iter.chk --checkpoint-every 50   --monitor-every 50   --field-csv cardd_dualtime4_field.csv --wall-csv cardd_dualtime4_wall.csv   > cardd_dualtime4_run.log 2>&1
echo DUALTIME_LONGRUN4_DONE exit=$? >> cardd_dualtime4_run.log
