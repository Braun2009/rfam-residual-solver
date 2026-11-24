#!/bin/bash
cd /root/realFluidAllMach
export RFAM_SOURCE_DIR=/root/realFluidAllMach RFAM_TEQP_DATA=/root/teqp/dev/fluiddata OMP_NUM_THREADS=36
./build-a406/rfam_v20_sector_smoke --sandiego --reacting   --allow-single-phase-reacting-diagnostic --les   --dual-time --dt-mult 40 --steps 800   --restart cardd_dualtime_iter.chk   --checkpoint cardd_dualtime2_iter.chk --checkpoint-every 100   --monitor-every 100   --field-csv cardd_dualtime2_field.csv --wall-csv cardd_dualtime2_wall.csv   > cardd_dualtime2_run.log 2>&1
echo DUALTIME_LONGRUN2_DONE exit=$? >> cardd_dualtime2_run.log
