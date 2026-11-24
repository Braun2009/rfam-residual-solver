#!/usr/bin/env bash
# Run the Vortex 3D cold-air SST commissioning journal natively in Fluent.
# Usage: run_vortex_3d_cold_air_sst_journal.sh [work-dir] [ranks]
set -euo pipefail

WORK_DIR="${1:-/root/vortex_fluent/gmsh_fixed_v2/solver_cold_air_sst_journal_20260815_r4_36phys}"
RANKS="${2:-36}"
JOURNAL="$WORK_DIR/vortex_3d_cold_air_sst_r4.jou"

if [[ ! -f "$JOURNAL" ]]; then
    echo "journal not found: $JOURNAL" >&2
    exit 2
fi
if [[ ! "$RANKS" =~ ^[1-9][0-9]*$ ]]; then
    echo "rank count must be a positive integer: $RANKS" >&2
    exit 2
fi

export ANSYSLMD_LICENSE_FILE="${ANSYSLMD_LICENSE_FILE:-1055@localhost}"
export ANSYSCL_FNP_PATH="${ANSYSCL_FNP_PATH:-$ANSYSLMD_LICENSE_FILE}"

cd "$WORK_DIR"
exec /ansys_inc/v261/fluent/bin/fluent 3ddp -g "-t${RANKS}" -i "$JOURNAL"
