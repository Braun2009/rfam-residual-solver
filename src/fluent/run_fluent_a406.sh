#!/usr/bin/env bash
set -euo pipefail

export ANSYSLMD_LICENSE_FILE="${ANSYSLMD_LICENSE_FILE:-1055@localhost}"
export ANSYSCL_FNP_PATH="${ANSYSCL_FNP_PATH:-$ANSYSLMD_LICENSE_FILE}"

FLUENT_BIN="${FLUENT_BIN:-/ansys_inc/v261/fluent/bin/fluent}"
exec "$FLUENT_BIN" "$@"
