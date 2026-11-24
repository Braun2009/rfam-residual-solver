#!/usr/bin/env bash
# build_gerg_udrgm.sh - teqp GERG-2008 UDRGM for Fluent 2026 R1 (A406).
set -Eeuo pipefail

FLUENT_ROOT="${AWP_ROOT261:-/ansys_inc/v261}"
FLUENT_DIR="${FLUENT_ROOT}/fluent/fluent26.1.0"
FLUENT_SRC="${FLUENT_DIR}/src"
SRCDIR="${SRCDIR:-/root/vortex_fluent/udf_src}"
OUTDIR="${OUTDIR:-/root/vortex_fluent/udf_src/lnamd64}"
mkdir -p "$OUTDIR"

TEQP_INC="/root/teqp/include:/root/teqp/externals/nlohmann_json:/root/teqp/externals/Eigen:/root/teqp/externals/autodiff:/root/teqp/externals/stbrumme-hashing:/root/teqp/externals/json-schema-validator/src:/root/teqp/externals/REFPROP-interop/include"

includes=(main addon-wrapper io species pbns numerics sphysics storage mphase bc models
  material amg util mesh udf ht dx turbulence parallel etc ue acoustics dpm dbns)
include_args=()
for inc in "${includes[@]}"; do include_args+=("-I${FLUENT_SRC}/${inc}"); done
include_args+=(
  "-I${FLUENT_DIR}/lnamd64/3d_node"
  "-I${FLUENT_DIR}/cortex/src"
  "-I${FLUENT_DIR}/client/src"
  "-I${FLUENT_DIR}/tgrid/src"
  "-I${FLUENT_DIR}/multiport/src"
  "-I${FLUENT_ROOT}/fluent/include"
)

echo "== 1/3 g++17 teqp kernel =="
g++ -std=c++20 -O2 -fPIC -c "$SRCDIR/gerg_kernel.cpp" -o "$OUTDIR/gerg_kernel.o" \
  -I"$SRCDIR" $(echo "$TEQP_INC" | tr ':' '\n' | sed 's/^/-I/')

echo "== 2/3 gcc udrgm shell =="
gcc -std=gnu11 -O3 -fPIC -fopenmp -c "$SRCDIR/udrgm_gerg.c" -o "$OUTDIR/udrgm_gerg.o" \
  -DUDF_EXPORTING -D_lnamd64 -Dlinux -Dunix -DSYSV -D_REENTRANT \
  -I"$SRCDIR" "${include_args[@]}"

echo "== 3/3 link =="
g++ -shared -fopenmp -o "$OUTDIR/librfam_gerg_udrgm.so" \
  "$OUTDIR/gerg_kernel.o" "$OUTDIR/udrgm_gerg.o" -lm

echo "BUILD_OK: $OUTDIR/librfam_gerg_udrgm.so"
ls -la "$OUTDIR/librfam_gerg_udrgm.so"
