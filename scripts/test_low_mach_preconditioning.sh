#!/bin/bash
# Quick test: v3 injector 0.6 kg/s with low-Mach preconditioning

set -e

RFAM=./build/solver_3d_tetra
MESH=work/path4_correct_domain/injector_r3.rfamtet

if [ ! -f "$RFAM" ]; then
    echo "ERROR: solver not built at $RFAM"
    echo "Build with: cd build && cmake .. && make -j"
    exit 1
fi

if [ ! -f "$MESH" ]; then
    echo "ERROR: mesh not found at $MESH"
    exit 1
fi

echo "=== Test 1: Baseline (no preconditioning) ==="
$RFAM \
    --mesh $MESH \
    --steps 200 \
    --cfl 1.2 \
    --lusgs \
    --lusgs-dt-scale 2.0 \
    --lusgs-relax 0.7 \
    --lox-mass-flow 0.6 \
    --lox-pinj 4.2e6 \
    --lox-temperature 90.0 \
    --outlet-pressure 4.0e6 \
    --guard-p-max 8e6 \
    --monitor-every 20 \
    --dump baseline_06kg.vtu \
    --dump-every 200 \
    > baseline_06kg.log 2>&1

echo "Baseline complete. Check baseline_06kg.log for convergence."

echo ""
echo "=== Test 2: With low-Mach preconditioning ==="
$RFAM \
    --mesh $MESH \
    --steps 200 \
    --cfl 1.2 \
    --lusgs \
    --lusgs-dt-scale 2.0 \
    --lusgs-relax 0.7 \
    --low-mach-preconditioning \
    --low-mach-cutoff 0.05 \
    --lox-mass-flow 0.6 \
    --lox-pinj 4.2e6 \
    --lox-temperature 90.0 \
    --outlet-pressure 4.0e6 \
    --guard-p-max 8e6 \
    --monitor-every 20 \
    --dump precond_06kg.vtu \
    --dump-every 200 \
    > precond_06kg.log 2>&1

echo "Preconditioning test complete. Check precond_06kg.log for convergence."

echo ""
echo "=== Comparison ==="
echo "Baseline:"
grep -E "step.*residual|CONVERGED|step 200" baseline_06kg.log | tail -5
echo ""
echo "Preconditioned:"
grep -E "step.*residual|CONVERGED|step 200" precond_06kg.log | tail -5

echo ""
echo "Expected: Preconditioning should converge in ~50-100 steps vs 200+ for baseline"
