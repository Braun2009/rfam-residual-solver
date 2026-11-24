#!/bin/bash
# run_convergence.sh - MUSCL grid convergence + transcritical stability study
# Run on server: bash run_convergence.sh
set -e
cd "$(dirname "$0")/.."
export SUN="${SUN:-/opt/sundials}"
INCALL="-Isrc -Isrc/thermo -Isrc/thermo/EquationOfState -Isrc/thermo/MixingRules -Isrc/thermo/RootSelection -Isrc/variables -Isrc/numerics/flux -Isrc/numerics/preconditioning -Isrc/numerics/reconstruction -Isrc/mesh -Isrc/combustion/FPV -Isrc/framework -Isrc/transport -Isrc/turbulence"
CORE="src/thermo/RealFluidState.C src/thermo/EquilibriumThermoModel.C src/thermo/EquationOfState/PengRobinsonEOS.C src/thermo/MixingRules/VdWOneFluidMixing.C src/thermo/MixingRules/SpeciesCriticalDB.C src/thermo/MixingRules/BinaryInteractionDB.C src/thermo/RootSelection/CubicRoots.C src/variables/ConservedVars.C src/variables/CompleteCellState.C src/numerics/flux/AUSMPlusUpFlux.C src/numerics/preconditioning/ExperimentalCflAcousticScaler.C src/combustion/FPV/FPVTable.C"

echo "============================================="
echo " RFAM Server Compute Tasks"
echo "============================================="

echo ""
echo "--- Task 3: MUSCL Sod (default 200 cells) ---"
g++ -O2 -std=c++17 -I. $INCALL tests/test_muscl_sod.cpp $CORE -o /tmp/rfam_muscl
/tmp/rfam_muscl

echo ""
echo "--- Task 1: Implicit stability (transcritical Tj=120K) ---"
g++ -O2 -std=c++17 -I. $INCALL tests/test_implicit_stability.cpp $CORE -o /tmp/rfam_impl
/tmp/rfam_impl

echo ""
echo "--- Task 2: Free jet k-omega (80x40, check spreading) ---"
g++ -O2 -std=c++17 -I. $INCALL tests/test_kw_freejet.cpp $CORE -o /tmp/rfam_kw
/tmp/rfam_kw

echo ""
echo "ALL_SERVER_TASKS_DONE"
