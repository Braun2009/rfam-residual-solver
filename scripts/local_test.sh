#!/bin/bash
#------------------------------------------------------------------------------
# RFAM 本地回归套件 (WSL g++ ; 纯 std C++, 无需 OpenFOAM/Eigen/GTest)
#   Phase 6 (CVODE) 需本地 SUNDIALS: F:\opt\sundials  (见 memory/reference/rfam-local-build-env)
# 用法: wsl -e bash -lc 'bash /mnt/f/realFluidAllMach/scripts/local_test.sh'
#------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/.." || exit 1
SUN=/mnt/f/opt/sundials
PASS=0; FAIL=0; MATCHED=0
STAGE5_QWALL_READY=0

# Headers use short project-local include names. Keep the shared search path in
# one place so newly introduced operators do not require editing every command.
export CPLUS_INCLUDE_PATH="$PWD/src/thermo:$PWD/src/numerics/reconstruction:${CPLUS_INCLUDE_PATH:-}"

INCT="-Isrc -Isrc/thermo -Isrc/thermo/EquationOfState -Isrc/thermo/MixingRules -Isrc/thermo/RootSelection"
THERMO="src/thermo/RealFluidState.C src/thermo/ThermoDataProvenance.C src/thermo/EquationOfState/PengRobinsonEOS.C src/thermo/MixingRules/VdWOneFluidMixing.C src/thermo/MixingRules/SpeciesCriticalDB.C src/thermo/MixingRules/BinaryInteractionDB.C src/thermo/RootSelection/CubicRoots.C"
PHASE="src/thermo/PhaseEquilibrium/PRMixtureKernel.C src/thermo/PhaseEquilibrium/PRTPFlash.C src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.C src/thermo/PhaseEquilibrium/PRPhaseCaloric.C src/thermo/PhaseEquilibrium/PRPHFlash.C src/thermo/PhaseEquilibrium/PRRhoEFlash.C src/thermo/PhaseEquilibrium/PRPhaseSoundSpeed.C"
EQUILIBRIUM_MODEL="$PHASE src/thermo/PhaseEquilibrium/PREquilibriumThermoModel.C src/thermo/PhaseEquilibrium/PRPhaseAuthorityIdentity.C"
CORE="$THERMO src/thermo/EquilibriumThermoModel.C src/variables/ConservedVars.C src/variables/CompleteCellState.C src/numerics/flux/AUSMPlusUpFlux.C src/numerics/preconditioning/ExperimentalCflAcousticScaler.C"
INCALL="$INCT -Isrc/variables -Isrc/numerics/flux -Isrc/numerics/preconditioning -Isrc/mesh -Isrc/transport -Isrc/turbulence -Isrc/combustion/FPV"
CHEM="src/combustion/FiniteRate/TwoStepCO.C"
SUNLIBS="-L$SUN/lib -Wl,-rpath,$SUN/lib -lsundials_cvode -lsundials_sunlinsoldense -lsundials_sunmatrixdense -lsundials_sunnonlinsolnewton -lsundials_nvecserial -lsundials_core"
CANTERA_GENERATOR="scripts/cantera_reference.py"
CANTERA_REFERENCE="tests/reference/cantera_gri30_propane_oxygen_v1.txt"
CANTERA_MECHANISM="tests/reference/cantera_gri30_v3_0.yaml"
CANTERA_GENERATOR_SHA=$(sha256sum "$CANTERA_GENERATOR" | awk '{print $1}')
CANTERA_REFERENCE_SHA=$(sha256sum "$CANTERA_REFERENCE" | awk '{print $1}')
CANTERA_MECHANISM_SHA=$(sha256sum "$CANTERA_MECHANISM" | awk '{print $1}')
CANTERA_CPP_DEFINES="-DRFAM_CANTERA_GENERATOR_SHA256=\\\"$CANTERA_GENERATOR_SHA\\\" -DRFAM_CANTERA_REFERENCE_SHA256=\\\"$CANTERA_REFERENCE_SHA\\\" -DRFAM_CANTERA_MECHANISM_SHA256=\\\"$CANTERA_MECHANISM_SHA\\\""

selected() {
  local name="$1"
  [ -z "${RFAM_TEST_FILTER:-}" ] || [[ "$name" =~ ${RFAM_TEST_FILTER} ]]
}

run() { # $1 name  $2 compile-cmd (without -o)  $3 expected  $4 run args
  local name="$1" cc="$2" exp="$3" args="${4:-}" bin="/tmp/rfam_$1"
  local -a run_argv=()
  if [ -n "$args" ]; then read -r -a run_argv <<< "$args"; fi
  selected "$name" || return
  MATCHED=$((MATCHED+1))
  if [ "$name" = "stage5_qwall" ]; then
    rm -f /tmp/rfam_qwall.json
    STAGE5_QWALL_READY=0
  fi
  if eval "$cc -o $bin" 2>"/tmp/rfam_err_$1"; then
    local out run_status
    out=$(timeout "${RFAM_TEST_TIMEOUT:-90}" "$bin" "${run_argv[@]}" 2>&1)
    run_status=$?
    if [ "$run_status" -eq 0 ] && echo "$out" | grep -qF "$exp"; then
      echo "  [PASS] $name"; PASS=$((PASS+1))
      if [ "$name" = "stage5_qwall" ]; then
        STAGE5_QWALL_READY=1
      fi
    elif [ "$run_status" -ne 0 ]; then
      echo "  [FAIL-run] $name (exit=$run_status)"
      [ -z "$out" ] || printf '%s\n' "$out"
      FAIL=$((FAIL+1))
    else
      echo "  [FAIL-run] $name (missing: '$exp')"
      [ -z "$out" ] || printf '%s\n' "$out"
      FAIL=$((FAIL+1))
    fi
  else echo "  [FAIL-compile] $name  (see /tmp/rfam_err_$1)"; FAIL=$((FAIL+1)); fi
}

echo "=============================================="
echo " RFAM Local Regression Suite"
echo "=============================================="

run stage0_thermo  "g++ -O2 -std=c++17 $INCT tests/test_cv_validation.cpp $THERMO"                                "Stage 0 COMPLETE"
run nasa7_thermo   "g++ -O2 -std=c++17 -Isrc -Isrc/thermo tests/test_nasa7_thermo.cpp"                             "NASA7_THERMO_PASS"
run pr_mixing      "g++ -O2 -std=c++17 $INCT tests/test_pr_mixing.cpp $THERMO"                                     "PR_MIXING_PASS"
run pr_state_scan  "g++ -O2 -std=c++17 $INCT tests/test_pr_state_scan.cpp $THERMO"                               "PR_STATE_SCAN_PASS"
run pr_tp_flash    "g++ -O2 -std=c++17 $INCT tests/thermo/test_pr_tp_flash.cpp $THERMO $PHASE"              "PR_TP_FLASH_PASS"
run pr_phase_caloric "g++ -O2 -std=c++17 $INCT tests/thermo/test_pr_phase_caloric.cpp $THERMO $PHASE"       "PR_PHASE_CALORIC_PASS"
run pr_ph_flash    "g++ -O2 -std=c++17 $INCT tests/thermo/test_pr_ph_flash.cpp $THERMO $PHASE"              "PR_PH_FLASH_PASS"
run pr_rho_e_flash "g++ -O2 -std=c++17 $INCT tests/thermo/test_pr_rho_e_flash.cpp $THERMO $PHASE"          "PR_RHO_E_FLASH_PASS"
run pr_equilibrium_domain "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework -Isrc/transport -Isrc/turbulence $INCALL tests/thermo/test_pr_equilibrium_domain.cpp $CORE $EQUILIBRIUM_MODEL" "PR_EQUILIBRIUM_DOMAIN_PASS"
run pr_equilibrium_short_step "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework -Isrc/transport -Isrc/turbulence $INCALL tests/thermo/test_pr_equilibrium_short_step.cpp $CORE $EQUILIBRIUM_MODEL" "PR_EQUILIBRIUM_SHORT_STEP_PASS"
run acoustic_state_routing "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework $INCALL tests/thermo/test_acoustic_state_routing.cpp $CORE $EQUILIBRIUM_MODEL" "ACOUSTIC_STATE_ROUTING_PASS"
run transport_state_routing "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework -Isrc/mesh -Isrc/transport -Isrc/turbulence $INCALL tests/thermo/test_transport_state_routing.cpp $CORE $EQUILIBRIUM_MODEL src/transport/CryogenicTransportTable.C" "TRANSPORT_STATE_ROUTING_PASS"
run transport_no_alloc "g++ -O2 -std=c++17 -I. -Isrc/framework -Isrc/transport $INCALL tests/thermo/test_transport_no_alloc.cpp $CORE" "TRANSPORT_NO_ALLOC_PASS"
run state_conversion_contract "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework $INCALL tests/thermo/test_state_conversion_contract.cpp $CORE" "STATE_CONVERSION_CONTRACT_PASS"
run tabulated_ideal_caloric "g++ -O2 -std=c++17 $INCT tests/thermo/test_tabulated_ideal_caloric.cpp src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.C" "TABULATED_IDEAL_CALORIC_PASS"
run two_phase_state_contract "g++ -O2 -std=c++17 $INCT tests/thermo/test_two_phase_state_contract.cpp"      "TWO_PHASE_STATE_ALGEBRA_PASS"
run model_capability_gate "g++ -O2 -std=c++17 -Isrc/framework tests/thermo/test_model_capability_gate.cpp"  "MODEL_CAPABILITY_GATE_BLOCKED"
run cantera_compare "g++ -O2 -std=c++17 $CANTERA_CPP_DEFINES $INCT -Isrc/variables -Isrc/combustion -Isrc/combustion/FiniteRate -I$SUN/include tests/test_cantera_compare.cpp src/thermo/ThermoDataProvenance.C $THERMO $CHEM $SUNLIBS" "CANTERA_KINETICS_GATE_BLOCKED" "$CANTERA_GENERATOR $CANTERA_REFERENCE $CANTERA_MECHANISM"
run solver2d        "g++ -O2 -std=c++17 -I. -Isrc/framework $INCALL tests/test_solver2d.cpp $CORE" "SOLVER2D_PASS"
run domain3d_lifecycle "g++ -O2 -std=c++17 -I. -Isrc/framework $INCALL tests/test_domain3d_lifecycle.cpp $CORE" "DOMAIN3D_LIFECYCLE_PASS"
run step_transaction "g++ -O2 -std=c++17 -I. -Isrc/framework -Isrc/transport $INCALL tests/test_step_transaction.cpp $CORE" "STEP_TRANSACTION_PASS"
run openmp_consistency "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework $INCALL tests/test_openmp_consistency.cpp $CORE" "OPENMP_CONSISTENCY_PASS"
run isothermal_projection "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework $INCALL tests/test_isothermal_projection.cpp $CORE" "ISOTHERMAL_PROJECTION_PASS"
run sector_metrics "g++ -O2 -std=c++17 -Isrc -Isrc/mesh tests/test_sector_metrics.cpp" "SECTOR_METRICS_PASS"
run sector_freestream "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework -Isrc/mesh $INCALL tests/test_sector_freestream.cpp $CORE" "SECTOR_FREESTREAM_PASS"
run sector_turbulence "g++ -O2 -std=c++17 -fopenmp -I. -Isrc/framework -Isrc/mesh -Isrc/turbulence -Isrc/transport -Isrc/thermo $INCALL tests/test_sector_turbulence.cpp $CORE src/transport/CryogenicTransportTable.C" "SECTOR_TURBULENCE_PASS"
run reacting2d      "g++ -O2 -std=c++17 -I. -Isrc/framework -Isrc/combustion -Isrc/combustion/FiniteRate -Isrc/transport -I$SUN/include $INCALL tests/test_2d_reacting.cpp $CORE $CHEM $SUNLIBS" "REACTING2D_PASS"
run cht2d           "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/cht -Isrc/transport tests/test_cht2d.cpp $CORE" "CHT2D_PASS"
run cht2d_conserv   "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/cht -Isrc/transport tests/test_cht2d_conservation.cpp $CORE" "CHT2D_CONSERVATION_PASS"
run cht_relax_cons  "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/cht -Isrc/transport tests/test_cht_relaxation_vs_conservative.cpp $CORE" "CHT_RELAX_VS_CONS_PASS"
run turbulence      "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/turbulence tests/test_turbulence.cpp $CORE" "TURBULENCE_PASS"
run sod_tube        "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework tests/test_sod_tube.cpp $CORE" "SOD_TUBE_PASS"
run vtk_output      "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/turbulence tests/test_vtk_output.cpp $CORE" "VTK_OUTPUT_PASS"
run mech_reader     "g++ -O2 -std=c++17 -I. $INCALL -Isrc/combustion -Isrc/combustion/FiniteRate tests/test_mechanism_reader.cpp $CHEM" "TWOSTEP_CO_CONTRACT_PASS"
run mayer_jet       "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework tests/test_mayer_jet.cpp $CORE" "MAYER_JET_PASS"
run stage5_qwall    "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/cht -Isrc/transport tests/test_stage5_qwall.cpp $CORE" "STAGE5_QWALL_PASS"
run tetra_low_mach_preconditioning "g++ -O2 -std=c++17 tests/test_tetra_low_mach_preconditioning.cpp" "TETRA_LOW_MACH_PRECONDITIONING_PASS"
run kw_freejet      "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/turbulence tests/test_kw_freejet.cpp $CORE" "KW_FREEJET_PASS"
run vortex_coldflow  "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework tests/test_vortex_coldflow.cpp $CORE" "VORTEX_COLDFLOW_PASS"
run vortex_cht       "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/cht tests/test_vortex_cht.cpp $CORE" "VORTEX_CHT_PASS"
run implicit_stab    "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework tests/test_implicit_stability.cpp $CORE" "IMPLICIT_STABILITY_PASS"
run muscl_sod        "g++ -O2 -std=c++17 -I. $INCALL -Isrc/numerics/reconstruction -Isrc/framework tests/test_muscl_sod.cpp $CORE" "MUSCL_SOD_PASS"
run nonreflecting   "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework tests/test_nonreflecting.cpp $CORE" "NONREFLECTING_PASS"
run diffusion      "g++ -O2 -std=c++17 -Isrc -Isrc/thermo -Isrc/framework -Isrc/transport $INCALL tests/test_diffusion.cpp $CORE" "DIFFUSION_PASS"
run experimental_cfl_scaler "g++ -O2 -std=c++17 -I. -Isrc/framework $INCALL tests/test_preconditioner.cpp $CORE" "EXPERIMENTAL_CFL_SCALER_CONTRACT_PASS"
run stage3_fpv     "g++ -O2 -std=c++17 -I. $INCALL -Isrc/framework -Isrc/combustion/FPV $CORE src/combustion/FPV/FPVTable.C tests/test_fpv_table.cpp" "Test Complete"
run phase5_rosen   "g++ -O2 -std=c++17 -Isrc -Isrc/combustion/FiniteRate tests/test_phase5_rosenbrock.cpp $CHEM"  "Rosenbrock solver stable"
run phase6_cvode   "g++ -O2 -std=c++17 -Isrc -Isrc/combustion/FiniteRate -I$SUN/include tests/test_phase6_cvode.cpp $CHEM $SUNLIBS" "all three conservation tests PASS"
run stage4_solid   "g++ -O2 -std=c++17 -Isrc tests/test_stage4_solid.cpp"                                         "SolidEnergyEqn VALIDATED"
run stage4_cht     "g++ -O2 -std=c++17 -Isrc tests/test_stage4_cht.cpp"                                           "CHTCoupler VALIDATED"
run mvp1d_frame    "g++ -O2 -std=c++17 -I. -Isrc/framework $INCALL apps/solver_mvp_1d.cpp $CORE src/combustion/FPV/FPVTable.C" "MVP_REPRODUCES_solver1d_fpv"
run frame_cfl_scaler "g++ -O2 -std=c++17 -I. -Isrc/framework $INCALL apps/solver_frame_precond.cpp $CORE src/combustion/FPV/FPVTable.C" "FRAME_CFL_SCALER_DIAGNOSTIC_PASS"
run stage15_species "g++ -O2 -std=c++17 -Isrc -Isrc/thermo -Isrc/combustion -Isrc/combustion/FiniteRate tests/test_stage15_species.cpp src/combustion/FiniteRate/TwoStepCO.C" "STAGE1.5_SPECIES_PASS"
run phase7_chem0d  "g++ -O2 -std=c++17 -I. -Isrc/framework -Isrc/combustion -Isrc/combustion/FiniteRate -I$SUN/include $INCALL tests/test_phase7_chem0d.cpp $CORE $CHEM $SUNLIBS" "PHASE7_CHEM0D_PASS"
run phase7_operator "g++ -O2 -std=c++17 -Isrc -Isrc/framework -Isrc/combustion -Isrc/combustion/FPV -Isrc/combustion/FiniteRate -I$SUN/include $INCALL tests/test_phase7_operator_compile.cpp $CORE src/combustion/FPV/FPVTable.C src/combustion/FiniteRate/TwoStepCO.C $SUNLIBS" "PHASE7_OPERATOR_COMPILE_OK"
# Phase 7 Step 2: species-aware EOS reback (setSpecies) + 1D reacting flow
run setspecies_mw  "g++ -O2 -std=c++17 $INCT tests/test_setspecies_mw.cpp $THERMO" "SETSPECIES_MW_PASS"
run phase7_1d      "g++ -O2 -std=c++17 -fopenmp -Isrc -Isrc/thermo -Isrc/framework -Isrc/combustion -Isrc/combustion/FiniteRate -I$SUN/include $INCALL tests/test_phase7_1d.cpp $CORE src/combustion/FiniteRate/TwoStepCO.C $SUNLIBS" "PHASE7_1D_PASS"

# Stage 5 RFAM->PinFin JSON bridge. It may consume only the qwall artifact
# produced successfully by this invocation, never a stale /tmp file.
if selected pinfin_bridge; then
  MATCHED=$((MATCHED+1))
  if [ "$STAGE5_QWALL_READY" -eq 1 ] && command -v python3 >/dev/null 2>&1; then
    if python3 tests/test_pinfin_bridge.py /tmp/rfam_qwall.json 2>/tmp/rfam_err_pinfin_bridge | grep -qF "STAGE5_BRIDGE_PASS"; then
      echo "  [PASS] pinfin_bridge"; PASS=$((PASS+1))
    else
      echo "  [FAIL-run] pinfin_bridge (see /tmp/rfam_err_pinfin_bridge)"; FAIL=$((FAIL+1))
    fi
  else
    echo "  [SKIP] pinfin_bridge (current stage5_qwall result unavailable)"
  fi
fi

# Vortex Bartz 1D heat flux validation (Python; validates q_throat~42 MW/m2)
if selected bartz_vortex; then
  MATCHED=$((MATCHED+1))
  if command -v python3 >/dev/null 2>&1; then
    if python3 tools/bartz_vortex.py 2>/tmp/rfam_err_bartz | grep -qF "BARTZ_VORTEX_PASS"; then
      echo "  [PASS] bartz_vortex"; PASS=$((PASS+1))
    else
      echo "  [FAIL-run] bartz_vortex (see /tmp/rfam_err_bartz)"; FAIL=$((FAIL+1))
    fi
  else
    echo "  [SKIP] bartz_vortex (python3 unavailable)"
  fi
fi

# SUS316 coolant reverse-solve + pin-fin/film assessment
if selected coolant_solve; then
  MATCHED=$((MATCHED+1))
  if command -v python3 >/dev/null 2>&1; then
    if python3 tools/coolant_reverse_solve.py 2>/tmp/rfam_err_coolant | grep -qF "COOLANT_SOLVE_PASS"; then
      echo "  [PASS] coolant_solve"; PASS=$((PASS+1))
    else
      echo "  [FAIL-run] coolant_solve (see /tmp/rfam_err_coolant)"; FAIL=$((FAIL+1))
    fi
  else
    echo "  [SKIP] coolant_solve (python3 unavailable)"
  fi
fi

if [ -n "${RFAM_TEST_FILTER:-}" ] && [ "$MATCHED" -eq 0 ]; then
  echo "  [FAIL-filter] no test matched RFAM_TEST_FILTER='$RFAM_TEST_FILTER'"
  FAIL=$((FAIL+1))
fi

echo "----------------------------------------------"
echo " RESULT: $PASS passed, $FAIL failed"
echo "=============================================="
exit "$FAIL"
