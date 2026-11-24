#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "usage: $0 [build-dir] [zero_d|one_d|uniform|shock|phase_foundation|phase|injector|long_3d]"
}

build_dir="${1:-build}"
through="${2:-shock}"

case "$through" in
    zero_d) target_rank=0 ;;
    one_d) target_rank=1 ;;
    uniform) target_rank=2 ;;
    shock) target_rank=3 ;;
    phase_foundation) target_rank=4 ;;
    phase) target_rank=5 ;;
    injector) target_rank=6 ;;
    long_3d) target_rank=7 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 64 ;;
esac

if [[ ! -f "$build_dir/CTestTestfile.cmake" ]]; then
    echo "VALIDATION_LADDER_ERROR: configure CMake build directory '$build_dir' first" >&2
    exit 66
fi

run_tier() {
    local rank="$1" name="$2" regex="$3"
    shift 3
    if (( target_rank < rank )); then
        return
    fi
    for required_test in "$@"; do
        if ! (cd "$build_dir" && ctest -N -R "^${required_test}$") \
            | grep -q "Test .*: ${required_test}$"; then
            echo "VALIDATION_LADDER_ERROR: missing test ${required_test}" >&2
            exit 67
        fi
    done
    echo "=== RFAM validation tier: $name ==="
    (cd "$build_dir" && ctest --output-on-failure -R "$regex")
    echo "VALIDATION_TIER_PASS: $name"
}

run_tier 0 zero_d \
    'test_pr_thermo_identities|test_pr_roots|test_pr_low_pressure_closure|test_cp_cv_positive|test_state_conversion_contract|test_real_fluid_chemistry_energy' \
    test_pr_thermo_identities test_pr_roots test_pr_low_pressure_closure \
    test_cp_cv_positive test_state_conversion_contract \
    test_real_fluid_chemistry_energy

run_tier 1 one_d \
    'test_phase7_1d|test_diffusion_real_cv' \
    test_phase7_1d test_diffusion_real_cv

run_tier 2 uniform_flow \
    'test_openmp_consistency|test_sector_freestream|test_sector_reacting_freestream' \
    test_openmp_consistency test_sector_freestream \
    test_sector_reacting_freestream

run_tier 3 shock \
    'test_sod_tube' test_sod_tube

run_tier 4 phase_foundation \
    'test_pr_tp_flash|test_pr_phase_caloric|test_pr_ph_flash|test_pr_rho_e_flash|test_pr_closed_volume_heating|test_pr_constant_pressure_heating|test_pr_two_phase_mixing_contract|test_pr_vapor_chemistry_reflash_contract|test_pr_equilibrium_domain|test_acoustic_state_routing|test_transport_state_routing|test_tabulated_ideal_caloric|test_two_phase_state_contract|test_model_capability_gate|test_restart_manifest|test_restart_manifest_codec|test_restart_generation_store|test_restart_atomic_publication_contract|test_restart_publication|test_restart_generation_store_sector_clock_preflight|test_sector_restart_owner_codec|test_sector_restart_publication_stagers' \
    test_pr_tp_flash test_pr_phase_caloric test_pr_ph_flash test_pr_rho_e_flash \
    test_pr_closed_volume_heating test_pr_constant_pressure_heating \
    test_pr_two_phase_mixing_contract \
    test_pr_vapor_chemistry_reflash_contract \
    test_pr_equilibrium_domain test_acoustic_state_routing \
    test_transport_state_routing \
    test_tabulated_ideal_caloric test_two_phase_state_contract \
    test_model_capability_gate \
    test_restart_manifest test_restart_manifest_codec \
    test_restart_generation_store test_restart_atomic_publication_contract \
    test_restart_publication \
    test_restart_generation_store_sector_clock_preflight \
    test_sector_restart_owner_codec test_sector_restart_publication_stagers

if (( target_rank >= 5 )); then
    cat >&2 <<'EOF'
VALIDATION_LADDER_BLOCKED: phase
    RFAM now routes its verified HEM sound speed through Q into first-order
    flux/CFL and spectral-radius consumers, but phase-aware face reconstruction,
    characteristic boundaries, dense cryogenic transport, a full flow Jacobian,
    and phase-aware chemistry remain incomplete. Foundation tests do not by
    themselves authorize production CFD.
Injector and long-3D tiers are intentionally unreachable until dedicated
phase-boundary tests exist and pass.
EOF
    exit 2
fi

echo "VALIDATION_LADDER_PASS: through $through"
