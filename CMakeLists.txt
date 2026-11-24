cmake_minimum_required(VERSION 3.16)
project(realFluidAllMach VERSION 0.1.0 LANGUAGES CXX)

# teqp GERG-2008 multifluid C++ integration (requires C++20 <concepts>,
# g++ 10+; A406 has g++ 11).  RFAM_TEQP_DATA points at the teqp fluiddata
# root at runtime (default /root/teqp/dev/fluiddata).
option(RFAM_HAVE_TEQP "teqp multifluid GERG-2008 mixture thermodynamics" OFF)
if(RFAM_HAVE_TEQP)
    add_compile_definitions(RFAM_HAVE_TEQP)
    include_directories(${RFAM_TEQP_INCLUDE})
    # teqp derivs.hpp uses C++20 concepts; -fconcepts in C++17 mode is
    # a compatible extension on g++ 10+ (A406 g++ 11).
    add_compile_options(-fconcepts)
endif()

set(RFAM_GIT_COMMIT "unknown" CACHE STRING "Source git commit embedded in checkpoints")
if(RFAM_GIT_COMMIT STREQUAL "unknown")
    execute_process(COMMAND git rev-parse --git-path HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE RFAM_GIT_HEAD_FILE OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    execute_process(COMMAND git symbolic-ref -q HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE RFAM_GIT_HEAD_REF OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(RFAM_GIT_HEAD_FILE)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${RFAM_GIT_HEAD_FILE}")
    endif()
    if(RFAM_GIT_HEAD_REF)
        execute_process(COMMAND git rev-parse --git-path ${RFAM_GIT_HEAD_REF}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE RFAM_GIT_REF_FILE OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(RFAM_GIT_REF_FILE)
            set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/${RFAM_GIT_REF_FILE}")
        endif()
    endif()
    execute_process(COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE RFAM_DETECTED_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(RFAM_DETECTED_GIT_COMMIT)
        set(RFAM_GIT_COMMIT ${RFAM_DETECTED_GIT_COMMIT})
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(RFAM_HAVE_COOLPROP "Enable native CoolProp 7.2.0 HEOS thermodynamics" OFF)
set(RFAM_COOLPROP_INCLUDE_DIR "" CACHE PATH "CoolProp include directory")
set(RFAM_COOLPROP_LIBRARY "" CACHE FILEPATH "Native CoolProp library")
set(RFAM_COOLPROP_PROVENANCE "unrecorded" CACHE STRING "Verified native CoolProp source identity")
option(RFAM_ENABLE_OPENMP "Enable OpenMP CPU parallelism" ON)
option(RFAM_ENABLE_CHEMISTRY "Build finite-rate SUNDIALS chemistry support" ON)
set(SUNDIALS_ROOT "" CACHE PATH "SUNDIALS installation prefix")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(rfam_core
    src/thermo/ThermoDataProvenance.C
    src/thermo/RealFluidState.C
    src/thermo/EquilibriumThermoModel.C
    src/thermo/CoolPropThermoModel.C
    src/thermo/EquationOfState/PengRobinsonEOS.C
    src/thermo/MixingRules/VdWOneFluidMixing.C
    src/thermo/MixingRules/SpeciesCriticalDB.C
    src/thermo/MixingRules/BinaryInteractionDB.C
    src/thermo/RootSelection/CubicRoots.C
    src/thermo/PhaseEquilibrium/PRMixtureKernel.C
    src/thermo/PhaseEquilibrium/PRTPFlash.C
    src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.C
    src/thermo/PhaseEquilibrium/PRPhaseCaloric.C
    src/thermo/PhaseEquilibrium/PRPhaseAuthorityIdentity.C
    src/thermo/PhaseEquilibrium/PRPhaseRhoERecovery.C
    src/thermo/PhaseEquilibrium/PRPhaseRootSignature.C
    src/thermo/PhaseEquilibrium/PRPhaseDisplacementAuthority.C
    src/thermo/PhaseEquilibrium/PRPhaseThermoOwnerAdapter.C
    src/thermo/PhaseEquilibrium/PRResidualHelmholtzAuthority.C
    src/thermo/PhaseEquilibrium/PRPHFlash.C
    src/thermo/PhaseEquilibrium/PRRhoEFlash.C
    src/thermo/PhaseEquilibrium/PRPhaseSoundSpeed.C
    src/thermo/PhaseEquilibrium/PREquilibriumThermoModel.C
    src/transport/CryogenicTransportTable.C
    src/combustion/CondensedCarbonModel.C
    src/combustion/CondensedCarbonOperator.C
    src/combustion/CondensedCarbonCheckpoint.C
    src/combustion/CondensedCarbonOutput.C
    src/variables/ConservedVars.C
    src/variables/CompleteCellState.C
    src/framework/RestartManifest.C
    src/framework/RestartManifestCodec.C
    src/framework/RestartGenerationStore.C
    src/framework/RestartPublication.C
    src/framework/RestartLiveRoot.C
    src/framework/SectorRestartPublicationStagers.C
    src/numerics/flux/AUSMPlusUpFlux.C
    src/numerics/preconditioning/ExperimentalCflAcousticScaler.C
    src/mesh/TetraMesh.C
    src/mesh/TetraFiniteVolume.C
    src/framework/JsonLite.C
    src/framework/StackSampler.C
    src/framework/RunTelemetry.C
)

if(RFAM_HAVE_COOLPROP)
    if(NOT EXISTS "${RFAM_COOLPROP_INCLUDE_DIR}/AbstractState.h" OR
       NOT EXISTS "${RFAM_COOLPROP_LIBRARY}")
        message(FATAL_ERROR "CoolProp requires RFAM_COOLPROP_INCLUDE_DIR and RFAM_COOLPROP_LIBRARY")
    endif()
    target_compile_definitions(rfam_core PRIVATE RFAM_HAVE_COOLPROP=1
        RFAM_COOLPROP_PROVENANCE="${RFAM_COOLPROP_PROVENANCE}")
    find_path(RFAM_COOLPROP_FMT_INCLUDE_DIR fmt/format.h
        HINTS "${RFAM_COOLPROP_INCLUDE_DIR}/../externals/fmtlib/include")
    if(NOT RFAM_COOLPROP_FMT_INCLUDE_DIR)
        message(FATAL_ERROR "CoolProp headers require fmt/format.h")
    endif()
    target_include_directories(rfam_core SYSTEM PRIVATE
        ${RFAM_COOLPROP_INCLUDE_DIR} ${RFAM_COOLPROP_FMT_INCLUDE_DIR})
    target_link_libraries(rfam_core PUBLIC ${RFAM_COOLPROP_LIBRARY} ${CMAKE_DL_LIBS})
endif()

if(RFAM_ENABLE_CHEMISTRY)
    find_path(RFAM_SUNDIALS_INCLUDE_DIR cvode/cvode.h
        HINTS ${SUNDIALS_ROOT} ENV SUNDIALS_ROOT PATH_SUFFIXES include)
    foreach(component cvode core nvecserial sunmatrixdense sunlinsoldense)
        string(TOUPPER ${component} component_upper)
        find_library(RFAM_SUNDIALS_${component_upper}_LIBRARY
            NAMES sundials_${component}
            HINTS ${SUNDIALS_ROOT} ENV SUNDIALS_ROOT PATH_SUFFIXES lib lib64)
        list(APPEND RFAM_SUNDIALS_LIBRARIES ${RFAM_SUNDIALS_${component_upper}_LIBRARY})
    endforeach()
    if(RFAM_SUNDIALS_INCLUDE_DIR AND NOT "${RFAM_SUNDIALS_LIBRARIES}" MATCHES "NOTFOUND")
        add_library(rfam_chemistry src/combustion/FiniteRate/TwoStepCO.C)
        target_include_directories(rfam_chemistry PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo
            ${CMAKE_CURRENT_SOURCE_DIR}/src/framework
            ${CMAKE_CURRENT_SOURCE_DIR}/src/combustion
            ${CMAKE_CURRENT_SOURCE_DIR}/src/combustion/FiniteRate
            ${CMAKE_CURRENT_SOURCE_DIR}/src/variables
            ${RFAM_SUNDIALS_INCLUDE_DIR})
        target_link_libraries(rfam_chemistry PUBLIC rfam_core ${RFAM_SUNDIALS_LIBRARIES})
        # A SUNDIALS built with MPI (Arch's system package, for one) drags
        # libmpi into the link line through its own shared libraries. Linking it
        # explicitly is required, not optional, or the final link fails with
        # "DSO missing from command line".
        find_package(MPI QUIET)
        if(MPI_CXX_FOUND)
            target_link_libraries(rfam_chemistry PUBLIC MPI::MPI_CXX)
            message(STATUS "RFAM chemistry: SUNDIALS requires MPI; linking MPI::MPI_CXX")
        endif()
        target_compile_definitions(rfam_chemistry PUBLIC RFAM_HAVE_CHEMISTRY=1)
        set(RFAM_HAVE_CHEMISTRY TRUE)
        message(STATUS "RFAM finite-rate chemistry enabled: ${RFAM_SUNDIALS_INCLUDE_DIR}")
    else()
        message(WARNING "SUNDIALS not found; finite-rate chemistry targets disabled")
    endif()
endif()

target_include_directories(rfam_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo
    ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/EquationOfState
    ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/MixingRules
    ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/RootSelection
    ${CMAKE_CURRENT_SOURCE_DIR}/src/variables
    ${CMAKE_CURRENT_SOURCE_DIR}/src/numerics/flux
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hotstart
    ${CMAKE_CURRENT_SOURCE_DIR}/src/numerics/preconditioning
    ${CMAKE_CURRENT_SOURCE_DIR}/src/numerics/reconstruction
    ${CMAKE_CURRENT_SOURCE_DIR}/src/framework
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/transport
    ${CMAKE_CURRENT_SOURCE_DIR}/src/turbulence
    ${CMAKE_CURRENT_SOURCE_DIR}/src/combustion
    ${CMAKE_CURRENT_SOURCE_DIR}/src/combustion/FPV
)
target_compile_definitions(rfam_core PUBLIC
    RFAM_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

if(RFAM_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(rfam_core PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(rfam_core PUBLIC RFAM_USE_OPENMP=1)
    else()
        message(WARNING "RFAM_ENABLE_OPENMP=ON but OpenMP was not found; building serial")
    endif()
endif()

# Run telemetry uses pthreads for the heartbeat/sampler threads and dladdr for
# stack symbolisation. Both are core Linux facilities; no external dependency.
find_package(Threads REQUIRED)
target_link_libraries(rfam_core PUBLIC Threads::Threads ${CMAKE_DL_LIBS})

if(MSVC)
    target_compile_options(rfam_core PRIVATE /W4)
else()
    target_compile_options(rfam_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_artifact_composite
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_artifact_composite.py)
    set_tests_properties(test_artifact_composite PROPERTIES
        LABELS "provenance;tooling"
        PASS_REGULAR_EXPRESSION "ARTIFACT_COMPOSITE_TEST_PASS"
        FAIL_REGULAR_EXPRESSION "FAILED;ARTIFACT_COMPOSITE.*FAIL"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME test_injector_flow_allocation_audit
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_injector_flow_allocation_audit.py)
    set_tests_properties(test_injector_flow_allocation_audit PROPERTIES
        LABELS "provenance;nonproduction"
        PASS_REGULAR_EXPRESSION "INJECTOR_FLOW_ALLOCATION_AUDIT_PASS"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME test_vortex_physics_contract
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_vortex_physics_contract.py)
    set_tests_properties(test_vortex_physics_contract PROPERTIES
        LABELS "provenance;nonproduction"
        PASS_REGULAR_EXPRESSION "VORTEX_PHYSICS_CONTRACT_TEST_PASS"
        FAIL_REGULAR_EXPRESSION "FAILED;VORTEX_PHYSICS_CONTRACT_FAIL"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
else()
    message(WARNING
        "Python3 not found; Python provenance/tooling gates are not registered")
endif()

function(add_rfam_test target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE rfam_core)
    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

add_rfam_test(test_cv_validation tests/test_cv_validation.cpp)
add_rfam_test(test_pr_mixing tests/test_pr_mixing.cpp)
add_rfam_test(test_pr_state_scan tests/test_pr_state_scan.cpp)
add_rfam_test(test_pr_thermo_identities tests/thermo/test_pr_thermo_identities.cpp)
add_rfam_test(test_pr_roots tests/thermo/test_pr_roots.cpp)
add_rfam_test(test_pr_low_pressure_closure tests/thermo/test_pr_low_pressure_closure.cpp)
add_rfam_test(test_cp_cv_positive tests/thermo/test_cp_cv_positive.cpp)
add_rfam_test(test_wale tests/thermo/test_wale.cpp)
add_rfam_test(test_cryogenic_transport_table
    tests/thermo/test_cryogenic_transport_table.cpp)
add_rfam_test(test_cryogenic_backend tests/thermo/test_cryogenic_backend.cpp)
add_rfam_test(test_droplet_evaporation tests/thermo/test_droplet_evaporation.cpp)
add_rfam_test(test_hybrid_helmholtz tests/thermo/test_hybrid_helmholtz.cpp)
add_rfam_test(test_hybrid_bubble tests/thermo/test_hybrid_bubble.cpp)
add_rfam_test(test_gerg_residual tests/thermo/test_gerg_residual.cpp)
add_rfam_test(test_realfluid_gerg tests/thermo/test_realfluid_gerg.cpp)
add_rfam_test(test_sandiego_nasa7_data
    tests/thermo/test_sandiego_nasa7_data.cpp)
set_tests_properties(test_sandiego_nasa7_data PROPERTIES
    LABELS "kinetics_data;nonproduction"
    PASS_REGULAR_EXPRESSION "SANDIEGO_NASA7_DATA_PASS"
    FAIL_REGULAR_EXPRESSION "SANDIEGO_NASA7_DATA_FAIL")
add_rfam_test(test_sandiego_kinetics_data
    tests/thermo/test_sandiego_kinetics_data.cpp)
set_tests_properties(test_sandiego_kinetics_data PROPERTIES
    LABELS "kinetics_data;nonproduction"
    PASS_REGULAR_EXPRESSION "SANDIEGO_KINETICS_DATA_PASS"
    FAIL_REGULAR_EXPRESSION "SANDIEGO_KINETICS_DATA_FAIL")
add_rfam_test(test_sandiego_kinetics_rate
    tests/thermo/test_sandiego_kinetics_rate.cpp)
set_tests_properties(test_sandiego_kinetics_rate PROPERTIES
    LABELS "kinetics_implementation;nonproduction"
    PASS_REGULAR_EXPRESSION "SANDIEGO_KINETICS_RATE_PASS"
    FAIL_REGULAR_EXPRESSION "SANDIEGO_KINETICS_RATE_FAIL")
if(RFAM_HAVE_CHEMISTRY)
    add_executable(test_sandiego_0d_ignition
        tests/thermo/test_sandiego_0d_ignition.cpp)
    target_link_libraries(test_sandiego_0d_ignition
        PRIVATE rfam_chemistry ${RFAM_SUNDIALS_LIBRARIES})
    target_include_directories(test_sandiego_0d_ignition PRIVATE
        ${RFAM_SUNDIALS_INCLUDE_DIR})
    add_test(NAME test_sandiego_0d_ignition COMMAND test_sandiego_0d_ignition)
    set_tests_properties(test_sandiego_0d_ignition PROPERTIES
        LABELS "kinetics_implementation;nonproduction"
        PASS_REGULAR_EXPRESSION "SANDIEGO_0D_IGNITION_PASS"
        FAIL_REGULAR_EXPRESSION "SANDIEGO_0D_IGNITION_FAIL")
    add_executable(test_sandiego_sector_bridge
        tests/thermo/test_sandiego_sector_bridge.cpp)
    target_link_libraries(test_sandiego_sector_bridge
        PRIVATE rfam_chemistry ${RFAM_SUNDIALS_LIBRARIES})
    target_include_directories(test_sandiego_sector_bridge PRIVATE
        ${RFAM_SUNDIALS_INCLUDE_DIR})
    add_test(NAME test_sandiego_sector_bridge COMMAND test_sandiego_sector_bridge)
    set_tests_properties(test_sandiego_sector_bridge PROPERTIES
        LABELS "kinetics_implementation;nonproduction"
        PASS_REGULAR_EXPRESSION "SANDIEGO_SECTOR_CHEM_BRIDGE_PASS"
        FAIL_REGULAR_EXPRESSION "SANDIEGO_SECTOR_CHEM_BRIDGE_FAIL")
endif()
if(RFAM_HAVE_TEQP)
    add_rfam_test(test_teqp_multifluid_mass_basis
        tests/thermo/test_teqp_multifluid_mass_basis.cpp)
    set_tests_properties(test_teqp_multifluid_mass_basis PROPERTIES
        LABELS "teqp;thermo;nonproduction"
        PASS_REGULAR_EXPRESSION "TEQP_MULTIFLUID_MASS_BASIS_PASS"
        FAIL_REGULAR_EXPRESSION "TEQP_MULTIFLUID_MASS_BASIS_FAIL")
endif()
add_rfam_test(test_hotstart_reconciler tests/thermo/test_hotstart_reconciler.cpp)
add_rfam_test(test_saturation_table tests/thermo/test_saturation_table.cpp)
add_rfam_test(test_hem_two_phase tests/thermo/test_hem_two_phase.cpp)
add_rfam_test(test_hem_hint tests/thermo/test_hem_hint.cpp)
add_rfam_test(test_hem_domain_hint tests/thermo/test_hem_domain_hint.cpp)
add_rfam_test(test_recovery_profile tests/thermo/test_recovery_profile.cpp)
add_rfam_test(test_hem_projection tests/thermo/test_hem_projection.cpp)
add_rfam_test(test_state_conversion_contract tests/thermo/test_state_conversion_contract.cpp)
add_rfam_test(test_diffusion_real_cv tests/thermo/test_diffusion_real_cv.cpp)
add_rfam_test(test_model_capability_gate tests/thermo/test_model_capability_gate.cpp)
add_rfam_test(test_pr_tp_flash tests/thermo/test_pr_tp_flash.cpp)
add_rfam_test(test_pr_phase_caloric tests/thermo/test_pr_phase_caloric.cpp)
add_rfam_test(test_pr_phase_thermo_owner_adapter
    tests/thermo/test_pr_phase_thermo_owner_adapter.cpp)
set_tests_properties(test_pr_phase_thermo_owner_adapter PROPERTIES
    LABELS "phase_chemistry;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_PHASE_THERMO_OWNER_ADAPTER_PASS"
    FAIL_REGULAR_EXPRESSION "PR_PHASE_THERMO_OWNER_ADAPTER_FAIL")
add_rfam_test(test_pr_phase_thermo_owner_displacement
    tests/thermo/test_pr_phase_thermo_owner_displacement.cpp)
set_tests_properties(test_pr_phase_thermo_owner_displacement PROPERTIES
    LABELS "phase_chemistry;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_PHASE_THERMO_OWNER_DISPLACEMENT_PASS"
    FAIL_REGULAR_EXPRESSION "PR_PHASE_THERMO_OWNER_DISPLACEMENT_FAIL")
add_rfam_test(test_pr_ph_flash tests/thermo/test_pr_ph_flash.cpp)
add_rfam_test(test_pr_rho_e_flash tests/thermo/test_pr_rho_e_flash.cpp)
add_rfam_test(test_pr_phase_displacement_authority
    tests/thermo/test_pr_phase_displacement_authority.cpp)
set_tests_properties(test_pr_phase_displacement_authority PROPERTIES
    LABELS "phase_chemistry;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_PHASE_DISPLACEMENT_AUTHORITY_PASS"
    FAIL_REGULAR_EXPRESSION "PR_PHASE_DISPLACEMENT_AUTHORITY_FAIL")
add_rfam_test(test_pr_residual_helmholtz_authority
    tests/thermo/test_pr_residual_helmholtz_authority.cpp)
target_include_directories(test_pr_residual_helmholtz_authority PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools)
set_tests_properties(test_pr_residual_helmholtz_authority PROPERTIES
    LABELS "thermo_authority;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_RESIDUAL_HELMHOLTZ_PASS"
    FAIL_REGULAR_EXPRESSION "PR_RESIDUAL_HELMHOLTZ_FAIL")
add_rfam_test(test_pr_closed_volume_heating
    tests/thermo/test_pr_closed_volume_heating.cpp)
set_tests_properties(test_pr_closed_volume_heating PROPERTIES
    LABELS "phase_box;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_CLOSED_VOLUME_HEATING_PASS"
    FAIL_REGULAR_EXPRESSION "PR_CLOSED_VOLUME_HEATING_FAIL")
add_rfam_test(test_pr_constant_pressure_heating
    tests/thermo/test_pr_constant_pressure_heating.cpp)
set_tests_properties(test_pr_constant_pressure_heating PROPERTIES
    LABELS "phase_box;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_CONSTANT_PRESSURE_HEATING_PASS"
    FAIL_REGULAR_EXPRESSION "PR_CONSTANT_PRESSURE_HEATING_FAIL")
add_rfam_test(test_pr_equilibrium_domain
    tests/thermo/test_pr_equilibrium_domain.cpp)
add_rfam_test(test_pr_equilibrium_short_step
    tests/thermo/test_pr_equilibrium_short_step.cpp)
set_tests_properties(test_pr_equilibrium_short_step PROPERTIES
    LABELS "phase_diagnostic;nonproduction")
add_rfam_test(test_acoustic_state_routing
    tests/thermo/test_acoustic_state_routing.cpp)
add_rfam_test(test_transport_state_routing
    tests/thermo/test_transport_state_routing.cpp)
add_rfam_test(test_transport_no_alloc
    tests/thermo/test_transport_no_alloc.cpp)
add_rfam_test(test_tabulated_ideal_caloric
    tests/thermo/test_tabulated_ideal_caloric.cpp)
add_rfam_test(test_two_phase_state_contract
    tests/thermo/test_two_phase_state_contract.cpp)
add_rfam_test(test_pr_two_phase_mixing_contract
    tests/thermo/test_pr_two_phase_mixing_contract.cpp)
set_tests_properties(test_pr_two_phase_mixing_contract PROPERTIES
    LABELS "phase_box;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_TWO_PHASE_MIXING_CONTRACT_PASS"
    FAIL_REGULAR_EXPRESSION "PR_TWO_PHASE_MIXING_CONTRACT_FAIL")
add_rfam_test(test_pr_vapor_chemistry_reflash_contract
    tests/thermo/test_pr_vapor_chemistry_reflash_contract.cpp)
set_tests_properties(test_pr_vapor_chemistry_reflash_contract PROPERTIES
    LABELS "phase_box;nonproduction"
    PASS_REGULAR_EXPRESSION "PR_VAPOR_CHEMISTRY_REFLASH_CONTRACT_PASS"
    FAIL_REGULAR_EXPRESSION "PR_VAPOR_CHEMISTRY_REFLASH_CONTRACT_FAIL")
add_rfam_test(test_preconditioner tests/test_preconditioner.cpp)
add_rfam_test(test_solver2d tests/test_solver2d.cpp)
add_rfam_test(test_implicit_stability tests/test_implicit_stability.cpp)
set_tests_properties(test_implicit_stability PROPERTIES
    LABELS "long;legacy_stability"
    TIMEOUT 300)
add_rfam_test(test_dual_time tests/test_dual_time.cpp)
set_tests_properties(test_dual_time PROPERTIES
    LABELS "framework;dual_time"
    TIMEOUT 300)
add_rfam_test(test_sod_tube tests/test_sod_tube.cpp)
add_rfam_test(test_domain3d_lifecycle tests/test_domain3d_lifecycle.cpp)
add_rfam_test(test_step_transaction tests/test_step_transaction.cpp)
add_rfam_test(test_restart_manifest tests/test_restart_manifest.cpp)
add_rfam_test(test_restart_manifest_codec tests/test_restart_manifest_codec.cpp)
add_rfam_test(test_restart_generation_store
    tests/test_restart_generation_store.cpp)
add_rfam_test(test_restart_atomic_publication_contract
    tests/test_restart_atomic_publication_contract.cpp)
add_rfam_test(test_restart_publication
    tests/test_restart_publication.cpp)
add_rfam_test(test_restart_generation_store_sector_clock_preflight
    tests/test_restart_generation_store_sector_clock_preflight.cpp)
add_rfam_test(test_sector_restart_owner_codec
    tests/test_sector_restart_owner_codec.cpp)
add_rfam_test(test_sector_restart_publication_stagers
    tests/test_sector_restart_publication_stagers.cpp)
set_tests_properties(
    test_restart_manifest
    test_restart_manifest_codec
    test_restart_generation_store
    test_restart_atomic_publication_contract
    test_restart_publication
    test_restart_generation_store_sector_clock_preflight
    test_sector_restart_owner_codec
    test_sector_restart_publication_stagers
    PROPERTIES LABELS "restart;nonproduction")
add_rfam_test(test_rhs_composition_preparation
    tests/test_rhs_composition_preparation.cpp)
set_tests_properties(test_rhs_composition_preparation PROPERTIES
    LABELS "kinetics_implementation;nonproduction"
    PASS_REGULAR_EXPRESSION "RHS_COMPOSITION_PREPARATION_PASS"
    FAIL_REGULAR_EXPRESSION "RHS_COMPOSITION_PREPARATION_FAIL")
add_rfam_test(test_condensed_carbon tests/test_condensed_carbon.cpp)
set_tests_properties(test_condensed_carbon PROPERTIES
    LABELS "condensed_carbon;nonproduction"
    PASS_REGULAR_EXPRESSION "CONDENSED_CARBON_V1_PASS"
    FAIL_REGULAR_EXPRESSION "CONDENSED_CARBON_V1_FAIL")
add_rfam_test(test_openmp_consistency tests/test_openmp_consistency.cpp)
add_rfam_test(test_isothermal_projection tests/test_isothermal_projection.cpp)
add_rfam_test(test_sector_metrics tests/test_sector_metrics.cpp)
add_rfam_test(test_sector_freestream tests/test_sector_freestream.cpp)
add_rfam_test(test_sector_turbulence tests/test_sector_turbulence.cpp)
add_rfam_test(test_sector_diffusion tests/test_sector_diffusion.cpp)
add_rfam_test(test_sector_diagnostics tests/test_sector_diagnostics.cpp)
add_rfam_test(test_sector_wall_boundary tests/test_sector_wall_boundary.cpp)
add_rfam_test(test_sector_injection tests/test_sector_injection.cpp)
add_rfam_test(test_sector_outlet tests/test_sector_outlet.cpp)
add_rfam_test(test_reacting_background tests/test_reacting_background.cpp)
add_rfam_test(test_sector_reacting_freestream tests/test_sector_reacting_freestream.cpp)
add_rfam_test(test_sector_performance tests/test_sector_performance.cpp)
add_rfam_test(test_sector_checkpoint tests/test_sector_checkpoint.cpp)
add_rfam_test(test_sector_checkpoint_remap tests/test_sector_checkpoint_remap.cpp)
add_rfam_test(test_sector_conservative_grid_transfer
    tests/test_sector_conservative_grid_transfer.cpp)
add_rfam_test(test_sector_coarse_correction
    tests/test_sector_coarse_correction.cpp)
add_rfam_test(test_steady_residual_evaluator
    tests/test_steady_residual_evaluator.cpp)
add_rfam_test(test_sector_fas_defect
    tests/test_sector_fas_defect.cpp)
add_rfam_test(test_sector_residual_field
    tests/test_sector_residual_field.cpp)
add_rfam_test(test_sector_galerkin_coarse_stencil
    tests/test_sector_galerkin_coarse_stencil.cpp)
add_rfam_test(test_sector_pressure_schur
    tests/test_sector_pressure_schur.cpp)
add_rfam_test(test_sector_structured_momentum_inverse
    tests/test_sector_structured_momentum_inverse.cpp)
add_rfam_test(test_sector_structured_pressure_schur
    tests/test_sector_structured_pressure_schur.cpp)
add_rfam_test(test_sector_output tests/test_sector_output.cpp)
add_rfam_test(test_sector_restart_equivalence_warm
    tests/test_sector_restart_equivalence_warm.cpp)
add_rfam_test(test_cht2d tests/test_cht2d.cpp)
target_sources(test_cht2d PRIVATE src/combustion/FPV/FPVTable.C)
add_rfam_test(test_residual_integrator tests/test_residual_integrator.cpp)
add_rfam_test(test_mach_residual_integrator tests/test_mach_residual_integrator.cpp)

option(RFAM_ENABLE_PRECICE "Build the preCICE 3 thermal coupling example" OFF)
if(RFAM_ENABLE_PRECICE)
    find_package(precice 3 CONFIG REQUIRED)
    add_executable(rfam_precice_thermal apps/rfam_precice_thermal.cpp)
    target_link_libraries(rfam_precice_thermal PRIVATE rfam_core precice::precice)
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_precice_residual_coupling
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_precice_residual_coupling.py
                $<TARGET_FILE:rfam_precice_thermal>
                ${CMAKE_CURRENT_SOURCE_DIR}/cases/precice_residual/precice-config.xml)
        set_tests_properties(test_precice_residual_coupling PROPERTIES TIMEOUT 90)
    endif()
endif()

add_executable(test_mvp1d_frame
    apps/solver_mvp_1d.cpp
    src/combustion/FPV/FPVTable.C)
target_link_libraries(test_mvp1d_frame PRIVATE rfam_core)
add_test(NAME test_mvp1d_frame COMMAND test_mvp1d_frame)
set_tests_properties(test_mvp1d_frame PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(RFAM_HAVE_CHEMISTRY)
    add_executable(test_phase7_operator_compile tests/test_phase7_operator_compile.cpp)
    target_link_libraries(test_phase7_operator_compile PRIVATE rfam_chemistry)
    add_test(NAME test_phase7_operator_compile COMMAND test_phase7_operator_compile)

    add_executable(test_phase7_1d tests/test_phase7_1d.cpp)
    target_link_libraries(test_phase7_1d PRIVATE rfam_chemistry)
    add_test(NAME test_phase7_1d COMMAND test_phase7_1d)

    add_executable(test_sector_chemistry tests/test_sector_chemistry.cpp)
    target_link_libraries(test_sector_chemistry PRIVATE rfam_chemistry)
    add_test(NAME test_sector_chemistry COMMAND test_sector_chemistry)

    add_executable(test_real_fluid_chemistry_energy
        tests/thermo/test_real_fluid_chemistry_energy.cpp)
    target_link_libraries(test_real_fluid_chemistry_energy PRIVATE rfam_chemistry)
    add_test(NAME test_real_fluid_chemistry_energy
        COMMAND test_real_fluid_chemistry_energy)

    add_executable(test_real_fluid_phase_owner_cvode
        tests/thermo/test_real_fluid_phase_owner_cvode.cpp)
    target_link_libraries(test_real_fluid_phase_owner_cvode PRIVATE rfam_chemistry)
    add_test(NAME test_real_fluid_phase_owner_cvode
        COMMAND test_real_fluid_phase_owner_cvode)
    set_tests_properties(test_real_fluid_phase_owner_cvode PROPERTIES
        LABELS "phase_chemistry;nonproduction"
        PASS_REGULAR_EXPRESSION "REAL_FLUID_PHASE_OWNER_CVODE_PASS"
        FAIL_REGULAR_EXPRESSION "REAL_FLUID_PHASE_OWNER_CVODE_FAIL")

    add_executable(test_2d_reacting tests/test_2d_reacting.cpp)
    target_link_libraries(test_2d_reacting PRIVATE rfam_chemistry)
    add_test(NAME test_2d_reacting COMMAND test_2d_reacting)

    add_executable(test_pinfin_cylinder tests/test_pinfin_cylinder.cpp)
    target_link_libraries(test_pinfin_cylinder PRIVATE rfam_core)
    add_test(NAME test_pinfin_cylinder COMMAND test_pinfin_cylinder)
    set_tests_properties(test_pinfin_cylinder PROPERTIES
        TIMEOUT 900
        LABELS "physical_validation;nonproduction")

    add_executable(test_1d_flame tests/test_1d_flame.cpp)
    target_link_libraries(test_1d_flame PRIVATE rfam_chemistry)
    add_test(NAME test_1d_flame COMMAND test_1d_flame)
    set_tests_properties(test_1d_flame PROPERTIES
        TIMEOUT 900
        LABELS "physical_validation;nonproduction")

    add_executable(test_stage15_species tests/test_stage15_species.cpp)
    target_link_libraries(test_stage15_species PRIVATE rfam_chemistry)
    add_test(NAME test_stage15_species COMMAND test_stage15_species)

    add_executable(test_mechanism_reader tests/test_mechanism_reader.cpp)
    target_link_libraries(test_mechanism_reader PRIVATE rfam_chemistry)
    add_test(NAME test_mechanism_reader COMMAND test_mechanism_reader)
    set_tests_properties(test_mechanism_reader PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set(RFAM_CANTERA_GENERATOR
        ${CMAKE_CURRENT_SOURCE_DIR}/scripts/cantera_reference.py)
    set(RFAM_CANTERA_REFERENCE
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/cantera_gri30_propane_oxygen_v1.txt)
    set(RFAM_CANTERA_MECHANISM
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/cantera_gri30_v3_0.yaml)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        ${RFAM_CANTERA_GENERATOR}
        ${RFAM_CANTERA_REFERENCE}
        ${RFAM_CANTERA_MECHANISM})
    file(SHA256 ${RFAM_CANTERA_GENERATOR}
        RFAM_CANTERA_GENERATOR_SHA256)
    file(SHA256 ${RFAM_CANTERA_REFERENCE}
        RFAM_CANTERA_REFERENCE_SHA256)
    file(SHA256 ${RFAM_CANTERA_MECHANISM}
        RFAM_CANTERA_MECHANISM_SHA256)

    add_executable(test_cantera_kinetics_gate tests/test_cantera_compare.cpp)
    target_link_libraries(test_cantera_kinetics_gate PRIVATE rfam_chemistry)
    target_compile_definitions(test_cantera_kinetics_gate PRIVATE
        RFAM_CANTERA_GENERATOR_SHA256="${RFAM_CANTERA_GENERATOR_SHA256}"
        RFAM_CANTERA_REFERENCE_SHA256="${RFAM_CANTERA_REFERENCE_SHA256}"
        RFAM_CANTERA_MECHANISM_SHA256="${RFAM_CANTERA_MECHANISM_SHA256}")
    add_test(NAME test_cantera_kinetics_gate
        COMMAND test_cantera_kinetics_gate
            ${RFAM_CANTERA_GENERATOR}
            ${RFAM_CANTERA_REFERENCE}
            ${RFAM_CANTERA_MECHANISM})
    set_tests_properties(test_cantera_kinetics_gate PROPERTIES
        LABELS "kinetics_diagnostic;nonproduction"
        PASS_REGULAR_EXPRESSION "CANTERA_KINETICS_GATE_BLOCKED"
        FAIL_REGULAR_EXPRESSION
            "CANTERA_COMPARE_PASS;CANTERA_COMPARE_FAIL;CANTERA_KINETICS_GATE_FAIL"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME test_cantera_kinetics_provenance_gate
        COMMAND ${CMAKE_COMMAND}
            -DRFAM_GATE_EXE=$<TARGET_FILE:test_cantera_kinetics_gate>
            -DRFAM_GENERATOR=${RFAM_CANTERA_GENERATOR}
            -DRFAM_REFERENCE=${RFAM_CANTERA_REFERENCE}
            -DRFAM_MECHANISM=${RFAM_CANTERA_MECHANISM}
            -DRFAM_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/cantera_kinetics_provenance_gate
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cantera_kinetics_provenance_gate.cmake)
    set_tests_properties(test_cantera_kinetics_provenance_gate PROPERTIES
        LABELS "kinetics_diagnostic;nonproduction"
        PASS_REGULAR_EXPRESSION "CANTERA_KINETICS_PROVENANCE_GATE_PASS")

    set(RFAM_R2_GENERATOR
        ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_r2_instantaneous_reference.py)
    set(RFAM_R2_REFERENCE
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/r2_instantaneous_reference_v1.json)
    set(RFAM_R2_CANTERA_MECHANISM
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/rfam_twostep_cantera_v1.yaml)
    set(RFAM_R2_THERMO
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/cantera_gri30_v3_0.yaml)
    set(RFAM_R2_RFAM_MECHANISM
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/twostep_co.yaml)
    set(RFAM_R2_SPECIES_REGISTRY
        ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/SpeciesRegistry.H)
    set(RFAM_R2_MASS_ACTION_THERMO
        ${CMAKE_CURRENT_SOURCE_DIR}/src/combustion/FiniteRate/MassActionThermo.H)
    set(RFAM_R2_NASA7_THERMO
        ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/NASA7Thermo.H)
    set(RFAM_R2_MECHANISM_READER
        ${CMAKE_CURRENT_SOURCE_DIR}/src/combustion/FiniteRate/MechanismReader.H)
    set(RFAM_R2_THERMO_PROVENANCE_HEADER
        ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/ThermoDataProvenance.H)
    set(RFAM_R2_THERMO_PROVENANCE_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo/ThermoDataProvenance.C)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        ${RFAM_R2_GENERATOR}
        ${RFAM_R2_REFERENCE}
        ${RFAM_R2_CANTERA_MECHANISM}
        ${RFAM_R2_THERMO}
        ${RFAM_R2_RFAM_MECHANISM}
        ${RFAM_R2_SPECIES_REGISTRY}
        ${RFAM_R2_MASS_ACTION_THERMO}
        ${RFAM_R2_NASA7_THERMO}
        ${RFAM_R2_MECHANISM_READER}
        ${RFAM_R2_THERMO_PROVENANCE_HEADER}
        ${RFAM_R2_THERMO_PROVENANCE_SOURCE})
    file(SHA256 ${RFAM_R2_GENERATOR} RFAM_R2_GENERATOR_SHA256)
    file(SHA256 ${RFAM_R2_REFERENCE} RFAM_R2_REFERENCE_SHA256)
    file(SHA256 ${RFAM_R2_CANTERA_MECHANISM}
        RFAM_R2_CANTERA_MECHANISM_SHA256)
    file(SHA256 ${RFAM_R2_THERMO} RFAM_R2_THERMO_SHA256)
    file(SHA256 ${RFAM_R2_RFAM_MECHANISM}
        RFAM_R2_RFAM_MECHANISM_SHA256)
    file(SHA256 ${RFAM_R2_SPECIES_REGISTRY}
        RFAM_R2_SPECIES_REGISTRY_SHA256)
    file(SHA256 ${RFAM_R2_MASS_ACTION_THERMO}
        RFAM_R2_MASS_ACTION_THERMO_SHA256)
    file(SHA256 ${RFAM_R2_NASA7_THERMO}
        RFAM_R2_NASA7_THERMO_SHA256)
    file(SHA256 ${RFAM_R2_MECHANISM_READER}
        RFAM_R2_MECHANISM_READER_SHA256)
    file(SHA256 ${RFAM_R2_THERMO_PROVENANCE_HEADER}
        RFAM_R2_THERMO_PROVENANCE_HEADER_SHA256)
    file(SHA256 ${RFAM_R2_THERMO_PROVENANCE_SOURCE}
        RFAM_R2_THERMO_PROVENANCE_SOURCE_SHA256)

    add_executable(test_r2_instantaneous_compare
        tests/test_r2_instantaneous_compare.cpp)
    target_link_libraries(test_r2_instantaneous_compare PRIVATE rfam_chemistry)
    target_compile_definitions(test_r2_instantaneous_compare PRIVATE
        RFAM_R2_GENERATOR_SHA256="${RFAM_R2_GENERATOR_SHA256}"
        RFAM_R2_REFERENCE_SHA256="${RFAM_R2_REFERENCE_SHA256}"
        RFAM_R2_CANTERA_MECHANISM_SHA256="${RFAM_R2_CANTERA_MECHANISM_SHA256}"
        RFAM_R2_THERMO_SHA256="${RFAM_R2_THERMO_SHA256}"
        RFAM_R2_RFAM_MECHANISM_SHA256="${RFAM_R2_RFAM_MECHANISM_SHA256}"
        RFAM_R2_SPECIES_REGISTRY_SHA256="${RFAM_R2_SPECIES_REGISTRY_SHA256}"
        RFAM_R2_MASS_ACTION_THERMO_SHA256="${RFAM_R2_MASS_ACTION_THERMO_SHA256}"
        RFAM_R2_NASA7_THERMO_SHA256="${RFAM_R2_NASA7_THERMO_SHA256}"
        RFAM_R2_MECHANISM_READER_SHA256="${RFAM_R2_MECHANISM_READER_SHA256}"
        RFAM_R2_THERMO_PROVENANCE_HEADER_SHA256="${RFAM_R2_THERMO_PROVENANCE_HEADER_SHA256}"
        RFAM_R2_THERMO_PROVENANCE_SOURCE_SHA256="${RFAM_R2_THERMO_PROVENANCE_SOURCE_SHA256}")

    add_executable(rfam_r2_0d_runner
        tools/r2_0d_rfam_runner.cpp)
    target_link_libraries(rfam_r2_0d_runner PRIVATE rfam_chemistry)

    set(RFAM_R2_LOW_CONCENTRATION_PROBE
        ${CMAKE_CURRENT_BINARY_DIR}/rfam_r2_low_concentration_regularization_probe)
    add_custom_command(
        OUTPUT ${RFAM_R2_LOW_CONCENTRATION_PROBE}
        COMMAND /usr/bin/g++
            -std=c++17 -O2 -Wall -Wextra -pedantic
            -Isrc -Isrc/thermo -Isrc/framework
            -Isrc/combustion -Isrc/combustion/FiniteRate -Isrc/variables
            -I${RFAM_SUNDIALS_INCLUDE_DIR}
            -Isrc/thermo/EquationOfState -Isrc/thermo/MixingRules
            -Isrc/thermo/RootSelection -Isrc/numerics/flux
            -Isrc/numerics/preconditioning -Isrc/numerics/reconstruction
            -Isrc/mesh -Isrc/transport -Isrc/turbulence
            -Isrc/combustion/FPV
            tools/r2_low_concentration_regularization_probe.cpp
            -o ${RFAM_R2_LOW_CONCENTRATION_PROBE}
        DEPENDS
            tools/r2_low_concentration_regularization_probe.cpp
            src/combustion/FiniteRate/ArrheniusRate.H
            src/combustion/FiniteRate/MassActionThermo.H
            src/combustion/FiniteRate/RhsCompositionPreparation.H
            src/combustion/SpeciesMapper.H
            src/thermo/NASA7Thermo.H
            src/thermo/SpeciesRegistry.H
            src/thermo/ThermoDataProvenance.H
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        VERBATIM)
    add_custom_target(rfam_r2_low_concentration_regularization_probe_artifact ALL
        DEPENDS ${RFAM_R2_LOW_CONCENTRATION_PROBE})

    add_executable(test_r2_thermo_envelope
        tests/test_r2_thermo_envelope.cpp)
    target_link_libraries(test_r2_thermo_envelope PRIVATE rfam_chemistry)
    add_test(NAME test_r2_thermo_envelope
        COMMAND test_r2_thermo_envelope)
    set_tests_properties(test_r2_thermo_envelope PROPERTIES
        LABELS "kinetics_implementation;nonproduction"
        PASS_REGULAR_EXPRESSION "R2_THERMO_ENVELOPE_TEST_PASS"
        FAIL_REGULAR_EXPRESSION "R2_THERMO_ENVELOPE_TEST_FAIL")

    if(Python3_Interpreter_FOUND)
        add_test(NAME test_r2_low_concentration_regularization
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_r2_low_concentration_regularization.py
                --root ${CMAKE_CURRENT_SOURCE_DIR}
                --build-dir ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(test_r2_low_concentration_regularization PROPERTIES
            LABELS "kinetics_implementation;nonproduction"
            PASS_REGULAR_EXPRESSION
                "R2_LOW_CONCENTRATION_REGULARIZATION_TEST_PASS"
            FAIL_REGULAR_EXPRESSION
                "R2_LOW_CONCENTRATION_REGULARIZATION_TEST_FAIL"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        add_test(NAME test_r2_same_mechanism_trajectory_contract
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_r2_same_mechanism_trajectory_contract.py)
        set_tests_properties(test_r2_same_mechanism_trajectory_contract PROPERTIES
            LABELS "kinetics_implementation;nonproduction"
            PASS_REGULAR_EXPRESSION
                "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_TEST_PASS"
            FAIL_REGULAR_EXPRESSION
                "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_TEST_FAIL"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        add_test(NAME test_r2_same_mechanism_trajectory
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_r2_same_mechanism_trajectory.py
                --root ${CMAKE_CURRENT_SOURCE_DIR}
                --runner $<TARGET_FILE:rfam_r2_0d_runner>)
        set_tests_properties(test_r2_same_mechanism_trajectory PROPERTIES
            LABELS "kinetics_implementation;nonproduction"
            PASS_REGULAR_EXPRESSION
                "R2_SAME_MECHANISM_TRAJECTORY_TEST_PASS"
            FAIL_REGULAR_EXPRESSION
                "R2_SAME_MECHANISM_TRAJECTORY_TEST_FAIL"
            TIMEOUT 300
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        add_test(NAME test_r2_0d_rfam_protocol
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_r2_0d_rfam_protocol.py
                --root ${CMAKE_CURRENT_SOURCE_DIR}
                --runner $<TARGET_FILE:rfam_r2_0d_runner>)
        set_tests_properties(test_r2_0d_rfam_protocol PROPERTIES
            LABELS "kinetics_implementation;nonproduction"
            PASS_REGULAR_EXPRESSION "R2_0D_RFAM_PROTOCOL_TEST_PASS"
            FAIL_REGULAR_EXPRESSION "R2_0D_RFAM_PROTOCOL_TEST_FAIL"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        add_test(NAME test_r2_instantaneous_equivalence
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_r2_instantaneous_equivalence.py
                --probe $<TARGET_FILE:test_r2_instantaneous_compare>
                --generator ${RFAM_R2_GENERATOR}
                --reference ${RFAM_R2_REFERENCE}
                --mechanism ${RFAM_R2_CANTERA_MECHANISM}
                --thermo ${RFAM_R2_THERMO}
                --rfam-mechanism ${RFAM_R2_RFAM_MECHANISM}
                --species-registry ${RFAM_R2_SPECIES_REGISTRY}
                --mass-action-thermo ${RFAM_R2_MASS_ACTION_THERMO}
                --nasa7-thermo ${RFAM_R2_NASA7_THERMO}
                --mechanism-reader ${RFAM_R2_MECHANISM_READER}
                --thermo-provenance-header ${RFAM_R2_THERMO_PROVENANCE_HEADER}
                --thermo-provenance-source ${RFAM_R2_THERMO_PROVENANCE_SOURCE})
        set_tests_properties(test_r2_instantaneous_equivalence PROPERTIES
            LABELS "kinetics_implementation;nonproduction"
            PASS_REGULAR_EXPRESSION "R2_INSTANTANEOUS_EQUIVALENCE_PASS"
            FAIL_REGULAR_EXPRESSION "R2_INSTANTANEOUS_EQUIVALENCE_FAIL"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        add_test(NAME test_r2_instantaneous_provenance_gate
            COMMAND ${CMAKE_COMMAND}
                -DRFAM_R2_PROBE=$<TARGET_FILE:test_r2_instantaneous_compare>
                -DRFAM_R2_GENERATOR=${RFAM_R2_GENERATOR}
                -DRFAM_R2_REFERENCE=${RFAM_R2_REFERENCE}
                -DRFAM_R2_CANTERA_MECHANISM=${RFAM_R2_CANTERA_MECHANISM}
                -DRFAM_R2_THERMO=${RFAM_R2_THERMO}
                -DRFAM_R2_RFAM_MECHANISM=${RFAM_R2_RFAM_MECHANISM}
                -DRFAM_R2_SPECIES_REGISTRY=${RFAM_R2_SPECIES_REGISTRY}
                -DRFAM_R2_MASS_ACTION_THERMO=${RFAM_R2_MASS_ACTION_THERMO}
                -DRFAM_R2_NASA7_THERMO=${RFAM_R2_NASA7_THERMO}
                -DRFAM_R2_MECHANISM_READER=${RFAM_R2_MECHANISM_READER}
                -DRFAM_R2_THERMO_PROVENANCE_HEADER=${RFAM_R2_THERMO_PROVENANCE_HEADER}
                -DRFAM_R2_THERMO_PROVENANCE_SOURCE=${RFAM_R2_THERMO_PROVENANCE_SOURCE}
                -DRFAM_R2_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/r2_instantaneous_provenance_gate
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/r2_instantaneous_provenance_gate.cmake)
        set_tests_properties(test_r2_instantaneous_provenance_gate PROPERTIES
            LABELS "kinetics_implementation;nonproduction"
            PASS_REGULAR_EXPRESSION "R2_INSTANTANEOUS_PROVENANCE_GATE_PASS")

        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c
                "import cantera as ct; raise SystemExit(0 if ct.__version__ == '3.1.0' else 1)"
            RESULT_VARIABLE RFAM_R2_CANTERA_AVAILABLE
            OUTPUT_QUIET
            ERROR_QUIET)
        if(RFAM_R2_CANTERA_AVAILABLE EQUAL 0)
            add_test(NAME test_r2_six_species_equilibrium_probe
                COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_r2_six_species_equilibrium_probe.py
                    --root ${CMAKE_CURRENT_SOURCE_DIR})
            set_tests_properties(test_r2_six_species_equilibrium_probe PROPERTIES
                LABELS "kinetics_physics;nonproduction"
                PASS_REGULAR_EXPRESSION "R2_SIX_SPECIES_EQUILIBRIUM_TEST_PASS"
                FAIL_REGULAR_EXPRESSION "R2_SIX_SPECIES_EQUILIBRIUM_TEST_FAIL"
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

            add_test(NAME test_r2_reference_regeneration_gate
                COMMAND ${CMAKE_COMMAND}
                    -DRFAM_PYTHON=${Python3_EXECUTABLE}
                    -DRFAM_R2_GENERATOR=${RFAM_R2_GENERATOR}
                    -DRFAM_R2_REFERENCE=${RFAM_R2_REFERENCE}
                    -DRFAM_R2_CANTERA_MECHANISM=${RFAM_R2_CANTERA_MECHANISM}
                    -DRFAM_R2_THERMO=${RFAM_R2_THERMO}
                    -DRFAM_R2_SPECIES_REGISTRY=${RFAM_R2_SPECIES_REGISTRY}
                    -DRFAM_R2_REGEN_DIR=${CMAKE_CURRENT_BINARY_DIR}/r2_reference_regeneration_gate
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/r2_reference_regeneration_gate.cmake)
            set_tests_properties(test_r2_reference_regeneration_gate PROPERTIES
                LABELS "kinetics_implementation;provenance;nonproduction"
                PASS_REGULAR_EXPRESSION "R2_REFERENCE_REGENERATION_GATE_PASS")
        else()
            message(WARNING
                "Pinned Cantera 3.1.0 unavailable; R2 byte-for-byte regeneration gate is not registered")
        endif()
    else()
        message(WARNING
            "Python3 not found; R2 instantaneous equivalence gate is not registered")
    endif()
endif()

add_executable(solver_1d_validation apps/solver_1d_validation.cpp)
target_link_libraries(solver_1d_validation PRIVATE rfam_core)
add_test(NAME test_1d_validation COMMAND solver_1d_validation)
set_tests_properties(test_1d_validation PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(solver_nozzle_1d apps/solver_nozzle_1d.cpp)
target_link_libraries(solver_nozzle_1d PRIVATE rfam_core)
add_test(NAME test_nozzle_allmach COMMAND solver_nozzle_1d)
set_tests_properties(test_nozzle_allmach PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(solver_mayer_o2 apps/solver_mayer_o2.cpp)
target_compile_definitions(solver_mayer_o2 PRIVATE
    RFAM_GIT_COMMIT=""
    RFAM_ALLOW_DIAGNOSTIC_HEM_ACOUSTICS)
target_link_libraries(solver_mayer_o2 PRIVATE rfam_core)

if(RFAM_BUILD_DIAGNOSTICS)
    add_executable(dbg119 _archive/diag/dbg119.cpp)
    target_link_libraries(dbg119 PRIVATE rfam_core)
    add_executable(dbg_qclass _archive/diag/dbg_qclass.cpp)
    target_compile_definitions(dbg_qclass PRIVATE RFAM_THERMO_DIAGNOSTICS)
    target_link_libraries(dbg_qclass PRIVATE rfam_core)
    add_executable(dbg_pr_fail _archive/diag/dbg_pr_fail.cpp)
    target_link_libraries(dbg_pr_fail PRIVATE rfam_core)
    if(RFAM_HAVE_TEQP)
        target_compile_definitions(dbg_pr_fail PRIVATE RFAM_HAVE_TEQP)
    endif()
endif()

add_executable(extract_hotstart_field tools/extract_hotstart_field.cpp)
target_link_libraries(extract_hotstart_field PRIVATE rfam_core)
add_executable(export_mesh_stl tools/export_mesh_stl.cpp)
target_link_libraries(export_mesh_stl PRIVATE rfam_core)
# Run telemetry contract: JSON round-trip, manifest sources, heartbeat under
# an idle solver, ledger accounting, and summary on a clean exit.
add_executable(test_run_telemetry tests/test_run_telemetry.cpp)
target_link_libraries(test_run_telemetry PRIVATE rfam_core)
add_test(NAME test_run_telemetry COMMAND test_run_telemetry)
set_tests_properties(test_run_telemetry PROPERTIES LABELS "telemetry;contract")

add_executable(test_tetra_pressure_outlet tests/test_tetra_pressure_outlet.cpp)
target_link_libraries(test_tetra_pressure_outlet PRIVATE rfam_core)
add_test(NAME test_tetra_pressure_outlet COMMAND test_tetra_pressure_outlet)
set_tests_properties(test_tetra_pressure_outlet PROPERTIES SKIP_RETURN_CODE 77)

# Single-owner face invariant guard. Deliberately does NOT skip: this test must
# run in every configuration, because the drift it guards against (a duplicated
# |u_n|+a formula with its own sound-speed floor) previously went unnoticed in
# exactly the builds where the tetra path was exercised.
add_executable(test_tetra_face_invariant tests/test_tetra_face_invariant.cpp)
target_link_libraries(test_tetra_face_invariant PRIVATE rfam_core)
add_test(NAME test_tetra_face_invariant COMMAND test_tetra_face_invariant)

# AUSM+-up splitting measurement: uniform-flow identity (P1) and acoustic-mode
# damping (P2), read directly off the production flux. No mesh, no EOS, no time
# integration, so it runs in every configuration. This is the gate a corrected
# M^4 form must pass; see the header of the test for the measured 2026-10-09
# result and why the M^4 branch is currently rejected.
add_executable(test_acoustic_decay_1d tests/test_acoustic_decay_1d.cpp)
target_link_libraries(test_acoustic_decay_1d PRIVATE rfam_core)
add_test(NAME test_acoustic_decay_1d COMMAND test_acoustic_decay_1d)

add_executable(test_tetra_inlet_allocation tests/test_tetra_inlet_allocation.cpp)
target_link_libraries(test_tetra_inlet_allocation PRIVATE rfam_core)
add_test(NAME test_tetra_inlet_allocation COMMAND test_tetra_inlet_allocation)

add_executable(test_tetra_blueprint_wall tests/test_tetra_blueprint_wall.cpp)
target_link_libraries(test_tetra_blueprint_wall PRIVATE rfam_core)
add_test(NAME test_tetra_blueprint_wall COMMAND test_tetra_blueprint_wall
    ${CMAKE_CURRENT_SOURCE_DIR}/work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet)
set_tests_properties(test_tetra_blueprint_wall PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_tetra_moving_slip_wall tests/test_tetra_moving_slip_wall.cpp)
target_link_libraries(test_tetra_moving_slip_wall PRIVATE rfam_core)
add_test(NAME test_tetra_moving_slip_wall COMMAND test_tetra_moving_slip_wall)

add_executable(test_tetra_inlet_psy tests/thermo/test_tetra_inlet_psy.cpp)
target_link_libraries(test_tetra_inlet_psy PRIVATE rfam_core)
add_test(NAME test_tetra_inlet_psy COMMAND test_tetra_inlet_psy)
set_tests_properties(test_tetra_inlet_psy PROPERTIES SKIP_RETURN_CODE 77)

add_executable(test_tetra_conservation tests/test_tetra_conservation.cpp)
target_link_libraries(test_tetra_conservation PRIVATE rfam_core)
add_test(NAME test_tetra_conservation COMMAND test_tetra_conservation)
add_executable(test_tetra_physical_retry tests/test_tetra_physical_retry.cpp)
target_link_libraries(test_tetra_physical_retry PRIVATE rfam_core)
add_test(NAME test_tetra_physical_retry COMMAND test_tetra_physical_retry)
add_executable(test_tetra_physical_checkpoint tests/test_tetra_physical_checkpoint.cpp)
target_link_libraries(test_tetra_physical_checkpoint PRIVATE rfam_core)
add_test(NAME test_tetra_physical_checkpoint COMMAND test_tetra_physical_checkpoint)

if(RFAM_HAVE_COOLPROP)
    add_executable(probe_transport tools/probe_transport.cpp)
    target_link_libraries(probe_transport PRIVATE ${RFAM_COOLPROP_LIBRARY} ${CMAKE_DL_LIBS})
    find_path(RFAM_COOLPROP_FMT_INCLUDE_DIR fmt/format.h
        HINTS "${RFAM_COOLPROP_INCLUDE_DIR}/../externals/fmtlib/include")
    target_include_directories(probe_transport PRIVATE
        ${RFAM_COOLPROP_INCLUDE_DIR} ${RFAM_COOLPROP_FMT_INCLUDE_DIR})
endif()

# B (benchmark/diagnostic) mesh: ~950 KB, small enough to track in-tree. It lives
# beside the source rather than under work/ so that a fresh clone can run the
# viscous, gradient and molecular-spectral tests. Declaring it once here means a
# future move is a one-line change instead of a 12-site search.
set(RFAM_TETRA_MESH_B "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/fluid_volume_B.rfamtet"
    CACHE PATH "Benchmark tetrahedral mesh used by the transport/gradient tests")

if(NOT EXISTS "${RFAM_TETRA_MESH_B}")
    message(STATUS
        "Tetra B-mesh tests SKIPPED: benchmark mesh not present at ${RFAM_TETRA_MESH_B}. "
        "Restore it (see docs/PROJECT_TREE.md section 10.8) before reporting coverage.")
endif()

add_executable(test_tetra_viscous_flux tests/test_tetra_viscous_flux.cpp)
target_link_libraries(test_tetra_viscous_flux PRIVATE rfam_core)
add_test(NAME test_tetra_viscous_flux COMMAND test_tetra_viscous_flux ${RFAM_TETRA_MESH_B})
set_tests_properties(test_tetra_viscous_flux PROPERTIES SKIP_RETURN_CODE 77 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_tetra_viscous_wall tests/test_tetra_viscous_wall.cpp)
target_link_libraries(test_tetra_viscous_wall PRIVATE rfam_core)
add_test(NAME test_tetra_viscous_wall COMMAND test_tetra_viscous_wall
    ${RFAM_TETRA_MESH_B})
set_tests_properties(test_tetra_viscous_wall PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_tetra_turbulence tests/test_tetra_turbulence.cpp)
target_link_libraries(test_tetra_turbulence PRIVATE rfam_core)
add_test(NAME test_tetra_turbulence COMMAND test_tetra_turbulence ${RFAM_TETRA_MESH_B})
set_tests_properties(test_tetra_turbulence PROPERTIES SKIP_RETURN_CODE 77
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_tetra_viscous_fv tests/test_tetra_viscous_fv.cpp)
target_link_libraries(test_tetra_viscous_fv PRIVATE rfam_core)
add_test(NAME test_tetra_viscous_fv COMMAND test_tetra_viscous_fv
    ${RFAM_TETRA_MESH_B})
set_tests_properties(test_tetra_viscous_fv PROPERTIES SKIP_RETURN_CODE 77
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_tetra_molecular_admission tests/test_tetra_molecular_admission.cpp)
target_link_libraries(test_tetra_molecular_admission PRIVATE rfam_core)
add_test(NAME test_tetra_molecular_admission COMMAND test_tetra_molecular_admission)

add_executable(test_tetra_molecular_spectral tests/test_tetra_molecular_spectral.cpp)
target_link_libraries(test_tetra_molecular_spectral PRIVATE rfam_core)
add_test(NAME test_tetra_molecular_spectral COMMAND test_tetra_molecular_spectral
    ${RFAM_TETRA_MESH_B})
if(RFAM_HAVE_COOLPROP)
    set(RFAM_TETRA_PHASE_TEST_TABLE ${CMAKE_CURRENT_SOURCE_DIR}/work/a406_rootcause/property_table_b.bin)
    if(NOT EXISTS "${RFAM_TETRA_PHASE_TEST_TABLE}")
        set(RFAM_TETRA_PHASE_TEST_TABLE ${CMAKE_CURRENT_SOURCE_DIR}/runs/property_table_b.bin)
    endif()
    if(EXISTS "${RFAM_TETRA_PHASE_TEST_TABLE}")
        add_test(NAME test_tetra_molecular_phase_gate COMMAND ${CMAKE_COMMAND}
            -DRFAM_TETRA_EXE=$<TARGET_FILE:rfam_v20_tetra>
            -DRFAM_TETRA_MESH=${RFAM_TETRA_MESH_B}
            -DRFAM_TETRA_TABLE=${RFAM_TETRA_PHASE_TEST_TABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/tetra_molecular_phase_gate.cmake)
        add_test(NAME test_tetra_periodic_phase_gate COMMAND ${CMAKE_COMMAND}
            -DRFAM_TETRA_EXE=$<TARGET_FILE:rfam_v20_tetra>
            -DRFAM_TETRA_MESH=${RFAM_TETRA_MESH_B}
            -DRFAM_TETRA_TABLE=${RFAM_TETRA_PHASE_TEST_TABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/tetra_periodic_phase_gate.cmake)
        add_test(NAME test_tetra_physical_retry_gate COMMAND ${CMAKE_COMMAND}
            -DRFAM_TETRA_EXE=$<TARGET_FILE:rfam_v20_tetra>
            -DRFAM_TETRA_MESH=${RFAM_TETRA_MESH_B}
            -DRFAM_TETRA_TABLE=${RFAM_TETRA_PHASE_TEST_TABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/tetra_physical_retry_gate.cmake)
        add_test(NAME test_tetra_physical_checkpoint_gate COMMAND ${CMAKE_COMMAND}
            -DRFAM_TETRA_EXE=$<TARGET_FILE:rfam_v20_tetra>
            -DRFAM_TETRA_MESH=${RFAM_TETRA_MESH_B}
            -DRFAM_TETRA_TABLE=${RFAM_TETRA_PHASE_TEST_TABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/tetra_physical_checkpoint_gate.cmake)
    endif()
endif()


add_executable(test_tetra_gradients tests/test_tetra_gradients.cpp)
target_link_libraries(test_tetra_gradients PRIVATE rfam_core)
add_test(NAME test_tetra_gradients COMMAND test_tetra_gradients ${RFAM_TETRA_MESH_B})
set_tests_properties(test_tetra_gradients PROPERTIES SKIP_RETURN_CODE 77 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_tetra_mesh tests/test_tetra_mesh.cpp)
target_link_libraries(test_tetra_mesh PRIVATE rfam_core)
add_executable(test_tetra_mesh_geometry tests/test_tetra_mesh_geometry.cpp)
target_link_libraries(test_tetra_mesh_geometry PRIVATE rfam_core)
add_executable(test_tetra_cell_balance tests/test_tetra_cell_balance.cpp)
target_link_libraries(test_tetra_cell_balance PRIVATE rfam_core)
add_executable(test_tetra_freestream tests/test_tetra_freestream.cpp)
target_link_libraries(test_tetra_freestream PRIVATE rfam_core)

# Full-resolution (A) mesh tests. RFAM_TETRA_MESH_A defaults to the delivery mesh
# recorded in docs/PROJECT_TREE.md; it is 145 MB and therefore NOT tracked.
#
# These four targets used to be registered unconditionally against a hard-coded
# path under work/. When that path disappeared (the meshes moved to
# work/archive/old_data/), the tests silently became "cannot open mesh"
# failures: a red run and a never-executed run look the same in a summary, and
# the kernel's principal claims were left unverified without anyone deciding to
# stop verifying them. Registering them only when the mesh is actually present
# keeps the signal honest, and the STATUS line makes the absence visible at
# configure time instead of at test time.
set(RFAM_TETRA_MESH_A "${CMAKE_CURRENT_SOURCE_DIR}/work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet"
    CACHE PATH "Full-resolution tetrahedral delivery mesh used by the mesh-level tests")

if(EXISTS "${RFAM_TETRA_MESH_A}")
    add_test(NAME test_tetra_mesh COMMAND test_tetra_mesh ${RFAM_TETRA_MESH_A})
    set_tests_properties(test_tetra_mesh PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME test_tetra_mesh_geometry COMMAND test_tetra_mesh_geometry ${RFAM_TETRA_MESH_A})
    set_tests_properties(test_tetra_mesh_geometry PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME test_tetra_face_table COMMAND test_tetra_cell_balance ${RFAM_TETRA_MESH_A})
    set_tests_properties(test_tetra_face_table PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME test_tetra_freestream COMMAND test_tetra_freestream ${RFAM_TETRA_MESH_A})
    set_tests_properties(test_tetra_freestream PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    message(STATUS "Tetra A-mesh tests ENABLED  (${RFAM_TETRA_MESH_A})")
else()
    message(STATUS
        "Tetra A-mesh tests SKIPPED: delivery mesh not present at ${RFAM_TETRA_MESH_A}. "
        "4 tests (test_tetra_mesh, _mesh_geometry, _face_table, _freestream) will NOT run. "
        "Restore it, or point -DRFAM_TETRA_MESH_A=<path> at a copy, before reporting that "
        "the tetra kernel is covered.")
endif()
if(RFAM_HAVE_CHEMISTRY)
    target_link_libraries(extract_hotstart_field PRIVATE rfam_chemistry)
endif()

add_executable(test_coolprop_backend tests/thermo/test_coolprop_backend.cpp)
target_link_libraries(test_coolprop_backend PRIVATE rfam_core)
add_test(NAME test_coolprop_backend COMMAND test_coolprop_backend)
add_executable(test_coolprop_transport tests/thermo/test_coolprop_transport.cpp)
target_link_libraries(test_coolprop_transport PRIVATE ${RFAM_COOLPROP_LIBRARY} ${CMAKE_DL_LIBS})
target_include_directories(test_coolprop_transport PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/transport
    ${RFAM_COOLPROP_INCLUDE_DIR} ${RFAM_COOLPROP_FMT_INCLUDE_DIR})
add_test(NAME test_coolprop_transport COMMAND test_coolprop_transport)

add_executable(test_recover_nearcritical tests/thermo/test_recover_nearcritical.cpp)
target_link_libraries(test_recover_nearcritical PRIVATE rfam_core)
add_test(NAME test_recover_nearcritical COMMAND test_recover_nearcritical)
set_tests_properties(test_recover_nearcritical PROPERTIES ENVIRONMENT "RFAM_STATE_ACCEPTANCE_REL=1e-2" SKIP_RETURN_CODE 77)
set_tests_properties(test_coolprop_backend PROPERTIES SKIP_RETURN_CODE 77)

add_executable(audit_tetra_thermo tools/audit_tetra_thermo.cpp)
target_link_libraries(audit_tetra_thermo PRIVATE rfam_core)
add_executable(audit_tetra_linearization tools/audit_tetra_linearization.cpp)
target_link_libraries(audit_tetra_linearization PRIVATE rfam_core)
add_executable(audit_tetra_wall_cell tools/audit_tetra_wall_cell.cpp)
target_link_libraries(audit_tetra_wall_cell PRIVATE rfam_core)
add_executable(audit_tetra_crossing_flux tools/audit_tetra_crossing_flux.cpp)
target_link_libraries(audit_tetra_crossing_flux PRIVATE rfam_core)
add_executable(audit_tetra_core_balance tools/audit_tetra_core_balance.cpp)
target_link_libraries(audit_tetra_core_balance PRIVATE rfam_core)
add_test(NAME test_tetra_binary_linearization COMMAND audit_tetra_linearization --check)
set_tests_properties(test_tetra_binary_linearization PROPERTIES SKIP_RETURN_CODE 77)
add_executable(probe_rhoe_recovery tools/probe_rhoe_recovery.cpp)
target_link_libraries(probe_rhoe_recovery PRIVATE rfam_core)
add_executable(build_property_table tools/build_property_table.cpp)
target_link_libraries(build_property_table PRIVATE rfam_core)
add_executable(bench_recovery tools/bench_recovery.cpp)
target_link_libraries(bench_recovery PRIVATE rfam_core)
if(RFAM_HAVE_CHEMISTRY)
    # Chemistry-only throughput: 6-species two-step against 59-species San Diego.
    add_executable(bench_chemistry tools/bench_chemistry.cpp)
    target_link_libraries(bench_chemistry PRIVATE rfam_chemistry)
    set_target_properties(bench_chemistry PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()
add_executable(check_table_states tools/check_table_states.cpp)
target_link_libraries(check_table_states PRIVATE rfam_core)
add_executable(verify_property_table tools/verify_property_table.cpp)
target_link_libraries(verify_property_table PRIVATE rfam_core)
if(RFAM_HAVE_COOLPROP)
    foreach(tool probe_rhoe_recovery build_property_table verify_property_table check_table_states bench_recovery)
        target_include_directories(${tool} SYSTEM PRIVATE
            ${RFAM_COOLPROP_INCLUDE_DIR} ${RFAM_COOLPROP_FMT_INCLUDE_DIR})
        target_compile_definitions(${tool} PRIVATE RFAM_HAVE_COOLPROP=1)
    endforeach()
else()
    foreach(tool probe_rhoe_recovery build_property_table verify_property_table)
        set_target_properties(${tool} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endforeach()
endif()

add_executable(rfam_v20_tetra apps/solver_3d_tetra.cpp)
target_compile_definitions(rfam_v20_tetra PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}")
if(RFAM_HAVE_COOLPROP)
    target_compile_definitions(rfam_v20_tetra PRIVATE RFAM_HAVE_COOLPROP=1)
    target_include_directories(rfam_v20_tetra SYSTEM PRIVATE
        ${RFAM_COOLPROP_INCLUDE_DIR} ${RFAM_COOLPROP_FMT_INCLUDE_DIR})
endif()
target_link_libraries(rfam_v20_tetra PRIVATE rfam_core)


# Reacting tetrahedral driver (hot-fire path: 30-species San Diego or 2-step CO).
# Same flow/turbulence machinery as rfam_v20_tetra plus finite-rate chemistry.
add_executable(rfam_v20_tetra_react rfam_work/solver_3d_tetra_react.cpp)
target_compile_definitions(rfam_v20_tetra_react PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}")
if(RFAM_HAVE_COOLPROP)
    target_compile_definitions(rfam_v20_tetra_react PRIVATE RFAM_HAVE_COOLPROP=1)
    target_include_directories(rfam_v20_tetra_react SYSTEM PRIVATE
        ${RFAM_COOLPROP_INCLUDE_DIR} ${RFAM_COOLPROP_FMT_INCLUDE_DIR})
endif()
if(RFAM_HAVE_CHEMISTRY)
    target_link_libraries(rfam_v20_tetra_react PRIVATE rfam_chemistry)
else()
    target_link_libraries(rfam_v20_tetra_react PRIVATE rfam_core)
endif()

# ---------------------------------------------------------------------------
# rfam_top: the btop-style run dashboard.
#
# Two build modes because two machines exist: one has ncurses headers, the WSL
# image this repository is developed on has only libncursesw.so.6. In shim mode
# the program still links the real ncurses runtime; only the declarations are
# local (see src/framework/CursesCompat.H).
# ---------------------------------------------------------------------------
find_package(Curses QUIET)
if(CURSES_FOUND)
    add_executable(rfam_top apps/rfam_top.cpp)
    target_compile_definitions(rfam_top PRIVATE RFAM_HAVE_CURSES_HEADERS=1)
    target_include_directories(rfam_top PRIVATE ${CURSES_INCLUDE_DIRS})
    target_link_libraries(rfam_top PRIVATE rfam_core ${CURSES_LIBRARIES})
    message(STATUS "rfam_top: ncurses headers found (${CURSES_LIBRARIES})")
else()
    find_file(RFAM_NCURSES_SO NAMES libncursesw.so.6 libncurses.so.6
        PATHS /usr/lib/x86_64-linux-gnu /lib/x86_64-linux-gnu /usr/lib /lib /usr/lib64)
    find_file(RFAM_TINFO_SO NAMES libtinfo.so.6 libtinfo.so
        PATHS /usr/lib/x86_64-linux-gnu /lib/x86_64-linux-gnu /usr/lib /lib /usr/lib64)
    if(RFAM_NCURSES_SO AND RFAM_TINFO_SO)
        add_executable(rfam_top apps/rfam_top.cpp)
        target_link_libraries(rfam_top PRIVATE rfam_core ${RFAM_NCURSES_SO} ${RFAM_TINFO_SO})
        message(STATUS "rfam_top: no ncurses headers; using the runtime ABI shim (${RFAM_NCURSES_SO})")
    else()
        message(WARNING "rfam_top: ncurses runtime not found; dashboard target skipped")
    endif()
endif()

add_executable(rfam_v20_sector_smoke apps/solver_3d_sector.cpp)
# This driver includes its own sources as "src/...", a spelling that only
# resolves with the repository root on the include path.
target_include_directories(rfam_v20_sector_smoke PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(rfam_v20_sector_smoke PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}"
    RFAM_ALLOW_DIAGNOSTIC_HEM_ACOUSTICS)
if(RFAM_HAVE_CHEMISTRY)
    target_link_libraries(rfam_v20_sector_smoke PRIVATE rfam_chemistry)
else()
    target_link_libraries(rfam_v20_sector_smoke PRIVATE rfam_core)
endif()

add_test(NAME test_sector_reacting_advance_gate
    COMMAND ${CMAKE_COMMAND}
        -DRFAM_SECTOR_EXE=$<TARGET_FILE:rfam_v20_sector_smoke>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/sector_reacting_advance_gate.cmake)

add_test(NAME test_sector_nonreacting_advance_gate
    COMMAND ${CMAKE_COMMAND}
        -DRFAM_SECTOR_EXE=$<TARGET_FILE:rfam_v20_sector_smoke>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/sector_nonreacting_advance_gate.cmake)

add_test(NAME test_sector_zero_step_status
    COMMAND rfam_v20_sector_smoke --steps 0)
set_tests_properties(test_sector_zero_step_status PROPERTIES
    PASS_REGULAR_EXPRESSION "SECTOR_3D_LOAD_OK"
    FAIL_REGULAR_EXPRESSION "SECTOR_3D_STEP_PASS;SECTOR_3D_PASS")

add_test(NAME test_sector_frozen_chemistry_mode_gate
    COMMAND ${CMAKE_COMMAND}
        -DRFAM_SECTOR_EXE=$<TARGET_FILE:rfam_v20_sector_smoke>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/sector_frozen_chemistry_mode_gate.cmake)

add_test(NAME test_sector_steady_local_time_mode_gate
    COMMAND ${CMAKE_COMMAND}
        -DRFAM_SECTOR_EXE=$<TARGET_FILE:rfam_v20_sector_smoke>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/sector_steady_local_time_mode_gate.cmake)

add_test(NAME test_sector_steady_lusgs_mode_gate
    COMMAND ${CMAKE_COMMAND}
        -DRFAM_SECTOR_EXE=$<TARGET_FILE:rfam_v20_sector_smoke>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/sector_steady_lusgs_mode_gate.cmake)

add_executable(rfam_sector_remap_checkpoint tools/remap_sector_checkpoint.cpp)
target_link_libraries(rfam_sector_remap_checkpoint PRIVATE rfam_core)
target_compile_definitions(rfam_sector_remap_checkpoint PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}")

add_executable(rfam_sector_apply_coarse_correction
    tools/apply_sector_coarse_correction.cpp)
target_link_libraries(rfam_sector_apply_coarse_correction PRIVATE rfam_core)
target_compile_definitions(rfam_sector_apply_coarse_correction PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}")

add_executable(rfam_sector_build_fas_source
    tools/build_sector_fas_source.cpp)
target_link_libraries(rfam_sector_build_fas_source PRIVATE rfam_core)
target_compile_definitions(rfam_sector_build_fas_source PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}")

add_executable(rfam_sector_build_galerkin_correction
    tools/build_sector_galerkin_correction.cpp)
target_link_libraries(rfam_sector_build_galerkin_correction PRIVATE rfam_core)
target_compile_definitions(rfam_sector_build_galerkin_correction PRIVATE
    RFAM_GIT_COMMIT="${RFAM_GIT_COMMIT}")

find_program(BASH_EXECUTABLE bash)
if(BASH_EXECUTABLE)
    add_custom_target(regression
        COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/local_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        USES_TERMINAL
    )
    if(RFAM_HAVE_CHEMISTRY)
        add_custom_target(validation_ladder_prephase
            COMMAND ${BASH_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validation_ladder.sh
                ${CMAKE_CURRENT_BINARY_DIR} shock
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            USES_TERMINAL
        )
        add_dependencies(validation_ladder_prephase
            test_pr_thermo_identities
            test_pr_roots
            test_pr_low_pressure_closure
            test_cp_cv_positive
            test_state_conversion_contract
            test_real_fluid_chemistry_energy
            test_phase7_1d
            test_diffusion_real_cv
            test_openmp_consistency
            test_sector_freestream
            test_sector_reacting_freestream
            test_sod_tube)

        add_custom_target(validation_ladder_phase_foundation
            COMMAND ${BASH_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validation_ladder.sh
                ${CMAKE_CURRENT_BINARY_DIR} phase_foundation
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            USES_TERMINAL
        )
        add_dependencies(validation_ladder_phase_foundation
            test_pr_thermo_identities
            test_pr_roots
            test_pr_low_pressure_closure
            test_cp_cv_positive
            test_state_conversion_contract
            test_real_fluid_chemistry_energy
            test_phase7_1d
            test_diffusion_real_cv
            test_openmp_consistency
            test_sector_freestream
            test_sector_reacting_freestream
            test_sod_tube
            test_pr_tp_flash
            test_pr_phase_caloric
            test_pr_phase_thermo_owner_adapter
            test_pr_phase_thermo_owner_displacement
            test_real_fluid_phase_owner_cvode
            test_pr_ph_flash
            test_pr_rho_e_flash
            test_pr_phase_displacement_authority
            test_pr_residual_helmholtz_authority
            test_restart_manifest
            test_restart_manifest_codec
            test_restart_generation_store
            test_restart_atomic_publication_contract
            test_restart_publication
            test_restart_generation_store_sector_clock_preflight
            test_sector_restart_owner_codec
            test_sector_restart_publication_stagers
            test_pr_closed_volume_heating
            test_pr_constant_pressure_heating
            test_pr_two_phase_mixing_contract
            test_pr_vapor_chemistry_reflash_contract
            test_pr_equilibrium_domain
            test_acoustic_state_routing
            test_transport_state_routing
            test_transport_no_alloc
            test_tabulated_ideal_caloric
            test_two_phase_state_contract
            test_model_capability_gate)
    endif()
endif()
