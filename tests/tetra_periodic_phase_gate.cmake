# A periodic dump must carry same-step EOS phase, boundary budget and status.
if(NOT RFAM_TETRA_EXE OR NOT RFAM_TETRA_MESH OR NOT RFAM_TETRA_TABLE
   OR NOT EXISTS "${RFAM_TETRA_TABLE}" OR NOT EXISTS "${RFAM_TETRA_MESH}")
    message(FATAL_ERROR "tetra periodic phase gate needs executable, mesh and table")
endif()
set(ENV{RFAM_TETRA_TABLE} "${RFAM_TETRA_TABLE}")
set(ENV{RFAM_STATE_ACCEPTANCE_REL} "1e-2")
set(ENV{RFAM_PROJECTION_ACCEPT} "1")
set(ENV{RFAM_SPECIES_REPAIR_REL} "1e-5")
set(ENV{OMP_NUM_THREADS} "8")
set(base "${CMAKE_CURRENT_BINARY_DIR}/tetra_periodic_phase_gate.rfamtetf")
execute_process(
    COMMAND "${RFAM_TETRA_EXE}" --mesh "${RFAM_TETRA_MESH}"
        --steps 1 --monitor-every 1 --dump-every 1 --dump-phase-every
        --p0 4.11e6 --T0 400 --yO2 0.784 --yC3H8 0.216
        --outlet-p 4.1116e6 --backflow-T 400 --backflow-yO2 0.784
        --lox-pinj 6e6 --fuel-pinj 4.5e6 --lox-T 130 --fuel-T 450
        --lox-mdot 2.8 --fuel-mdot 1.05 --inlet-supply-cap
        --inlet-face-supply-cap --lusgs --lusgs-dt-scale 1 --lusgs-relax 0.2
        --dump "${base}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 160)
if(NOT "${result}" STREQUAL "0" OR NOT "${output}" MATCHES
       "TETRA_COLD_FLOW_ITERATION_LIMIT accepted_steps=1")
    message(FATAL_ERROR "periodic dump failed: result=${result}\n${output}\n${error}")
endif()
foreach(ext "" ".phase.csv" ".boundary.csv")
    if(NOT EXISTS "${base}${ext}" OR NOT EXISTS "${base}.s1${ext}")
        message(FATAL_ERROR "missing paired periodic artifact ${ext}")
    endif()
    file(SHA256 "${base}${ext}" old_hash)
    file(SHA256 "${base}.s1${ext}" step_hash)
    if(NOT "${old_hash}" STREQUAL "${step_hash}")
        message(FATAL_ERROR "periodic artifact differs from same-step final ${ext}")
    endif()
endforeach()
if(NOT EXISTS "${base}.s1.status.json")
    message(FATAL_ERROR "missing periodic status")
endif()
file(READ "${base}.s1.status.json" status)
foreach(expected "\"steps\":1" "\"phase_sidecar_available\":true"
                 "\"boundary_budget_available\":true" "\"converged\":false"
                 "\"production_validated\":false" "\"projection_accept\":true"
                 "\"state_acceptance_rel\":0.01")
    string(FIND "${status}" "${expected}" index)
    if(index EQUAL -1)
        message(FATAL_ERROR "periodic status missing ${expected}: ${status}")
    endif()
endforeach()
message(STATUS "TETRA_PERIODIC_PHASE_GATE_PASS")
