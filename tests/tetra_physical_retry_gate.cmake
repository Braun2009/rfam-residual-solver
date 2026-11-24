if(NOT RFAM_TETRA_EXE OR NOT RFAM_TETRA_MESH OR NOT RFAM_TETRA_TABLE
   OR NOT EXISTS "${RFAM_TETRA_TABLE}" OR NOT EXISTS "${RFAM_TETRA_MESH}")
    message(FATAL_ERROR "physical retry gate needs executable, mesh and table")
endif()
set(ENV{RFAM_TETRA_TABLE} "${RFAM_TETRA_TABLE}")
unset(ENV{RFAM_STATE_ACCEPTANCE_REL})
unset(ENV{RFAM_PROJECTION_ACCEPT})
unset(ENV{RFAM_SPECIES_REPAIR_REL})
set(ENV{OMP_NUM_THREADS} "8")
set(base "${CMAKE_CURRENT_BINARY_DIR}/tetra_physical_retry_gate.rfamtetf")
execute_process(
    COMMAND "${RFAM_TETRA_EXE}" --mesh "${RFAM_TETRA_MESH}"
        --steps 3 --monitor-every 1 --p0 4e6 --T0 400 --physical-time
        --physical-dt 2e-7 --outlet-p 4e6 --backflow-T 400
        --backflow-yO2 .784 --lox-pinj 6e6 --fuel-pinj 4.5e6
        --lox-T 130 --fuel-T 450 --lox-mdot 2.7021 --fuel-mdot 1.0134
        --physical-inner-iters 100 --dump "${base}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 300)
if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "physical retry run failed: result=${result}\n${output}\n${error}")
endif()
set(log "${output}\n${error}")
string(REGEX MATCHALL "PHYSICAL_STEP_RETRY step=3 retry=[12] " retries "${log}")
list(LENGTH retries retry_count)
if(NOT retry_count EQUAL 2)
    message(FATAL_ERROR "expected two trial-boundary retries, got ${retry_count}:\n${log}")
endif()
if(NOT log MATCHES "STEP 3[^\n]*physical_time_s=4\\.5e-07 physical_dt_s=5e-08[^\n]*physical_step_retries=2 physical_boundary_trial_retries=2")
    message(FATAL_ERROR "accepted step did not report halved dt, unchanged rollback clock, and retry counts:\n${log}")
endif()
if(NOT EXISTS "${base}.status.json")
    message(FATAL_ERROR "missing physical retry status JSON")
endif()
file(READ "${base}.status.json" status)
foreach(expected "\"steps\":3" "\"physical_time_s\":4." 
                 "\"physical_step_retries\":2" 
                 "\"physical_boundary_trial_retries\":2" 
                 "\"physical_retry_count_total\":2" 
                 "\"physical_boundary_trial_retry_count_total\":2" 
                 "\"physical_inner_converged\":true" 
                 "\"production_validated\":false")
    string(FIND "${status}" "${expected}" index)
    if(index EQUAL -1)
        message(FATAL_ERROR "physical retry status missing ${expected}: ${status}")
    endif()
endforeach()
message(STATUS "TETRA_PHYSICAL_RETRY_GATE_PASS")

set(initial_base "${CMAKE_CURRENT_BINARY_DIR}/tetra_physical_initial_boundary_reject.rfamtetf")
execute_process(
    COMMAND "${RFAM_TETRA_EXE}" --mesh "${RFAM_TETRA_MESH}"
        --steps 0 --monitor-every 1 --p0 4e6 --T0 400 --physical-time
        --physical-dt 2e-7 --outlet-p 4e6 --backflow-T 400
        --backflow-yO2 .784 --lox-pinj 6e6 --fuel-pinj 4.05e6
        --lox-T 130 --fuel-T 450 --lox-mdot 2.7021 --fuel-mdot 1.0134
        --dump "${initial_base}"
    RESULT_VARIABLE initial_result OUTPUT_VARIABLE initial_output
    ERROR_VARIABLE initial_error TIMEOUT 120)
if("${initial_result}" STREQUAL "0")
    message(FATAL_ERROR "incompatible initial boundary unexpectedly succeeded:\n${initial_output}\n${initial_error}")
endif()
set(initial_log "${initial_output}\n${initial_error}")
if(NOT initial_log MATCHES "PHYSICAL_BOUNDARY_INPUT_REJECTED initial")
    message(FATAL_ERROR "initial boundary rejection was not classified explicitly:\n${initial_log}")
endif()
if(NOT EXISTS "${initial_base}.status.json")
    message(FATAL_ERROR "missing initial boundary rejection status JSON")
endif()
file(READ "${initial_base}.status.json" initial_status)
foreach(expected "\"steps\":0" "\"physical_time_s\":0" 
                 "\"physical_boundary_input_rejected\":true"
                 "\"boundary_state_valid\":false" "\"update_rejected\":true")
    string(FIND "${initial_status}" "${expected}" index)
    if(index EQUAL -1)
        message(FATAL_ERROR "initial boundary status missing ${expected}: ${initial_status}")
    endif()
endforeach()
message(STATUS "TETRA_PHYSICAL_INITIAL_BOUNDARY_REJECT_PASS")
