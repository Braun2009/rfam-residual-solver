if(NOT RFAM_TETRA_EXE OR NOT RFAM_TETRA_MESH OR NOT RFAM_TETRA_TABLE
   OR NOT EXISTS "${RFAM_TETRA_EXE}" OR NOT EXISTS "${RFAM_TETRA_MESH}"
   OR NOT EXISTS "${RFAM_TETRA_TABLE}")
    message(FATAL_ERROR "tetra checkpoint gate needs executable, mesh and table")
endif()
set(ENV{RFAM_TETRA_TABLE} "${RFAM_TETRA_TABLE}")
unset(ENV{RFAM_STATE_ACCEPTANCE_REL})
unset(ENV{RFAM_PROJECTION_ACCEPT})
unset(ENV{RFAM_SPECIES_REPAIR_REL})
set(ENV{OMP_NUM_THREADS} "8")
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef run_id)
set(base "${CMAKE_CURRENT_BINARY_DIR}/tetra_checkpoint_gate_${run_id}")
set(common --mesh "${RFAM_TETRA_MESH}" --monitor-every 1 --p0 4e6 --T0 400
    --physical-time --physical-dt 2e-7 --physical-inner-iters 100
    --outlet-p 4e6 --backflow-T 400 --backflow-yO2 .784
    --lox-pinj 6e6 --fuel-pinj 4.5e6 --lox-T 130 --fuel-T 450
    --lox-mdot 2.7021 --fuel-mdot 1.0134)
foreach(part continuous first resumed phase_first phase_resumed)
    if(part STREQUAL "continuous")
        set(extra --steps 3 --dump "${base}_continuous.rfamtetf")
    elseif(part STREQUAL "first")
        set(extra --steps 1 --physical-checkpoint "${base}_checkpoint"
            --dump "${base}_first.rfamtetf")
    elseif(part STREQUAL "resumed")
        set(extra --steps 2 --physical-resume "${base}_checkpoint.s1"
            --dump "${base}_resumed.rfamtetf")
    elseif(part STREQUAL "phase_first")
        set(extra --steps 2 --physical-checkpoint "${base}_phase_checkpoint"
            --dump "${base}_phase_first.rfamtetf")
    else()
        set(extra --steps 1 --physical-resume "${base}_phase_checkpoint.s2"
            --dump "${base}_phase_resumed.rfamtetf")
    endif()
    execute_process(COMMAND "${RFAM_TETRA_EXE}" ${common} ${extra}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 300)
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "${part} physical checkpoint run failed: ${result}\n${output}\n${error}")
    endif()
    if(part STREQUAL "resumed" AND NOT output MATCHES
       "PHYSICAL_CHECKPOINT_RESTORED steps=1 physical_time_s=2e-07")
        message(FATAL_ERROR "resume did not restore physical clock: ${output}")
    endif()
    if(part STREQUAL "phase_first" AND NOT output MATCHES "STEP 2[^\n]*two_phase_cells=37")
        message(FATAL_ERROR "two-phase checkpoint fixture did not reach two-phase cells: ${output}")
    endif()
    if(part STREQUAL "phase_resumed" AND NOT output MATCHES
       "PHYSICAL_CHECKPOINT_RESTORED steps=2 physical_time_s=4e-07")
        message(FATAL_ERROR "two-phase resume did not restore physical clock: ${output}")
    endif()
endforeach()
file(SHA256 "${base}_continuous.rfamtetf" continuous_sha)
file(SHA256 "${base}_resumed.rfamtetf" resumed_sha)
file(SHA256 "${base}_phase_resumed.rfamtetf" phase_resumed_sha)
if(NOT continuous_sha STREQUAL resumed_sha OR NOT continuous_sha STREQUAL phase_resumed_sha)
    message(FATAL_ERROR "continuous, single-phase and two-phase resumes differ")
endif()
file(READ "${base}_resumed.rfamtetf.status.json" status)
foreach(expected "\"steps\":3" "\"physical_time_s\":4.4999999999999998e-07"
                 "\"exact_restart\":true" "\"production_validated\":false")
    string(FIND "${status}" "${expected}" at)
    if(at EQUAL -1)
        message(FATAL_ERROR "checkpoint status missing ${expected}: ${status}")
    endif()
endforeach()
execute_process(COMMAND "${RFAM_TETRA_EXE}" ${common} --steps 0
    --fuel-mdot 0.9 --physical-resume "${base}_checkpoint.s1"
    RESULT_VARIABLE rejected OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
if("${rejected}" STREQUAL "0" OR NOT "${stderr}" MATCHES "tetra checkpoint identity mismatch")
    message(FATAL_ERROR "changed inlet flow was not rejected: ${rejected}\n${stdout}\n${stderr}")
endif()
message(STATUS "TETRA_PHYSICAL_CHECKPOINT_GATE_PASS sha256=${continuous_sha}")
