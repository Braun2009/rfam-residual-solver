# A two-phase state must fail the single-phase molecular transport gate.
if(NOT RFAM_TETRA_EXE OR NOT RFAM_TETRA_MESH OR NOT RFAM_TETRA_TABLE
   OR NOT EXISTS "${RFAM_TETRA_TABLE}" OR NOT EXISTS "${RFAM_TETRA_MESH}")
    message(FATAL_ERROR "tetra molecular phase gate needs executable, mesh and existing table")
endif()
set(ENV{RFAM_TETRA_TABLE} "${RFAM_TETRA_TABLE}")
set(ENV{RFAM_STATE_ACCEPTANCE_REL} "1e-2")
set(ENV{RFAM_PROJECTION_ACCEPT} "1")
set(ENV{OMP_NUM_THREADS} "8")
execute_process(
    COMMAND "${RFAM_TETRA_EXE}" --mesh "${RFAM_TETRA_MESH}"
        --steps 0 --p0 4e6 --T0 170 --yO2 0.9 --yC3H8 0.1
        --lox-pinj 6e6 --fuel-pinj 4.5e6 --lox-T 130 --fuel-T 450
        --lox-mdot 2.8 --fuel-mdot 1.05 --inlet-supply-cap --inlet-face-supply-cap
        --molecular-viscous --lusgs
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error
    TIMEOUT 120)
if("${result}" STREQUAL "0" OR NOT "${output}${error}" MATCHES
        "TETRA_COLD_FLOW_FAILED: molecular-viscous two-phase closure not admitted")
    message(FATAL_ERROR "two-phase transport did not fail closed: result=${result}\n${output}\n${error}")
endif()
message(STATUS "TETRA_MOLECULAR_PHASE_GATE_PASS: exit=${result}")
