if(NOT DEFINED RFAM_GATE_EXE OR NOT EXISTS "${RFAM_GATE_EXE}")
    message(FATAL_ERROR "RFAM_GATE_EXE is required")
endif()
if(NOT DEFINED RFAM_GENERATOR OR NOT EXISTS "${RFAM_GENERATOR}")
    message(FATAL_ERROR "RFAM_GENERATOR is required")
endif()
if(NOT DEFINED RFAM_REFERENCE OR NOT EXISTS "${RFAM_REFERENCE}")
    message(FATAL_ERROR "RFAM_REFERENCE is required")
endif()
if(NOT DEFINED RFAM_MECHANISM OR NOT EXISTS "${RFAM_MECHANISM}")
    message(FATAL_ERROR "RFAM_MECHANISM is required")
endif()
if(NOT DEFINED RFAM_TEST_DIR)
    message(FATAL_ERROR "RFAM_TEST_DIR is required")
endif()

file(REMOVE_RECURSE "${RFAM_TEST_DIR}")
file(MAKE_DIRECTORY "${RFAM_TEST_DIR}")

function(require_rejected case_name expected_pattern)
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${RFAM_TEST_DIR}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE stdout
        ERROR_VARIABLE stderr)
    set(combined "${stdout}\n${stderr}")
    if(result EQUAL 0)
        message(FATAL_ERROR
            "${case_name}: unexpectedly passed\n${combined}")
    endif()
    if(NOT combined MATCHES "${expected_pattern}")
        message(FATAL_ERROR
            "${case_name}: missing '${expected_pattern}'\n${combined}")
    endif()
    if(combined MATCHES "CANTERA_KINETICS_GATE_BLOCKED")
        message(FATAL_ERROR
            "${case_name}: emitted blocked success marker\n${combined}")
    endif()
endfunction()

require_rejected(
    no_arguments
    "provenance-arguments-required"
    "${RFAM_GATE_EXE}")

require_rejected(
    missing_generator
    "provenance-open-failed"
    "${RFAM_GATE_EXE}"
    "${RFAM_TEST_DIR}/missing-generator.py"
    "${RFAM_REFERENCE}"
    "${RFAM_MECHANISM}")

require_rejected(
    missing_reference
    "provenance-open-failed"
    "${RFAM_GATE_EXE}"
    "${RFAM_GENERATOR}"
    "${RFAM_TEST_DIR}/missing-reference.txt"
    "${RFAM_MECHANISM}")

require_rejected(
    missing_mechanism
    "provenance-open-failed"
    "${RFAM_GATE_EXE}"
    "${RFAM_GENERATOR}"
    "${RFAM_REFERENCE}"
    "${RFAM_TEST_DIR}/missing-mechanism.yaml")

file(READ "${RFAM_GENERATOR}" generator_contents)
file(WRITE "${RFAM_TEST_DIR}/stale-generator.py"
    "${generator_contents}\n# stale\n")
require_rejected(
    stale_generator
    "generator-sha256-mismatch"
    "${RFAM_GATE_EXE}"
    "${RFAM_TEST_DIR}/stale-generator.py"
    "${RFAM_REFERENCE}"
    "${RFAM_MECHANISM}")

file(READ "${RFAM_REFERENCE}" reference_contents)
file(WRITE "${RFAM_TEST_DIR}/stale-reference.txt"
    "${reference_contents}\n# stale\n")
require_rejected(
    stale_reference
    "reference-sha256-mismatch"
    "${RFAM_GATE_EXE}"
    "${RFAM_GENERATOR}"
    "${RFAM_TEST_DIR}/stale-reference.txt"
    "${RFAM_MECHANISM}")

file(READ "${RFAM_MECHANISM}" mechanism_contents)
file(WRITE "${RFAM_TEST_DIR}/stale-mechanism.yaml"
    "${mechanism_contents}\n# stale\n")
require_rejected(
    stale_mechanism
    "mechanism-sha256-mismatch"
    "${RFAM_GATE_EXE}"
    "${RFAM_GENERATOR}"
    "${RFAM_REFERENCE}"
    "${RFAM_TEST_DIR}/stale-mechanism.yaml")

message("CANTERA_KINETICS_PROVENANCE_GATE_PASS")
