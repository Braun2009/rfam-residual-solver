foreach(required_variable IN ITEMS
    RFAM_R2_PROBE
    RFAM_R2_GENERATOR
    RFAM_R2_REFERENCE
    RFAM_R2_CANTERA_MECHANISM
    RFAM_R2_THERMO
    RFAM_R2_RFAM_MECHANISM
    RFAM_R2_SPECIES_REGISTRY
    RFAM_R2_MASS_ACTION_THERMO
    RFAM_R2_NASA7_THERMO
    RFAM_R2_MECHANISM_READER
    RFAM_R2_THERMO_PROVENANCE_HEADER
    RFAM_R2_THERMO_PROVENANCE_SOURCE
    RFAM_R2_TEST_DIR)
    if(NOT DEFINED ${required_variable})
        message(FATAL_ERROR "${required_variable} is required")
    endif()
endforeach()

foreach(required_file IN ITEMS
    RFAM_R2_PROBE
    RFAM_R2_GENERATOR
    RFAM_R2_REFERENCE
    RFAM_R2_CANTERA_MECHANISM
    RFAM_R2_THERMO
    RFAM_R2_RFAM_MECHANISM
    RFAM_R2_SPECIES_REGISTRY
    RFAM_R2_MASS_ACTION_THERMO
    RFAM_R2_NASA7_THERMO
    RFAM_R2_MECHANISM_READER
    RFAM_R2_THERMO_PROVENANCE_HEADER
    RFAM_R2_THERMO_PROVENANCE_SOURCE)
    if(NOT EXISTS "${${required_file}}")
        message(FATAL_ERROR "${required_file} does not exist")
    endif()
endforeach()

file(REMOVE_RECURSE "${RFAM_R2_TEST_DIR}")
file(MAKE_DIRECTORY "${RFAM_R2_TEST_DIR}")

set(valid_state
    1200.0
    309298.0093953005
    0.9640864999999998
    0.0
    0.6638159542738126
    0.2905351335175837
    0.0456489122086037
    0.0
    0.0)

function(require_rejected case_name expected_pattern)
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${RFAM_R2_TEST_DIR}"
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
    if(combined MATCHES "R2_INSTANTANEOUS_PROBE_OK")
        message(FATAL_ERROR
            "${case_name}: emitted success marker\n${combined}")
    endif()
endfunction()

function(run_rejected case_name expected_pattern generator reference mechanism thermo rfam_mechanism species_registry mass_action_thermo nasa7_thermo mechanism_reader thermo_provenance_header thermo_provenance_source)
    require_rejected(
        "${case_name}"
        "${expected_pattern}"
        "${RFAM_R2_PROBE}"
        "${generator}"
        "${reference}"
        "${mechanism}"
        "${thermo}"
        "${rfam_mechanism}"
        "${species_registry}"
        "${mass_action_thermo}"
        "${nasa7_thermo}"
        "${mechanism_reader}"
        "${thermo_provenance_header}"
        "${thermo_provenance_source}"
        ${valid_state})
endfunction()

require_rejected(
    no_arguments
    "provenance-arguments-required"
    "${RFAM_R2_PROBE}")

set(labels generator reference cantera-mechanism thermo rfam-mechanism species-registry mass-action-thermo nasa7-thermo mechanism-reader thermo-provenance-header thermo-provenance-source)
set(paths
    "${RFAM_R2_GENERATOR}"
    "${RFAM_R2_REFERENCE}"
    "${RFAM_R2_CANTERA_MECHANISM}"
    "${RFAM_R2_THERMO}"
    "${RFAM_R2_RFAM_MECHANISM}"
    "${RFAM_R2_SPECIES_REGISTRY}"
    "${RFAM_R2_MASS_ACTION_THERMO}"
    "${RFAM_R2_NASA7_THERMO}"
    "${RFAM_R2_MECHANISM_READER}"
    "${RFAM_R2_THERMO_PROVENANCE_HEADER}"
    "${RFAM_R2_THERMO_PROVENANCE_SOURCE}")

foreach(index RANGE 0 10)
    list(GET labels ${index} label)
    list(GET paths ${index} original)
    set(missing "${RFAM_R2_TEST_DIR}/missing-${label}")
    set(arguments ${paths})
    list(REMOVE_AT arguments ${index})
    list(INSERT arguments ${index} "${missing}")
    list(GET arguments 0 generator)
    list(GET arguments 1 reference)
    list(GET arguments 2 mechanism)
    list(GET arguments 3 thermo)
    list(GET arguments 4 rfam_mechanism)
    list(GET arguments 5 species_registry)
    list(GET arguments 6 mass_action_thermo)
    list(GET arguments 7 nasa7_thermo)
    list(GET arguments 8 mechanism_reader)
    list(GET arguments 9 thermo_provenance_header)
    list(GET arguments 10 thermo_provenance_source)
    run_rejected(
        "missing_${label}"
        "provenance-open-failed"
        "${generator}"
        "${reference}"
        "${mechanism}"
        "${thermo}"
        "${rfam_mechanism}"
        "${species_registry}"
        "${mass_action_thermo}"
        "${nasa7_thermo}"
        "${mechanism_reader}"
        "${thermo_provenance_header}"
        "${thermo_provenance_source}")

    file(READ "${original}" contents)
    set(stale "${RFAM_R2_TEST_DIR}/stale-${label}")
    file(WRITE "${stale}" "${contents}\n# stale\n")
    set(arguments ${paths})
    list(REMOVE_AT arguments ${index})
    list(INSERT arguments ${index} "${stale}")
    list(GET arguments 0 generator)
    list(GET arguments 1 reference)
    list(GET arguments 2 mechanism)
    list(GET arguments 3 thermo)
    list(GET arguments 4 rfam_mechanism)
    list(GET arguments 5 species_registry)
    list(GET arguments 6 mass_action_thermo)
    list(GET arguments 7 nasa7_thermo)
    list(GET arguments 8 mechanism_reader)
    list(GET arguments 9 thermo_provenance_header)
    list(GET arguments 10 thermo_provenance_source)
    run_rejected(
        "stale_${label}"
        "${label}-sha256-mismatch"
        "${generator}"
        "${reference}"
        "${mechanism}"
        "${thermo}"
        "${rfam_mechanism}"
        "${species_registry}"
        "${mass_action_thermo}"
        "${nasa7_thermo}"
        "${mechanism_reader}"
        "${thermo_provenance_header}"
        "${thermo_provenance_source}")
endforeach()

require_rejected(
    invalid_state
    "mass-fraction-sum-mismatch"
    "${RFAM_R2_PROBE}"
    "${RFAM_R2_GENERATOR}"
    "${RFAM_R2_REFERENCE}"
    "${RFAM_R2_CANTERA_MECHANISM}"
    "${RFAM_R2_THERMO}"
    "${RFAM_R2_RFAM_MECHANISM}"
    "${RFAM_R2_SPECIES_REGISTRY}"
    "${RFAM_R2_MASS_ACTION_THERMO}"
    "${RFAM_R2_NASA7_THERMO}"
    "${RFAM_R2_MECHANISM_READER}"
    "${RFAM_R2_THERMO_PROVENANCE_HEADER}"
    "${RFAM_R2_THERMO_PROVENANCE_SOURCE}"
    1200.0 309298.0 0.9640865
    0.0 0.5 0.2 0.04 0.0 0.0)

message("R2_INSTANTANEOUS_PROVENANCE_GATE_PASS")
