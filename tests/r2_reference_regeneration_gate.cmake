foreach(required_variable IN ITEMS
    RFAM_PYTHON
    RFAM_R2_GENERATOR
    RFAM_R2_REFERENCE
    RFAM_R2_CANTERA_MECHANISM
    RFAM_R2_THERMO
    RFAM_R2_SPECIES_REGISTRY
    RFAM_R2_REGEN_DIR)
    if(NOT DEFINED ${required_variable})
        message(FATAL_ERROR "${required_variable} is required")
    endif()
endforeach()

file(REMOVE_RECURSE "${RFAM_R2_REGEN_DIR}")
file(MAKE_DIRECTORY "${RFAM_R2_REGEN_DIR}")
set(generated "${RFAM_R2_REGEN_DIR}/r2_instantaneous_reference_v1.json")

execute_process(
    COMMAND
        "${RFAM_PYTHON}"
        "${RFAM_R2_GENERATOR}"
        --mechanism "${RFAM_R2_CANTERA_MECHANISM}"
        --thermo "${RFAM_R2_THERMO}"
        --species-registry "${RFAM_R2_SPECIES_REGISTRY}"
        --output "${generated}"
    RESULT_VARIABLE generator_result
    OUTPUT_VARIABLE generator_stdout
    ERROR_VARIABLE generator_stderr)
if(NOT generator_result EQUAL 0)
    message(FATAL_ERROR
        "R2 reference regeneration failed\n${generator_stdout}\n${generator_stderr}")
endif()
if(NOT generator_stdout MATCHES "R2_REFERENCE_WRITTEN cases=12")
    message(FATAL_ERROR
        "R2 reference regeneration marker missing\n${generator_stdout}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${generated}" "${RFAM_R2_REFERENCE}"
    RESULT_VARIABLE compare_result)
if(NOT compare_result EQUAL 0)
    file(SHA256 "${generated}" generated_sha256)
    file(SHA256 "${RFAM_R2_REFERENCE}" reference_sha256)
    message(FATAL_ERROR
        "R2 reference regeneration mismatch: generated=${generated_sha256} reference=${reference_sha256}")
endif()

message("R2_REFERENCE_REGENERATION_GATE_PASS")
