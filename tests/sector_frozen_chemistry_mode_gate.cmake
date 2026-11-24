if(NOT DEFINED RFAM_SECTOR_EXE)
    message(FATAL_ERROR "RFAM_SECTOR_EXE is required")
endif()

# The mode is only meaningful in the 9-species reacting state space.
execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --freeze-chemistry --steps 0
    RESULT_VARIABLE missing_reacting_result
    OUTPUT_VARIABLE missing_reacting_stdout
    ERROR_VARIABLE missing_reacting_stderr
    TIMEOUT 30)
set(missing_reacting "${missing_reacting_stdout}\n${missing_reacting_stderr}")
if(NOT missing_reacting_result EQUAL 2
   OR NOT missing_reacting MATCHES "FREEZE_CHEMISTRY_REQUIRES_REACTING")
    message(FATAL_ERROR
        "freeze chemistry accepted a non-reacting state space\n${missing_reacting}")
endif()

# Zero-step load proves the 9-species state is constructed while no chemistry
# operator is selected. Advancement remains diagnostic-gated separately.
execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --reacting --freeze-chemistry --les --steps 0
    RESULT_VARIABLE load_result
    OUTPUT_VARIABLE load_stdout
    ERROR_VARIABLE load_stderr
    TIMEOUT 30)
set(load_output "${load_stdout}\n${load_stderr}")
if(NOT load_result EQUAL 0)
    message(FATAL_ERROR "frozen-chemistry zero-step load failed\n${load_output}")
endif()
if(NOT load_output MATCHES "Chemistry: FROZEN")
    message(FATAL_ERROR "frozen chemistry marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "Mode: frozen-chemistry-les")
    message(FATAL_ERROR "frozen checkpoint mode missing\n${load_output}")
endif()
if(load_output MATCHES "San Diego 59-species|EXPLICIT substeps|chemistry active/failures")
    message(FATAL_ERROR "a chemistry source operator was selected\n${load_output}")
endif()
if(NOT load_output MATCHES "SECTOR_3D_LOAD_OK")
    message(FATAL_ERROR "frozen chemistry zero-step did not load cleanly\n${load_output}")
endif()

message(STATUS "frozen chemistry mode contract passed")
