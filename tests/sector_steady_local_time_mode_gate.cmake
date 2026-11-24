if(NOT DEFINED RFAM_SECTOR_EXE)
    message(FATAL_ERROR "RFAM_SECTOR_EXE is required")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --dual-time --steady-local-time --steps 0
    RESULT_VARIABLE conflict_result
    OUTPUT_VARIABLE conflict_stdout
    ERROR_VARIABLE conflict_stderr
    TIMEOUT 30)
set(conflict_output "${conflict_stdout}\n${conflict_stderr}")
if(NOT conflict_result EQUAL 2
   OR NOT conflict_output MATCHES "TIME_INTEGRATOR_CONFLICT")
    message(FATAL_ERROR
        "dual/steady integrator conflict was not rejected\n${conflict_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --reacting --steady-local-time --steps 0
    RESULT_VARIABLE reacting_result
    OUTPUT_VARIABLE reacting_stdout
    ERROR_VARIABLE reacting_stderr
    TIMEOUT 30)
set(reacting_output "${reacting_stdout}\n${reacting_stderr}")
if(NOT reacting_result EQUAL 2
   OR NOT reacting_output MATCHES
       "STEADY_REACTING_REQUIRES_FROZEN_CHEMISTRY")
    message(FATAL_ERROR
        "finite-rate reacting steady mode was not rejected\n${reacting_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --reacting --freeze-chemistry --les --steady-local-time --steps 0
    RESULT_VARIABLE load_result
    OUTPUT_VARIABLE load_stdout
    ERROR_VARIABLE load_stderr
    TIMEOUT 30)
set(load_output "${load_stdout}\n${load_stderr}")
if(NOT load_result EQUAL 0)
    message(FATAL_ERROR "steady frozen zero-step load failed\n${load_output}")
endif()
if(NOT load_output MATCHES "STEADY local pseudo-time")
    message(FATAL_ERROR "steady local-time marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "physical clock frozen")
    message(FATAL_ERROR "steady physical-clock contract marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "Mode: frozen-chemistry-les")
    message(FATAL_ERROR "steady frozen checkpoint mode missing\n${load_output}")
endif()
if(NOT load_output MATCHES "SECTOR_3D_LOAD_OK")
    message(FATAL_ERROR "steady frozen zero-step did not load cleanly\n${load_output}")
endif()

message(STATUS "steady local-time mode contract passed")
