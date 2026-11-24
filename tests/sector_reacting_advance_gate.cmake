if(NOT DEFINED RFAM_SECTOR_EXE)
    message(FATAL_ERROR "RFAM_SECTOR_EXE is required")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --rans --reacting --steps 1
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
    TIMEOUT 30)

set(combined "${stdout}\n${stderr}")
if(NOT result EQUAL 3)
    message(FATAL_ERROR
        "reacting advance gate returned ${result}, expected 3\n${combined}")
endif()
if(NOT combined MATCHES "REACTING_ADVANCE_BLOCKED")
    message(FATAL_ERROR
        "reacting advance gate did not emit its block marker\n${combined}")
endif()
if(combined MATCHES "SECTOR_3D_STEP_PASS|SECTOR_3D_PASS")
    message(FATAL_ERROR
        "blocked reacting run emitted an advancement pass marker\n${combined}")
endif()

message(STATUS "reacting advance blocked before solver construction")
