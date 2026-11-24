if(NOT DEFINED RFAM_SECTOR_EXE)
    message(FATAL_ERROR "RFAM_SECTOR_EXE is required")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --rans --steps 1
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
    TIMEOUT 30)

set(combined "${stdout}\n${stderr}")
if(NOT result EQUAL 4)
    message(FATAL_ERROR
        "nonreacting advance gate returned ${result}, expected 4\n${combined}")
endif()
if(NOT combined MATCHES "SECTOR_ADVANCE_BLOCKED")
    message(FATAL_ERROR
        "nonreacting advance gate did not emit its block marker\n${combined}")
endif()
if(combined MATCHES "SECTOR_3D_STEP_PASS|SECTOR_3D_PASS")
    message(FATAL_ERROR
        "blocked sector run emitted an advancement pass marker\n${combined}")
endif()

message(STATUS "sector advance blocked before solver construction")
