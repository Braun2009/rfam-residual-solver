if(NOT DEFINED RFAM_SECTOR_EXE)
    message(FATAL_ERROR "RFAM_SECTOR_EXE is required")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --steady-local-time --steady-lusgs --steps 0
    RESULT_VARIABLE conflict_result
    OUTPUT_VARIABLE conflict_stdout
    ERROR_VARIABLE conflict_stderr
    TIMEOUT 30)
set(conflict_output "${conflict_stdout}\n${conflict_stderr}")
if(NOT conflict_result EQUAL 2
   OR NOT conflict_output MATCHES "TIME_INTEGRATOR_CONFLICT")
    message(FATAL_ERROR
        "steady integrator conflict was not rejected\n${conflict_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --lusgs-anderson --steps 0
    RESULT_VARIABLE anderson_result
    OUTPUT_VARIABLE anderson_stdout
    ERROR_VARIABLE anderson_stderr
    TIMEOUT 30)
set(anderson_output "${anderson_stdout}\n${anderson_stderr}")
if(NOT anderson_result EQUAL 2
   OR NOT anderson_output MATCHES "LUSGS_ANDERSON_REQUIRES_LUSGS")
    message(FATAL_ERROR
        "orphan LU-SGS Anderson mode was not rejected\n${anderson_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --lusgs-frozen-euler-block --steps 0
    RESULT_VARIABLE block_result
    OUTPUT_VARIABLE block_stdout
    ERROR_VARIABLE block_stderr
    TIMEOUT 30)
set(block_output "${block_stdout}\n${block_stderr}")
if(NOT block_result EQUAL 2
   OR NOT block_output MATCHES "LUSGS_BLOCK_REQUIRES_LUSGS")
    message(FATAL_ERROR
        "orphan LU-SGS block mode was not rejected\n${block_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --lusgs-hydro-diagonal-block --steps 0
    RESULT_VARIABLE hydro_block_result
    OUTPUT_VARIABLE hydro_block_stdout
    ERROR_VARIABLE hydro_block_stderr
    TIMEOUT 30)
set(hydro_block_output "${hydro_block_stdout}\n${hydro_block_stderr}")
if(NOT hydro_block_result EQUAL 2
   OR NOT hydro_block_output MATCHES "LUSGS_HYDRO_BLOCK_REQUIRES_LUSGS")
    message(FATAL_ERROR
        "orphan LU-SGS hydro block mode was not rejected\n${hydro_block_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --lusgs-pressure-schur-least-squares-gradient --steps 0
    RESULT_VARIABLE schur_ls_result
    OUTPUT_VARIABLE schur_ls_stdout
    ERROR_VARIABLE schur_ls_stderr
    TIMEOUT 30)
set(schur_ls_output "${schur_ls_stdout}\n${schur_ls_stderr}")
if(NOT schur_ls_result EQUAL 2
   OR NOT schur_ls_output MATCHES
       "LUSGS_PRESSURE_SCHUR_LS_REQUIRES_PRESSURE_SCHUR")
    message(FATAL_ERROR
        "orphan pressure Schur LS gradient was not rejected\n${schur_ls_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --steady-lusgs --lusgs-pressure-schur-fallback
        --lusgs-pressure-schur-nonorthogonal-correction --steps 0
    RESULT_VARIABLE schur_nonorth_result
    OUTPUT_VARIABLE schur_nonorth_stdout
    ERROR_VARIABLE schur_nonorth_stderr
    TIMEOUT 30)
set(schur_nonorth_output
    "${schur_nonorth_stdout}\n${schur_nonorth_stderr}")
if(NOT schur_nonorth_result EQUAL 2
   OR NOT schur_nonorth_output MATCHES
       "LUSGS_PRESSURE_SCHUR_NONORTH_REQUIRES_LS")
    message(FATAL_ERROR
        "orphan pressure Schur nonorthogonal correction was not rejected\n${schur_nonorth_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --lusgs-pressure-schur-structured-momentum --steps 0
    RESULT_VARIABLE schur_structured_result
    OUTPUT_VARIABLE schur_structured_stdout
    ERROR_VARIABLE schur_structured_stderr
    TIMEOUT 30)
set(schur_structured_output
    "${schur_structured_stdout}\n${schur_structured_stderr}")
if(NOT schur_structured_result EQUAL 2
   OR NOT schur_structured_output MATCHES
       "LUSGS_PRESSURE_SCHUR_STRUCTURED_REQUIRES_PRESSURE_SCHUR")
    message(FATAL_ERROR
        "orphan structured pressure Schur was not rejected\n${schur_structured_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --steady-lusgs --lusgs-pressure-schur-fallback
        --lusgs-pressure-schur-least-squares-gradient
        --lusgs-pressure-schur-structured-momentum --steps 0
    RESULT_VARIABLE schur_variant_conflict_result
    OUTPUT_VARIABLE schur_variant_conflict_stdout
    ERROR_VARIABLE schur_variant_conflict_stderr
    TIMEOUT 30)
set(schur_variant_conflict_output
    "${schur_variant_conflict_stdout}\n${schur_variant_conflict_stderr}")
if(NOT schur_variant_conflict_result EQUAL 2
   OR NOT schur_variant_conflict_output MATCHES
       "LUSGS_PRESSURE_SCHUR_VARIANT_CONFLICT")
    message(FATAL_ERROR
        "pressure Schur variant conflict was not rejected\n${schur_variant_conflict_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --steady-lusgs --lusgs-segregated-energy-fallback --steps 0
    RESULT_VARIABLE energy_fallback_result
    OUTPUT_VARIABLE energy_fallback_stdout
    ERROR_VARIABLE energy_fallback_stderr
    TIMEOUT 30)
set(energy_fallback_output
    "${energy_fallback_stdout}\n${energy_fallback_stderr}")
if(NOT energy_fallback_result EQUAL 2
   OR NOT energy_fallback_output MATCHES
       "LUSGS_SEGREGATED_FALLBACK_REQUIRES_BIDIRECTIONAL_LINE_SEARCH")
    message(FATAL_ERROR
        "orphan LU-SGS energy fallback was not rejected\n${energy_fallback_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --lusgs-residual-line-search --steps 0
    RESULT_VARIABLE line_search_result
    OUTPUT_VARIABLE line_search_stdout
    ERROR_VARIABLE line_search_stderr
    TIMEOUT 30)
set(line_search_output "${line_search_stdout}\n${line_search_stderr}")
if(NOT line_search_result EQUAL 2
   OR NOT line_search_output MATCHES "LUSGS_LINE_SEARCH_REQUIRES_LUSGS")
    message(FATAL_ERROR
        "orphan LU-SGS residual line search was not rejected\n${line_search_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --fas-source missing.bin --steps 0
    RESULT_VARIABLE fas_mode_result
    OUTPUT_VARIABLE fas_mode_stdout
    ERROR_VARIABLE fas_mode_stderr
    TIMEOUT 30)
set(fas_mode_output "${fas_mode_stdout}\n${fas_mode_stderr}")
if(NOT fas_mode_result EQUAL 2
   OR NOT fas_mode_output MATCHES "FAS_SOURCE_MODE_MISMATCH")
    message(FATAL_ERROR
        "orphan FAS source mode was not rejected\n${fas_mode_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --grid coarse --steady-lusgs --reacting --freeze-chemistry --les
        --fas-source missing.bin --steps 0
    RESULT_VARIABLE fas_restart_result
    OUTPUT_VARIABLE fas_restart_stdout
    ERROR_VARIABLE fas_restart_stderr
    TIMEOUT 30)
set(fas_restart_output "${fas_restart_stdout}\n${fas_restart_stderr}")
if(NOT fas_restart_result EQUAL 2
   OR NOT fas_restart_output MATCHES "FAS_SOURCE_MODE_MISMATCH")
    message(FATAL_ERROR
        "FAS source without restart was not rejected\n${fas_restart_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --steady-lusgs --lusgs-bidirectional-line-search --steps 0
    RESULT_VARIABLE bidirectional_result
    OUTPUT_VARIABLE bidirectional_stdout
    ERROR_VARIABLE bidirectional_stderr
    TIMEOUT 30)
set(bidirectional_output
    "${bidirectional_stdout}\n${bidirectional_stderr}")
if(NOT bidirectional_result EQUAL 2
   OR NOT bidirectional_output MATCHES
       "LUSGS_BIDIRECTIONAL_REQUIRES_LINE_SEARCH")
    message(FATAL_ERROR
        "orphan bidirectional LU-SGS mode was not rejected\n${bidirectional_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}" --reacting --steady-lusgs --steps 0
    RESULT_VARIABLE reacting_result
    OUTPUT_VARIABLE reacting_stdout
    ERROR_VARIABLE reacting_stderr
    TIMEOUT 30)
set(reacting_output "${reacting_stdout}\n${reacting_stderr}")
if(NOT reacting_result EQUAL 2
   OR NOT reacting_output MATCHES
       "STEADY_REACTING_REQUIRES_FROZEN_CHEMISTRY")
    message(FATAL_ERROR
        "finite-rate reacting LU-SGS was not rejected\n${reacting_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --reacting --freeze-chemistry --les --steady-lusgs
        --lusgs-anderson --lusgs-hydro-diagonal-block
        --lusgs-residual-line-search --lusgs-bidirectional-line-search
        --lusgs-pressure-schur-fallback
        --lusgs-pressure-schur-least-squares-gradient
        --lusgs-pressure-schur-nonorthogonal-correction
        --lusgs-segregated-acoustic-fallback
        --lusgs-segregated-flow-fallback
        --lusgs-segregated-energy-fallback --steps 0
    RESULT_VARIABLE load_result
    OUTPUT_VARIABLE load_stdout
    ERROR_VARIABLE load_stderr
    TIMEOUT 30)
set(load_output "${load_stdout}\n${load_stderr}")
if(NOT load_result EQUAL 0)
    message(FATAL_ERROR "steady LU-SGS frozen zero-step load failed\n${load_output}")
endif()
if(NOT load_output MATCHES "STEADY LU-SGS")
    message(FATAL_ERROR "steady LU-SGS marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "Anderson ON")
    message(FATAL_ERROR "steady LU-SGS Anderson marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "hydro diagonal block ON")
    message(FATAL_ERROR "steady LU-SGS hydro block marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "pressure Schur fallback ON")
    message(FATAL_ERROR "steady LU-SGS pressure Schur marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "pressure Schur least-squares gradient ON")
    message(FATAL_ERROR
        "steady LU-SGS pressure Schur LS marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "pressure Schur nonorthogonal correction ON")
    message(FATAL_ERROR
        "steady LU-SGS pressure Schur nonorthogonal marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "segregated acoustic fallback ON")
    message(FATAL_ERROR "steady LU-SGS acoustic fallback marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "segregated flow fallback ON")
    message(FATAL_ERROR "steady LU-SGS flow fallback marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "segregated energy fallback ON")
    message(FATAL_ERROR "steady LU-SGS energy fallback marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "residual line search ON")
    message(FATAL_ERROR "steady LU-SGS residual line-search marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "bidirectional ON")
    message(FATAL_ERROR "steady LU-SGS bidirectional marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "physical clock frozen")
    message(FATAL_ERROR "LU-SGS physical-clock contract marker missing\n${load_output}")
endif()
if(NOT load_output MATCHES "Mode: frozen-chemistry-les")
    message(FATAL_ERROR "LU-SGS frozen checkpoint mode missing\n${load_output}")
endif()
if(NOT load_output MATCHES "SECTOR_3D_LOAD_OK")
    message(FATAL_ERROR "LU-SGS frozen zero-step did not load cleanly\n${load_output}")
endif()

execute_process(
    COMMAND "${RFAM_SECTOR_EXE}"
        --reacting --freeze-chemistry --les --steady-lusgs
        --lusgs-residual-line-search --lusgs-bidirectional-line-search
        --lusgs-pressure-schur-fallback
        --lusgs-pressure-schur-structured-momentum --steps 0
    RESULT_VARIABLE structured_load_result
    OUTPUT_VARIABLE structured_load_stdout
    ERROR_VARIABLE structured_load_stderr
    TIMEOUT 30)
set(structured_load_output
    "${structured_load_stdout}\n${structured_load_stderr}")
if(NOT structured_load_result EQUAL 0
   OR NOT structured_load_output MATCHES
       "pressure Schur structured momentum ON"
   OR NOT structured_load_output MATCHES "SECTOR_3D_LOAD_OK")
    message(FATAL_ERROR
        "structured pressure Schur zero-step load failed\n${structured_load_output}")
endif()

message(STATUS "steady LU-SGS mode contract passed")
