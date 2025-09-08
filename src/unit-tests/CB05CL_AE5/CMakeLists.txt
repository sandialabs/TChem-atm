# Copy the 'references' directory to the current binary directory
FILE(COPY references
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Set the test name and executable name
SET(TEST_NAME AtmosphericChemistry)
SET(EXE_NAME ../../examples/TChem_${TEST_NAME}.x)

# Set the input path for the example run
SET(INPUT_PATH ../../examples/runs/atmospheric_chemistry/CB05CL_AE5)

# Set the validation script path
SET(VALIDATION_SCRIPT ../scripts/compare_tchem.py)

# Threshold error for the passing test
SET(PASSING_ERROR 1e-8)

# Set the output file name
SET(OUTPUT_NAME_TEST full_gas.dat)
SET(OUTPUT_NAME_REF full_gas.hdf5)

# Add a test to execute the example code
ADD_TEST(run_CB05CL_AE5_${TEST_NAME} ${EXE_NAME}
    --chemfile=${INPUT_PATH}/config_full_gas.yaml
    --outputfile=${OUTPUT_NAME_TEST}
    --outputfile_times=${TEST_NAME}_wall_times.json
    --time-iterations-per-interval=10
    --dtmin=1e-20
    --dtmax=10
    --tend=600
    --atol-newton=1e-18
    --rtol-newton=1e-8
    --max-newton-iterations=20
    --max-time-iterations=2
    )

# Add a test to validate the output data with the reference data
ADD_TEST(validate_CB05CL_AE5_${TEST_NAME} python3 ${VALIDATION_SCRIPT}
    -r_file=references/${OUTPUT_NAME_REF}
    -t_file=${OUTPUT_NAME_TEST}
    -check_norms=True
    -error=${PASSING_ERROR})
