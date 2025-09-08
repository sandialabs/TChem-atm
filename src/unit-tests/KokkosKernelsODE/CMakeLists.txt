# Copy the 'references' directory to the current binary directory
FILE(COPY references
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Set the test name and executable name
SET(TEST_NAME AerosolChemistry_KokkosKernels)
SET(EXE_NAME ../../examples/TChem_${TEST_NAME}.x)

# Set the input path for the example run
SET(INPUT_PATH ../../examples/runs/atmospheric_chemistry/Simpol)

# Set the validation script path
SET(VALIDATION_SCRIPT ../scripts/compare_tchem.py)

# Threshold error for the passing test
SET(PASSING_ERROR 1e-8)

# Set the output file name
SET(OUTPUT_NAME_TEST full_gas.dat)
SET(OUTPUT_NAME_REF full_gas.hdf5)

# Add a test to execute the example code
ADD_TEST(run_${TEST_NAME} ${EXE_NAME} 
    --chemfile=${INPUT_PATH}/config_gas.yaml 
    --aerofile=${INPUT_PATH}/test_SIMPOL_phase_transfer.yaml 
    --inputfile_particles=${INPUT_PATH}/scenario_conditions_particle.yaml 
    --outputfile=${OUTPUT_NAME_TEST} 
    --outputfile_times=${TEST_NAME}_wall_times.json
    --dtmin=1 
    --tend=5)

# Add a test to validate the output data with the reference data
ADD_TEST(validate_${TEST_NAME} python3 ${VALIDATION_SCRIPT} 
    -r_file=references/${OUTPUT_NAME_REF}
    -t_file=${OUTPUT_NAME_TEST}
    -check_norms=True
    -error=${PASSING_ERROR})


## Unit test with batched test

# Set the output file name
SET(TEST_NAME batched_AerosolChemistry_KokkosKernels)
SET(OUTPUT_NAME_TEST batched_full_gas.dat)
SET(OUTPUT_NAME_REF batched_full_gas.hdf5)

# Add a test to execute the example code
ADD_TEST(run_${TEST_NAME} ${EXE_NAME} 
    --chemfile=${INPUT_PATH}/config_gas.yaml 
    --aerofile=${INPUT_PATH}/test_SIMPOL_phase_transfer.yaml 
    --inputfile_particles=${INPUT_PATH}/scenario_conditions_particle.yaml 
    --outputfile=${OUTPUT_NAME_TEST} 
    --outputfile_times=${TEST_NAME}_wall_times.json
    --dtmin=1 
    --tend=5
    --use_cloned_samples=true 
    --batch_size=3 )

# Add a test to validate the output data with the reference data
ADD_TEST(validate_${TEST_NAME} python3 ${VALIDATION_SCRIPT} 
    -r_file=references/${OUTPUT_NAME_REF}
    -t_file=${OUTPUT_NAME_TEST}
    -check_norms=True
    -error=${PASSING_ERROR})


# 
# Set the output file name
SET(TEST_NAME batched_CB5SimpolAerosolChemistry_KokkosKernels)
SET(OUTPUT_NAME_TEST batched_cb_full_gas.dat)
SET(OUTPUT_NAME_REF batched_cb_full_gas.hdf5)

# Set the input path for the example run
SET(INPUT_PATH ../../examples/runs/atmospheric_chemistry/CB05CL_AE5_w_simpolSOA)


# Add a test to execute the example code
ADD_TEST(run_${TEST_NAME} ${EXE_NAME} 
    --chemfile=${INPUT_PATH}/config_full_gas.yaml 
    --aerofile=${INPUT_PATH}/mechanism_aero.yaml 
    --inputfile_particles=${INPUT_PATH}/scenario_conditions_particle.yaml 
    --outputfile=${OUTPUT_NAME_TEST} 
    --outputfile_times=${TEST_NAME}_wall_times.json
    --dtmin=1e-3 
    --tend=5
    --max-time-iterations=2
    --use_cloned_samples=true 
    --batch_size=3 )

# Add a test to validate the output data with the reference data
ADD_TEST(validate_${TEST_NAME} python3 ${VALIDATION_SCRIPT} 
    -r_file=references/${OUTPUT_NAME_REF}
    -t_file=${OUTPUT_NAME_TEST}
    -check_norms=True
    -error=${PASSING_ERROR})

