# to handle GPU compiling
SET(FILE_NAME TChem_Driver.c)
string(FIND "${CMAKE_CXX_COMPILER}" "nvcc" FOUND_NVCC)
string(FIND "${CMAKE_CXX_COMPILER}" "hipcc" FOUND_HIPCC)
IF(NOT FOUND_NVCC EQUAL -1 OR NOT FOUND_HIPCC EQUAL -1)
  FILE(COPY TChem_Driver.c DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  FILE(RENAME ${CMAKE_CURRENT_BINARY_DIR}/TChem_Driver.c ${CMAKE_CURRENT_BINARY_DIR}/TChem_Driver.cpp)
  SET(FILE_NAME TChem_Driver.cpp)
ENDIF()

SET_SOURCE_FILES_PROPERTIES(${FILE_NAME} PROPERTIES LANGUAGE CXX)
ADD_EXECUTABLE(TChem_Driver.x ${FILE_NAME})
TARGET_LINK_LIBRARIES(TChem_Driver.x ${TCHEM_ATM_LINK_LIBRARIES})
INSTALL(TARGETS TChem_Driver.x)

FILE(COPY config_solver.yaml
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the 'references' directory to the current binary directory
FILE(COPY references
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Set the test name and executable name
SET(TEST_NAME DriverInterface)
SET(EXE_NAME ./TChem_Driver.x)

# Set the input path for the example run
SET(INPUT_PATH ../../examples/runs/atmospheric_chemistry/CB05CL_AE5_w_simpolSOA)

# Set the validation script path
FILE(COPY compare_files.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
SET(VALIDATION_SCRIPT ./compare_files.py)

# Add a test to execute the example code
ADD_TEST(run_${TEST_NAME} ${EXE_NAME}
         ${INPUT_PATH}/config_full_gas.yaml
         ${INPUT_PATH}/mechanism_aero.yaml 
         ./config_solver.yaml)

# Add a test to validate the output data with the reference data
ADD_TEST(validate_${TEST_NAME} python3 ${VALIDATION_SCRIPT})
