# Check Tines installation
FILE(GLOB_RECURSE TINES_FOUND_CMAKE_FILE "${TINES_INSTALL_PATH}/TinesConfig.cmake")

MESSAGE(STATUS "Tines install path : ${TINES_INSTALL_PATH}")
MESSAGE(STATUS "Tines found cmake config : ${TINES_FOUND_CMAKE_FILE}")

IF (TINES_FOUND_CMAKE_FILE) 
  INCLUDE(${TINES_FOUND_CMAKE_FILE})
  SET(TINES_FOUND ON)
ELSE()
  MESSAGE(FATAL_ERROR "-- Tines is not found at ${TINES_INSTALL_PATH}")
ENDIF()