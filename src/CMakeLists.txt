# Kokkos requires 3.16
CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

# Project TChem_atm
PROJECT(TChem_atm)

# Built-in modules
INCLUDE(GNUInstallDirs)
INCLUDE(CMakePackageConfigHelpers)

# CMAKE setup
SET(CMAKE_VERBOSE_MAKEFILE ON) # see all output
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

# Define BUILD_TESTING in cmake
INCLUDE(CTest)

# Installation
IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}" CACHE PATH "FOO install prefix" FORCE)
ENDIF()

IF (NOT CMAKE_INSTALL_LIBDIR)
  SET(CMAKE_INSTALL_LIBDIR lib)
ENDIF()
MESSAGE(STATUS "CMAKE_INSTALL_LIBDIR is set ${CMAKE_INSTALL_LIBDIR}")

SET(TCHEM_ATM_INSTALL_LIB_PATH ${CMAKE_INSTALL_LIBDIR})
SET(TCHEM_ATM_INSTALL_INCLUDE_PATH include/tchem_atm)
SET(TCHEM_ATM_INSTALL_BIN_PATH     bin)
SET(TCHEM_ATM_INSTALL_TEST_PATH    unit-tests)
SET(TCHEM_ATM_INSTALL_EXAMPLE_PATH examples)

LIST(APPEND CMAKE_MODULE_PATH "${TChem_atm_SOURCE_DIR}/cmake/modules")
MESSAGE(STATUS "TChem_ATM uses CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}")
#
# Options : use TCHEM_ATM prefix
#
OPTION(TCHEM_ATM_ENABLE_TEST "Flag to enable unit tests" OFF)
OPTION(TCHEM_ATM_ENABLE_REAL_TYPE "Flag to determine working precision; float, double(default)" "double")
OPTION(TCHEM_ATM_ENABLE_EXAMPLE "Flag to enable unit examples" ON)
OPTION(TCHEM_ATM_ENABLE_VERBOSE "Flag to enable TChem_atm verbose flag" OFF)
OPTION(TCHEM_ATM_ENABLE_DEBUG "Flag to enable TChem_atm debug flag" OFF)
OPTION(TCHEM_ATM_ENABLE_INSTALL "Flag to enable install exports for TCHEM_ATM" ON)
OPTION(TCHEM_ATM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY "Flag to enable sacado jacobian for atmospheric chemistry" OFF)
OPTION(TCHEM_ATM_ENABLE_SKYWALKER "Flag to enable Skywalker" OFF)
OPTION(TCHEM_ATM_ENABLE_KOKKOSKERNELS "Flag to enable Kokkos-kernels" OFF)
OPTION(TCHEM_ATM_ENABLE_COVERAGE  "Enable code coverage instrumentation" OFF)

OPTION(KOKKOS_INSTALL_PATH "Path to Kokkos installation")
OPTION(TINES_INSTALL_PATH "Path to Tines installation")
OPTION(GTEST_INSTALL_PATH "Path to gtest installation")
OPTION(SKYWALKER_INSTALL_PATH "Path to Skywalker installation")
OPTION(SUNDIALS_INSTALL_PATH "Path to SUNDIALS installation")

INCLUDE_DIRECTORIES(${TChem_atm_BINARY_DIR})

# Working precision
IF (TCHEM_ATM_ENABLE_REAL_TYPE STREQUAL "float")
  SET(TCHEM_ATM_ENABLE_REAL_TYPE_SINGLE_PRECISION ON)
  MESSAGE(STATUS "TChem_atm uses real type SINGLE PRECISION")
ELSE()
  SET(TCHEM_ATM_ENABLE_REAL_TYPE_DOUBLE_PRECISION ON)
  MESSAGE(STATUS "TChem_atm uses real type DOUBLE PRECISION")
ENDIF()

# TChem_atm link libraries for tests
SET(TCHEM_ATM_LINK_LIBRARIES tchem_atm)

#
# Kokkos
#
IF(TARGET Kokkos::kokkos)
  # check to see if kokkos was loaded in another project
  message(STATUS "Found Kokkos::kokkos target")
  SET(TCHEM_ATM_INTERNAL_KOKKOS_TARGET "Kokkos::kokkos")
  SET(TCHEM_ATM_ENABLE_TPL_KOKKOS ON CACHE BOOL "Enable Kokkos")
  LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_KOKKOS_TARGET})
ELSE()
  IF (KOKKOS_INSTALL_PATH)
    MESSAGE(STATUS "TChem_atm attempts to use Kokkos installed at ${KOKKOS_INSTALL_PATH}")
    FIND_PACKAGE(Kokkos REQUIRED)
    SET(TCHEM_ATM_INTERNAL_KOKKOS_TARGET "Kokkos::kokkos")
    SET(TCHEM_ATM_ENABLE_TPL_KOKKOS ON CACHE BOOL "Enable Kokkos")
    LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_KOKKOS_TARGET})
  ELSE()
    MESSAGE(FATAL_ERROR "-- Kokkos install path is not provided")
  ENDIF()
ENDIF()

# Kokkos-kernels
#
IF(TARGET Kokkos::kokkoskernels)
  # check to see if kokkos-kernels was loaded in another project
  message(STATUS "Found Kokkos::kokkoskernels target")
  SET(TCHEM_ATM_INTERNAL_KOKKOSKERNELS_TARGET "Kokkos::kokkoskernels")
  SET(TCHEM_ATM_ENABLE_TPL_KOKKOSKERNELS ON CACHE BOOL "Enable Kokkos-kernels")
  LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_KOKKOSKERNELS_TARGET})
ELSE()
  IF (TCHEM_ATM_ENABLE_KOKKOSKERNELS)
  IF (KOKKOSKERNELS_INSTALL_PATH)
    MESSAGE(STATUS "TChem_atm attempts to use Kokkos-kernels installed at ${KOKKOSKERNELS_INSTALL_PATH}")
    FIND_PACKAGE(kokkoskernels REQUIRED)
    SET(TCHEM_ATM_INTERNAL_KOKKOSKERNELS_TARGET "Kokkos::kokkoskernels")
    SET(TCHEM_ATM_ENABLE_TPL_KOKKOSKERNELS ON CACHE BOOL "Enable Kokkos-kernels")
    LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_KOKKOSKERNELS_TARGET})
  ELSE()
    MESSAGE(FATAL_ERROR "-- Kokkos-kernels install path is not provided")
  ENDIF()
  ENDIF()
ENDIF()

#
# Tines
#
IF(TARGET Tines::tines)
  # check to see if tines was loaded in another project
  message(STATUS "Found Tines::tines target")
  SET(TCHEM_ATM_INTERNAL_TINES_TARGET "Tines::tines")
  SET(TCHEM_ATM_ENABLE_TPL_TINES ON CACHE BOOL "Enable Tines")
  LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_TINES_TARGET})
ELSE()
  IF (TINES_INSTALL_PATH)
    MESSAGE(STATUS "TChem_atm attempts to use Tines installed at ${TINES_INSTALL_PATH}")
    FIND_PACKAGE(Tines REQUIRED)
    SET(TCHEM_ATM_INTERNAL_TINES_TARGET "Tines::tines")
    SET(TCHEM_ATM_ENABLE_TPL_TINES ON CACHE BOOL "Enable Tines")
    LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_TINES_TARGET})
  ELSE()
    MESSAGE(FATAL_ERROR "-- Tines install path is not provided")
  ENDIF()
ENDIF()


IF (SUNDIALS_INSTALL_PATH)
  # Interface installed SUNDIALS
  FIND_PACKAGE(SUNDIALS)
  SET(TCHEM_ATM_INTERNAL_SUNDIALS_TARGET "SUNDIALS::cvode")
  SET(TCHEM_ATM_ENABLE_TPL_SUNDIALS ON CACHE BOOL "Enable Sundials")
  LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_SUNDIALS_TARGET})
ELSEIF(TCHEM_ATM_SUNDIALS_WARNING)
  MESSAGE(WARNING "-- SUNDIALS_INSTALL_PATH is not provided; SUNDIALS TPL is not enabled")
ENDIF()

#
# gtest
#
IF (TCHEM_ATM_ENABLE_TEST)
  IF (GTEST_INSTALL_PATH)
    # Interface installed gtest
    MESSAGE("-- TChem_atm uses gtest installed at ${GTEST_INSTALL_PATH}")
    FIND_PACKAGE(Gtest REQUIRED)
    SET(TCHEM_ATM_INTERNAL_GTEST_TARGET "GTest::gtest_main")
  ELSE()
    MESSAGE(FATAL_ERROR "-- Test is enabled but Gtest install path is not provided")
  ENDIF()
  LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_GTEST_TARGET})
ENDIF()

## skywalker

if (TCHEM_ATM_ENABLE_SKYWALKER)
IF(TARGET skywalker)
# check to see if tines was loaded in another project
  message(STATUS "Found skywalker target")
  SET(TCHEM_ATM_INTERNAL_SKYWALKER_TARGET "skywalker")
  SET(TCHEM_ATM_ENABLE_TPL_SKYWALKER ON CACHE BOOL "Enable skywalker")
  LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_ENABLE_TPL_SKYWALKER})
ELSE()
IF (SKYWALKER_INSTALL_PATH)
    MESSAGE(STATUS "TChem_atm attempts to use Skywalker installed at ${SKYWALKER_INSTALL_PATH}")
    FIND_PACKAGE(skywalker REQUIRED)
    SET(TCHEM_ATM_INTERNAL_SKYWALKER_TARGET "skywalker")
    SET(TCHEM_ATM_ENABLE_TPL_SKYWALKER ON CACHE BOOL "Enable skywalker")
    LIST(APPEND TCHEM_ATM_LINK_LIBRARIES ${TCHEM_ATM_INTERNAL_SKYWALKER_TARGET})
  ELSE()
    MESSAGE(FATAL_ERROR "-- Skywalker install path is not provided")
  ENDIF()

ENDIF()
ENDIF()


# Code coverage

if (TCHEM_ATM_ENABLE_COVERAGE)
  MESSAGE(STATUS "Enabling code coverage instrumentation")
  FIND_PROGRAM(LCOV_EXE lcov DOC "Lcov code coverage tool")
  if (LCOV_EXE MATCHES "NOTFOUND")
    MESSAGE(FATAL_ERROR "Could not find lcov for code coverage reporting!")
  endif()

  # Add code coverage compiler/linker flags
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

  # Add a "make coverage" target.
  add_custom_target(coverage ${LCOV_EXE} --capture --directory . -o unfiltered_coverage.info
    COMMAND ${LCOV_EXE} --remove unfiltered_coverage.info -o coverage.info '*/external/*' '*/tests/*' '*/validation/*'
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating coverage report (coverage.info)")
endif()


#
# Linking libraries used in example and test
#
MESSAGE("-- TCHEM_ATM_LINK_LIBRARIES = ${TCHEM_ATM_LINK_LIBRARIES}")

#
# Define configure-related preprocessor macros
#
CONFIGURE_FILE (
  ${TChem_atm_SOURCE_DIR}/core/TChem_atm_Config.hpp.in
  ${TChem_atm_BINARY_DIR}/TChem_atm_Config.hpp
  )
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/TChem_atm_Config.hpp
  DESTINATION "${CMAKE_INSTALL_PREFIX}/${TCHEM_ATM_INSTALL_INCLUDE_PATH}")

#
# Process subdirectories
#
ADD_SUBDIRECTORY (core)
IF (TCHEM_ATM_ENABLE_TEST)
  ADD_SUBDIRECTORY (unit-tests)
ENDIF()
IF (TCHEM_ATM_ENABLE_EXAMPLE)
  ADD_SUBDIRECTORY (examples)
ENDIF()
IF (TCHEM_ATM_ENABLE_SKYWALKER)
  ADD_SUBDIRECTORY (verification)
ENDIF()

# add in target include for build to enable finding TChem_atm_Config
TARGET_INCLUDE_DIRECTORIES(tchem_atm
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
