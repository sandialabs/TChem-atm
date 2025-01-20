#!/bin/bash

#=======================================================================================
# example script to download, build, install, & test:
#     tines, tchem
#=======================================================================================
# User configuration  -- begin

# example:
# set CUDA="ON" to use GPU, else set to "OFF"
# JFLAG is for makefile compilation on multiple cores, here 4
# if CUDA="ON", make sure nvcc is in your system PATH

MY_CC=gcc
MY_CXX=g++
MY_FC=gfortran
JFLAG="-j 40"
CUDA="OFF"
SACADO="OFF"

## update this path.
ROOT=/path/to/tchem
INSTALL_BASE_HOST=$ROOT/build_tchem/HOST/install
BUILD_TYPE=DEBUG
#BUILD_TYPE=RELEASE

# User configuration  -- end

#=======================================================================================
build_install_tines(){
echo "Building tines:"
mkdir ${TINES_BUILD_PATH}
mkdir ${TINES_INSTALL_PATH}
cd ${TINES_BUILD_PATH}
#    -D TINES_ENABLE_SHARED_BUILD=ON \
cmake \
    -D CMAKE_INSTALL_PREFIX="${TINES_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${KOKKOS_CXX_COMPILER}" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_BUILD_TYPE=$BUILD_TYPE \
    -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
    -D TINES_ENABLE_DEBUG=OFF \
    -D TINES_ENABLE_VERBOSE=OFF \
    -D TINES_ENABLE_TEST=ON \
    -D TINES_ENABLE_EXAMPLE=ON \
    -D SUNDIALS_INSTALL_PATH="${SUNDIALS_INSTALL_PATH}" \
    -D YAML_INSTALL_PATH="${YAML_INSTALL_PATH}" \
    -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
    -D GTEST_INSTALL_PATH="${GTEST_INSTALL_PATH}" \
    -D OPENBLAS_INSTALL_PATH="${OPENBLAS_INSTALL_PATH}" \
    ${TINES_REPOSITORY_PATH}/src
make ${JFLAG} install
# if using e.g. yum or apt-get installed openblas, may need this:
# Not needed when using above installed openblas and macport version on osx
#    -D LAPACKE_INSTALL_PATH="${LAPACKE_INSTALL_PATH}" \
}
build_install_tchem_atm(){
echo "Building TChem_atm:"
mkdir ${TCHEM_BUILD_PATH}
mkdir ${TCHEM_INSTALL_PATH}
cd ${TCHEM_BUILD_PATH}
#    -D TCHEM_ENABLE_SHARED_BUILD=ON \
cmake \
    -D CMAKE_INSTALL_PREFIX="${TCHEM_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${KOKKOS_CXX_COMPILER}" \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D TCHEM_ATM_ENABLE_VERBOSE=OFF \
    -D CMAKE_BUILD_TYPE=$BUILD_TYPE \
    -D TCHEM_ATM_ENABLE_TEST=ON \
    -D TCHEM_ATM_ENABLE_REAL_TYPE="double" \
    -D TCHEM_ATM_ENABLE_EXAMPLE=ON \
    -D TCHEM_ATM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY=${SACADO}\
    -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
    -D KOKKOSKERNELS_INSTALL_PATH="${KOKKOSKERNELS_INSTALL_PATH}" \
    -D TINES_INSTALL_PATH="${TINES_INSTALL_PATH}" \
    -D TCHEM_ATM_ENABLE_KOKKOSKERNELS=ON \
    -D TCHEM_ATM_ENABLE_SKYWALKER=ON \
    -D SKYWALKER_INSTALL_PATH=${SKYWALKER_INSTALL_PATH} \
    -D GTEST_INSTALL_PATH="${GTEST_INSTALL_PATH}" \
    ${TCHEM_REPOSITORY_PATH}/src
make ${JFLAG} install
}
#=======================================================================================
#=======================================================================================
# main

## tpls; assumes tpls_bls.sh was executed first.
SKYWALKER_INSTALL_PATH=${INSTALL_BASE_HOST}/skywalker
OPENBLAS_INSTALL_PATH=${INSTALL_BASE_HOST}/openblas
GTEST_INSTALL_PATH=${INSTALL_BASE_HOST}/gtest
YAML_INSTALL_PATH=${INSTALL_BASE_HOST}/yaml
SUNDIALS_INSTALL_PATH=${INSTALL_BASE_HOST}/sundials
KOKKOSKERNELS_INSTALL_PATH=${INSTALL_BASE_HOST}/kokkos-kernels

if [ "${CUDA}" = "ON" ]; then
  BUILD_BASE=${PWD}/DEVICE/$BUILD_TYPE/build
  INSTALL_BASE=${PWD}/DEVICE/$BUILD_TYPE/install
  INSTALL_BASE_TPL_DEVICE=${PWD}/$BUILD_TYPE/install
  KOKKOS_REPOSITORY_PATH=${TINES_REPOSITORY_PATH}/ext/kokkos
  KOKKOS_INSTALL_PATH=${INSTALL_BASE_TPL_DEVICE}/kokkos
  KOKKOS_CXX_COMPILER="${KOKKOS_INSTALL_PATH}/bin/nvcc_wrapper"
  KOKKOS_CXX_REPO_COMPILER="${KOKKOS_REPOSITORY_PATH}/bin/nvcc_wrapper"
  KOKKOSKERNELS_INSTALL_PATH=$INSTALL_BASE_TPL_DEVICE/kokkos-kernels
  #NOTE! we do not use sundials in GPUs.
  SUNDIALS_INSTALL_PATH=""
else
  BUILD_BASE=${PWD}/HOST/$BUILD_TYPE/build
  INSTALL_BASE=${PWD}/HOST/$BUILD_TYPE/install
  KOKKOS_INSTALL_PATH=${INSTALL_BASE_HOST}/kokkos
  KOKKOS_CXX_COMPILER="${MY_CXX}"
  KOKKOS_CXX_REPO_COMPILER="${MY_CXX}"
fi

mkdir -p ${BUILD_BASE}
mkdir -p ${INSTALL_BASE}

TCHEM_REPOSITORY_PATH=$ROOT/TChem-atm

TINES_REPOSITORY_PATH=$TCHEM_REPOSITORY_PATH/external/Tines
TINES_BUILD_PATH=${BUILD_BASE}/tines
TINES_INSTALL_PATH=${INSTALL_BASE}/tines
build_install_tines

TCHEM_BUILD_PATH=${BUILD_BASE}/tchem_atm
TCHEM_INSTALL_PATH=${INSTALL_BASE}/tchem_atm
build_install_tchem_atm

exit
