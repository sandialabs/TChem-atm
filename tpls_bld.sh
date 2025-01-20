#!/bin/bash

#=======================================================================================
# example script to download, build, install, & test:
#     openblas, kokkos, gtest, sundials, skywalker
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

# will build under BUILD_BASE
# will install under INSTALL_BASE
# example: as follows under .
ROOT=/path/to/tchem/
# User configuration  -- end
#=======================================================================================
# OpenBLAS
# nb. to make sure this openblas gets used when running outside this script
# add in your .bashrc/.bash_profile :
# on mac need this
#    export LIBRARY_PATH="${LIBRARY_PATH}:${OPENBLAS_INSTALL_PATH}/lib"
# on linux need this
#    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${OPENBLAS_INSTALL_PATH}/lib"

get_submodules() {
  echo "Getting submodules:"
  run_this="cd $TCHEM_BASE;git submodule update --init --recursive;cd -"
  eval $run_this
}

build_openblas (){
echo "Building OpenBLAS:"
cd ${OPENBLAS_REPOSITORY_PATH}
make CC=${MY_CC} FC=${MY_FC} HOSTCC=${MY_CC} USE_OPENMP=1
}

install_openblas (){
echo "Installing OpenBLAS:"
cd ${OPENBLAS_REPOSITORY_PATH}
make PREFIX=${OPENBLAS_INSTALL_PATH} install
}

build_install_kokkoskernels(){
echo "Building kokkos kernels:"
mkdir ${KOKKOSKERNELS_BUILD_PATH}
mkdir ${KOKKOSKERNELS_INSTALL_PATH}
cd ${KOKKOSKERNELS_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX="${KOKKOSKERNELS_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${KOKKOS_CXX_COMPILER}" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
    -D CMAKE_BUILD_TYPE=$BUILD_TYPE \
    -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
    -D Kokkos_ROOT=$KOKKOS_INSTALL_PATH \
    ${KOKKOSKERNELS_REPOSITORY_PATH}
make ${JFLAG} install
}
build_install_kokkos(){
echo "Building kokkos:"
mkdir ${KOKKOS_BUILD_PATH}
mkdir ${KOKKOS_INSTALL_PATH}
cd ${KOKKOS_BUILD_PATH}
# kokkos install bin directory is not yet created
cmake \
    -D CMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER=${KOKKOS_CXX_REPO_COMPILER} \
    -D CMAKE_CXX_FLAGS="-fopenmp -g" \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_CUDA=${CUDA} \
    -D Kokkos_ENABLE_CUDA_CONSTEXPR=${CUDA} \
    -D Kokkos_ENABLE_CUDA_LAMBDA=${CUDA} \
    ${KOKKOS_REPOSITORY_PATH}
make ${JFLAG} install
}
#=======================================================================================

build_install_gtest(){
echo "Building gtest:"
mkdir ${GTEST_BUILD_PATH}
mkdir ${GTEST_INSTALL_PATH}
cd ${GTEST_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX="${GTEST_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    ${GTEST_REPOSITORY_PATH}
make ${JFLAG} install
}

build_install_yaml(){
echo "Building gtest:"
mkdir ${YAML_BUILD_PATH}
mkdir ${YAML_INSTALL_PATH}
cd ${YAML_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX="${YAML_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_CXX_FLAGS="-g -c" \
    -D CMAKE_EXE_LINKER_FLAGS="" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    ${YAML_REPOSITORY_PATH}
make ${JFLAG} install
}
get_sundials (){
echo "get sundials:"
if [ -d "${SUNDIALS_REPOSITORY_PATH}" ] && [ "$(ls -A ${SUNDIALS_REPOSITORY_PATH})" ]; then
  echo "${SUNDIALS_REPOSITORY_PATH} exists and is not empty ... aborting clone"; return
fi
git clone https://github.com/LLNL/sundials.git ${SUNDIALS_REPOSITORY_PATH}
}

build_install_sundials(){
echo "Building gtest:"
mkdir ${SUNDIALS_BUILD_PATH}
mkdir ${SUNDIALS_INSTALL_PATH}
cd ${SUNDIALS_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D CMAKE_C_FLAGS="-g" \
    -D CMAKE_EXE_LINKER_FLAGS="" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D ENABLE_CALIPER:BOOL=OFF \
    ${SUNDIALS_REPOSITORY_PATH}
make ${JFLAG} install
}

build_install_skywalker(){
echo "Building skywalker:"
mkdir ${SKYWALKER_BUILD_PATH}
mkdir ${SKYWALKER_INSTALL_PATH}
cd ${SKYWALKER_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX=${SKYWALKER_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D SKYWALKER_PRECISION=double \
    -D CMAKE_BUILD_TYPE=RELEASE \
    ${SKYWALKER_REPOSITORY_PATH}
make ${JFLAG} install
}
#=======================================================================================

#=======================================================================================
# main

TCHEM_BASE=$ROOT/TChem-atm/external
REPO_BASE=$TCHEM_BASE/Tines/ext

if [ "${CUDA}" = "ON" ]; then
    BUILD_BASE=${PWD}/DEVICE/build
    INSTALL_BASE=${PWD}/DEVICE/install
else
    BUILD_BASE=${PWD}/HOST/build
    INSTALL_BASE=${PWD}/HOST/install
fi

mkdir -p ${BUILD_BASE}
mkdir -p ${INSTALL_BASE}
# Note: git submodule update --init --recursive
get_submodules

#only build for host
if [ "${CUDA}" = "OFF" ]; then
  # clone tpls
  OPENBLAS_REPOSITORY_PATH=${REPO_BASE}/OpenBLAS
  OPENBLAS_INSTALL_PATH=${INSTALL_BASE}/openblas
  build_openblas
  install_openblas

  GTEST_REPOSITORY_PATH=${REPO_BASE}/gtest
  GTEST_BUILD_PATH=${BUILD_BASE}/gtest
  GTEST_INSTALL_PATH=${INSTALL_BASE}/gtest
  build_install_gtest
  #
  YAML_REPOSITORY_PATH=${REPO_BASE}/yaml
  YAML_BUILD_PATH=${BUILD_BASE}/yaml
  YAML_INSTALL_PATH=${INSTALL_BASE}/yaml
  build_install_yaml

  SUNDIALS_REPOSITORY_PATH=$TCHEM_BASE/Sundials
  SUNDIALS_BUILD_PATH=${BUILD_BASE}/sundials
  SUNDIALS_INSTALL_PATH=${INSTALL_BASE}/sundials
  build_install_sundials

  SKYWALKER_REPOSITORY_PATH=$TCHEM_BASE/Skywalker
  SKYWALKER_BUILD_PATH=${BUILD_BASE}/skywalker
  SKYWALKER_INSTALL_PATH=${INSTALL_BASE}/skywalker
  build_install_skywalker
fi

KOKKOS_REPOSITORY_PATH=${REPO_BASE}/kokkos
KOKKOS_BUILD_PATH=${BUILD_BASE}/kokkos
KOKKOS_INSTALL_PATH=${INSTALL_BASE}/kokkos
if [ "${CUDA}" = "ON" ]; then
    KOKKOS_CXX_COMPILER="${KOKKOS_INSTALL_PATH}/bin/nvcc_wrapper"
    KOKKOS_CXX_REPO_COMPILER="${KOKKOS_REPOSITORY_PATH}/bin/nvcc_wrapper"
else
    KOKKOS_CXX_COMPILER="${MY_CXX}"
    KOKKOS_CXX_REPO_COMPILER="${MY_CXX}"
fi
build_install_kokkos

KOKKOSKERNELS_REPOSITORY_PATH=$ROOT/TChem-atm/external/kokkos-kernels
KOKKOSKERNELS_BUILD_PATH=${BUILD_BASE}/kokkos-kernels
KOKKOSKERNELS_INSTALL_PATH=${INSTALL_BASE}/kokkos-kernels
build_install_kokkoskernels

exit
