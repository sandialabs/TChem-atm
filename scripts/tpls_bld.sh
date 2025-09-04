#!/bin/bash

#=======================================================================================
# example script to download, build, install, & test:
#     openblas, kokkos, gtest, sundials, skywalker
#=======================================================================================
# User configuration  -- begin

# example:


# set CUDA="OFF" and HIP= OFF to use openMP
# JFLAG is for makefile compilation on multiple cores, here 4
# if CUDA="ON", make sure nvcc is in your system PATH

MY_CC=gcc
MY_CXX=g++
MY_FC=gfortran
JFLAG="-j 10"
# Kokkos backends
# set CUDA="ON" to use a NVIDIA GPU
CUDA="OFF"
# set HIP="ON" to use a AMD GPU
HIP="OFF"
#Note that both CUDA and HIP cannot be enabled simultaneously.
OPENMP="ON"

# Note: that we are getting a runtime error if OpenMP is ON while HIP is also ON.
if [ "${HIP}" = "ON" ]; then
  OPENMP="OFF"
fi

# Path to TChem repository.
# Where we git clone TChem-atm. Hint: look for the external folder.
TCHEM_REPOSITORY_PATH=/path/to/tchem/

# Whether to install OpenBLAS using this script.
# Some clusters already have a version of OpenBLAS installed.
INSTALL_OPENBLAS="ON"
# This script will install/build TPLs in the current working directory.
# will build under BUILD_BASE
# will install under INSTALL_BASE
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
  run_this="cd $TCHEM_REPOSITORY_PATH;git submodule update --init --recursive;cd -"
  eval $run_this
}

build_openblas (){
echo "Building OpenBLAS:"
cd ${OPENBLAS_REPOSITORY_PATH}
make CC=${MY_CC} FC=${MY_FC} HOSTCC=${MY_CC} USE_OPENMP=1
}
clean_openblas (){
echo "Cleaning OpenBLAS:"
cd ${OPENBLAS_REPOSITORY_PATH}
make clean
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
# -D CMAKE_CXX_FLAGS="-g" \
# -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
cmake \
    -D CMAKE_INSTALL_PREFIX="${KOKKOSKERNELS_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${KOKKOS_CXX_COMPILER}" \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D KokkosKernels_ENABLE_TPL_CUBLAS=OFF \
    -D KokkosKernels_ENABLE_TPL_CUSOLVER=OFF \
    -D KokkosKernels_ENABLE_TPL_CUSPARSE=OFF \
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
cmake \
    -D CMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER=${KOKKOS_CXX_REPO_COMPILER} \
    -D CMAKE_CXX_STANDARD="17" \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_HIP=${HIP} \
    -D Kokkos_ENABLE_OPENMP=${OPENMP} \
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
echo "Building Sundials:"
mkdir ${SUNDIALS_BUILD_PATH}
mkdir ${SUNDIALS_INSTALL_PATH}
cd ${SUNDIALS_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D SUNDIALS_ENABLE_ERROR_CHECKS=OFF \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D ENABLE_KOKKOS=ON \
    -D Kokkos_DIR=${KOKKOS_INSTALL_PATH} \
    -D ENABLE_KOKKOS_KERNELS=ON \
    -D KokkosKernels_DIR=${KOKKOSKERNELS_INSTALL_PATH} \
    -D CMAKE_BUILD_TYPE=RELEASE \
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

REPO_BASE_EXTERNAL=$TCHEM_REPOSITORY_PATH/external/

if [ "${CUDA}" = "ON" ]; then
    BUILD_BASE=${PWD}/CUDA/build
    INSTALL_BASE=${PWD}/CUDA/install
elif [ "${HIP}" = "ON" ]; then
    BUILD_BASE=${PWD}/HIP/build
    INSTALL_BASE=${PWD}/HIP/install
else
    BUILD_BASE=${PWD}/HOST/build
    INSTALL_BASE=${PWD}/HOST/install
fi

mkdir -p ${BUILD_BASE}
mkdir -p ${INSTALL_BASE}

if [ "${HIP}" = "ON" ]; then
    MY_CC="hipcc"
    MY_CXX="hipcc"
fi

# Note: git submodule update --init --recursive
get_submodules

#only build for host
if [[ "${CUDA}" = "OFF" &&  "${HIP}" = "OFF" ]]; then
  # clone tpls
  if [ "${INSTALL_OPENBLAS}" = "ON" ]; then
    OPENBLAS_REPOSITORY_PATH=${REPO_BASE_EXTERNAL}Tines/ext/OpenBLAS
    OPENBLAS_INSTALL_PATH=${INSTALL_BASE}/openblas
    clean_openblas
    build_openblas
    install_openblas
  fi


  GTEST_REPOSITORY_PATH=${REPO_BASE_EXTERNAL}Tines/ext/gtest
  GTEST_BUILD_PATH=${BUILD_BASE}/gtest
  GTEST_INSTALL_PATH=${INSTALL_BASE}/gtest
  build_install_gtest
  #
  YAML_REPOSITORY_PATH=${REPO_BASE_EXTERNAL}Tines/ext/yaml
  YAML_BUILD_PATH=${BUILD_BASE}/yaml
  YAML_INSTALL_PATH=${INSTALL_BASE}/yaml
  build_install_yaml


  SKYWALKER_REPOSITORY_PATH=$REPO_BASE_EXTERNAL/Skywalker
  SKYWALKER_BUILD_PATH=${BUILD_BASE}/skywalker
  SKYWALKER_INSTALL_PATH=${INSTALL_BASE}/skywalker
  build_install_skywalker
fi

KOKKOS_REPOSITORY_PATH=${REPO_BASE_EXTERNAL}Tines/ext/kokkos
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

KOKKOSKERNELS_REPOSITORY_PATH=$REPO_BASE_EXTERNAL/kokkos-kernels
KOKKOSKERNELS_BUILD_PATH=${BUILD_BASE}/kokkos-kernels
KOKKOSKERNELS_INSTALL_PATH=${INSTALL_BASE}/kokkos-kernels

build_install_kokkoskernels
SUNDIALS_REPOSITORY_PATH=$REPO_BASE_EXTERNAL/Sundials
SUNDIALS_BUILD_PATH=${BUILD_BASE}/sundials
SUNDIALS_INSTALL_PATH=${INSTALL_BASE}/sundials
build_install_sundials

exit
