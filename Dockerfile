FROM ubuntu:22.04

ARG BUILD_TYPE=RELEASE
ARG SACADO=ON
RUN echo "BUILD TYPE:" ${BUILD_TYPE}
RUN echo "SACADO:" ${SACADO}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autoconf \
        clang-format \
        cmake \
        gcc \
        g++ \
        gfortran \
        git \
        lcov \
        make \
        python3 \
        python3-pip \
        libopenblas-dev \
        pkg-config \
        ca-certificates

RUN pip install numpy h5py

COPY . /tchem_dir/

RUN cmake -S /tchem_dir/external/Tines/ext/kokkos -B /build/kokkos_build \
          -DCMAKE_INSTALL_PREFIX="/install/kokkos_install" \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-fopenmp -g" \
          -DKokkos_ENABLE_SERIAL=ON \
          -DKokkos_ENABLE_OPENMP=ON \
          -DKokkos_ENABLE_CUDA=OFF \
          -DKokkos_ENABLE_CUDA_CONSTEXPR=OFF \
          -DKokkos_ENABLE_CUDA_LAMBDA=OFF
WORKDIR /build/kokkos_build/
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/kokkos-kernels -B /build/kokkoskernels_build \
          -DCMAKE_INSTALL_PREFIX="/install/kokkoskernels_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_EXE_LINKER_FLAGS="-lgfortran" \
          -DKokkosKernels_ENABLE_EXAMPLES=OFF \
          -DKokkosKernels_ENABLE_EXPERIMENTAL=OFF \
          -DKokkosKernels_ENABLE_TESTS=OFF \
          -DKokkosKernels_ENABLE_COMPONENT_BLAS=OFF \
          -DKokkosKernels_ENABLE_COMPONENT_GRAPH=OFF \
          -DKokkosKernels_ENABLE_COMPONENT_LAPACK=OFF \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DKokkos_ROOT=/install/kokkos_install
WORKDIR /build/kokkoskernels_build
RUN make -j8 \
    && make install

RUN cmake -S /tchem_dir/external/Tines/ext/gtest -B /build/gtest_build \
          -DCMAKE_INSTALL_PREFIX="/install/gtest_install" \
          -DCMAKE_CXX_COMPILER=g++
WORKDIR /build/gtest_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Tines/ext/yaml -B /build//yaml_build \
          -DCMAKE_INSTALL_PREFIX="/install/yaml_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_CXX_FLAGS="-g -c" \
          -DCMAKE_EXE_LINKER_FLAGS="" \
          -DCMAKE_BUILD_TYPE=RELEASE
WORKDIR /build/yaml_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Skywalker -B /build/skywalker_build \
          -DCMAKE_INSTALL_PREFIX="/install/skywalker_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DSKYWALKER_PRECISION=double \
          -DCMAKE_BUILD_TYPE=RELEASE
WORKDIR /build/skywalker_build

RUN make -j \
    && make install
RUN cmake -S /tchem_dir/external/Sundials -B /build/sundials_build \
          -DCMAKE_INSTALL_PREFIX="/install/sundials_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_C_FLAGS="-g" \
          -DCMAKE_EXE_LINKER_FLAGS="" \
          -DENABLE_CALIPER:BOOL=OFF \
          -DENABLE_KOKKOS=ON \
          -DKokkos_DIR=/install/kokkos_install \
          -DENABLE_KOKKOS_KERNELS=ON \
          -DKokkosKernels_DIR=/install/kokkoskernels_install \
          -DCMAKE_BUILD_TYPE=RELEASE
WORKDIR /build/sundials_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Tines/src -B /build/tines_build \
          -DCMAKE_INSTALL_PREFIX="/install/tines_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_EXE_LINKER_FLAGS="-lgfortran" \
          -DTINES_ENABLE_DEBUG=OFF \
          -DTINES_ENABLE_VERBOSE=OFF \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DTINES_ENABLE_TEST=OFF \
          -DTINES_ENABLE_EXAMPLE=OFF \
          -DSUNDIALS_INSTALL_PATH=/install/sundials_install \
          -DYAML_INSTALL_PATH=/install/yaml_install \
          -DKOKKOS_INSTALL_PATH=/install/kokkos_install \
          -DOPENBLAS_INSTALL_PATH=`/usr/lib64` \
          -DGTEST_INSTALL_PATH=/install/gtest_install
WORKDIR /build/tines_build
RUN make -j \
    && make install


RUN cmake -S /tchem_dir/src -B /tchem_build \
          -DCMAKE_INSTALL_PREFIX=/tchem_install \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_EXE_LINKER_FLAGS="-lgfortran" \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DTCHEM_ATM_ENABLE_VERBOSE=OFF \
          -DTCHEM_ATM_ENABLE_TEST=ON \
          -DTCHEM_ATM_ENABLE_EXAMPLE=ON \
          -DTCHEM_ATM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY=${SACADO} \
          -DTCHEM_ATM_ENABLE_COVERAGE=ON \
          -DKOKKOS_INSTALL_PATH=/install/kokkos_install \
          -DTINES_INSTALL_PATH=/install/tines_install \
          -DSUNDIALS_INSTALL_PATH=/install/sundials_install \
          -DTCHEM_ATM_ENABLE_SKYWALKER=ON \
          -DTCHEM_ATM_ENABLE_REAL_TYPE="double" \
          -DSKYWALKER_INSTALL_PATH=/install/skywalker_install \
          -DTCHEM_ATM_ENABLE_KOKKOSKERNELS=ON \
          -DKOKKOSKERNELS_INSTALL_PATH=/install/kokkoskernels_install \
          -DGTEST_INSTALL_PATH=/install/gtest_install
WORKDIR /tchem_build
RUN make -j \
    && make install

WORKDIR /tchem_install
ENV TCHEM_INSTALL_PATH=/tchem_install/
