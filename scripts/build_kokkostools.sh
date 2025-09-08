
MY_CC=gcc
MY_CXX=g++
JFLAG="-j 10"

REPO_BASE=${PWD}/kokkos-tools
BUILD_BASE=${PWD}/kokkos-tools
INSTALL_BASE=${PWD}/kokkos-tools

get_kokkostools (){
echo "get kokkostools:"
if [ -d "${KOKKOSTOOLS_REPOSITORY_PATH}" ] && [ "$(ls -A ${KOKKOSTOOLS_REPOSITORY_PATH})" ]; then
  echo "${KOKKOSTOOLS_REPOSITORY_PATH} exists and is not empty ... aborting clone"; return
fi
git clone https://github.com/kokkos/kokkos-tools.git ${KOKKOSTOOLS_REPOSITORY_PATH}
}

build_install_kokkostools(){
echo "Building kokkos tools:"
mkdir ${KOKKOSTOOLS_BUILD_PATH}
cd ${KOKKOSTOOLS_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX=${KOKKOSTOOLS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    ${KOKKOSTOOLS_REPOSITORY_PATH}
make ${JFLAG} install
}
KOKKOSTOOLS_REPOSITORY_PATH=${REPO_BASE}/main
KOKKOSTOOLS_BUILD_PATH=${BUILD_BASE}/build
KOKKOSTOOLS_INSTALL_PATH=${INSTALL_BASE}/install

get_kokkostools
build_install_kokkostools

exit
