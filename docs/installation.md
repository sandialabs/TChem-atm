### **Obtaning TChem-atm**
<!--We will need to update the github link afther the library is released.-->
TChem-atm is open-source code available on [GitHub](https://github.com/PCLAeroParams/TChem-atm).
It can be downloaded (cloned) with the terminal command
```bash
git clone https://github.com/PCLAeroParams/TChem-atm
```

TChem-atm requires the following third-party libraries:

  * [Tines](https://github.com/sandialabs/Tines.git)
    * [Sacado](https://docs.trilinos.org/dev/packages/sacado/doc/html/index.html)
    * [BLAS/LAPACK](http://www.openblas.net)
    * [YAML](https://github.com/jbeder/yaml-cpp)
    * [Sundials](https://github.com/LLNL/sundials.git)
  * [Gtest](https://github.com/google/googletest.git)
  * [Skywalker](https://github.com/eagles-project/skywalker.git)
  * [Kokkos-kernels](https://github.com/kokkos/kokkos-kernels)

These third-party libraries are submodules in TChem-atm.
Thus, we can initialize them (i.e., download and update) using the following command:

```bash
git submodule update --init --recursive
```

### **Building and installing third-party libraries**

The script `scripts/tpls_bld.sh` clones, builds, and installs TChem-atm's third-party libraries. To use this script, one must provide compiler and configuration information at the top of the script as follows:

```bash
MY_CC=gcc # C++ compiler
MY_CXX=g++ # C++ compiler
MY_FC=gfortran # Fortran compiler
```

To build/install with CUDA `ON` or `OFF`, use the flag:

```bash
CUDA="ON" # Set to ON/OFF to compile TChem-atm with or without NVIDIA-GPUs support.
```

To build/install with HIP `ON` or `OFF`, use the flag:

```bash
HIP="ON" # Set to ON/OFF to compile TChem-atm with or without AMD-GPUs support.
```

Note that both CUDA and HIP cannot be enabled simultaneously.

To enable OpenMP.

```bash
OPENMP="ON"
```

The following path specifies the location of the TChem-atm source code.

```bash
TCHEM_REPOSITORY_PATH=/path/to/tchem-atm/.
```

Determine whether to install OpenBLAS using this script.
```bash
INSTALL_OPENBLAS="ON"
```

This script initializes submodules and manages the installation/building of third-party libraries within the `$PWD/HOST`, `$PWD/CUDA`, or `$PWD/HIP` directories, based on the selected settings: `CUDA/HIP=OFF`, `CUDA=ON`, or `HIP=ON`. Given that Kokkos and Kokkos-kernels support installation/building with both `CUDA` and `HIP` enabled (`ON`) or disabled (`OFF`), the script ensures their deployment in both `$PWD/HOST` and the relevant GPU-specific directories (`$PWD/CUDA` or `$PWD/HIP`). Other libraries, however, are exclusively installed in `$PWD/HOST`. To accommodate all third-party libraries, execute this script initially with `CUDA=OFF`. For systems equipped with NVIDIA or AMD GPUs, a subsequent run with CUDA=`ON` or `HIP=ON`, respectively, is necessary to leverage GPU-specific installations.


## **Building and installing TChem-atm and Tines**

The script `scripts/build_script.sh` builds and installs TChem-atm and Tines. Similarly to the script for third-party libraries, in `scripts/build_script.sh`, one must provide the same compiler information and CUDA/HIP flag.

In addition, this script adds the option to turn `SACADO="ON"` or `SACADO="OFF"` to enable or disable the SACADO library's automatic differentiation capability.

If OpenBLAS is included as part of a module; `USE_THIS_OPENBLAS="MODULE"`.

The installation path for the third-party libraries is specified by `INSTALL_BASE_HOST`, and `ROOT=/path/to/tchem-atm` is where TChem-atm and Tines are installed. If CUDA is `ON`, libraries will be installed in the `ROOT/CUDA` directory. If HIP is `ON`, libraries will be installed in the `ROOT/HIP` directory. Otherwise, they will be installed in the `ROOT/HOST` directory.
 
Finally, To select the build type `BUILD_TYPE`, choose either `DEBUG` or `RELEASE`.

