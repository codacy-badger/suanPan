***PLEASE NOTE THIS PROJECT IS UNDER DEVELOPMENT SO CURRENTLY YOU MAY NOT BE ABLE TO FIND SUFFICIENT INFORMATION***

***PLEASE CHECK EXAMPLE FILES UNDER THE EXAMPLE FOLDER***

# <img src="Resource/suanPan.svg" width="150" align="middle"/> suanPan

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2556012.svg)](https://doi.org/10.5281/zenodo.2556012)
[![License: GPL v3](https://img.shields.io/github/license/TLCFEM/suanPan.svg?color=44cc11)](https://www.gnu.org/licenses/gpl-3.0)
[![release](https://img.shields.io/github/release-pre/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/releases)
[![download](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11)](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11)
[![AppVeyor branch](https://img.shields.io/appveyor/ci/TLCFEM/suanPan/dev.svg?label=master&logo=appveyor)](https://ci.appveyor.com/project/TLCFEM/suanpan/branch/master)
[![AppVeyor branch](https://img.shields.io/appveyor/ci/TLCFEM/suanPan/dev.svg?label=dev&logo=appveyor)](https://ci.appveyor.com/project/TLCFEM/suanpan/branch/dev)
[![Travis (.org) branch](https://img.shields.io/travis/TLCFEM/suanPan/dev.svg?label=master&logo=travis)](https://travis-ci.org/TLCFEM/suanPan)
[![Travis (.org) branch](https://img.shields.io/travis/TLCFEM/suanPan/dev.svg?label=dev&logo=travis)](https://travis-ci.org/TLCFEM/suanPan)
[![codecov](https://codecov.io/gh/TLCFEM/suanPan/branch/dev/graph/badge.svg)](https://codecov.io/gh/TLCFEM/suanPan)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7cb47e58d7dc4c1680c2205c4ba02e72)](https://www.codacy.com/app/TLCFEM/suanPan?utm_source=github.com&utm_medium=referral&utm_content=TLCFEM/suanPan&utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/tlcfem/suanpan/badge)](https://www.codefactor.io/repository/github/tlcfem/suanpan)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/TLCFEM/suanPan.svg?logo=lgtm)](https://lgtm.com/projects/g/TLCFEM/suanPan/context:cpp)
[![language](https://img.shields.io/github/languages/count/TLCFEM/suanPan.svg?color=44cc11)]()
[![GitHub top language](https://img.shields.io/github/languages/top/TLCFEM/suanPan.svg?color=44cc11&logo=c%2B%2B)]()
[![code-size](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)
[![GitHub issues](https://img.shields.io/github/issues/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/issues)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/TLCFEM/suanPan.svg?logo=lgtm)](https://lgtm.com/projects/g/TLCFEM/suanPan/alerts/)


Intro
-----

[**suanPan**](https://suanpan.tk) is a finite element method (FEM) simulation platform for applications in solid mechanics, civil/structural/seismic engineering, etc. The name **suanPan** (in some places such as suffix also it is abbreviated as **suPan**) comes from the term *Suan Pan* (算盤), which is [Chinese abacus](https://en.wikipedia.org/wiki/Suanpan). **suanPan** is written in high quality C++ code and is targeted to provide an efficient, concise and reliable FEM simulation platform.

**suanPan** is partially influenced by popular (non-)commercial FEA packages, such as [ABAQUS UNIFIED FEA](https://www.3ds.com/products-services/simulia/products/abaqus/), [ANSYS](http://www.ansys.com/) and [OpenSees](http://opensees.berkeley.edu/).

[Documentation](https://docs.suanpan.tk) is available.

Features
--------

The highlights of **suanPan** are

-   **suanPan** is fast, both memory and thread safe.
-   **suanPan** is designed based on the [shared memory](https://en.wikipedia.org/wiki/Shared_memory) model and supports parallelism on heterogeneous architectures, for example multi-threaded CPU + optional GPU.
-   **suanPan** is open source and easy to be expanded to incorporate user-defined elements, materials, etc.
-   **suanPan** separates the FEA model part from the linear algebra operation part, which significantly reduces the complexity of development.
-   **suanPan** utilizes the new language features shipped with the latest standards (C++11, C++14 and later ones), such as new STL containers, smart pointers and many other features.
-   **suanPan** supports simple visualization via VTK.

How to Compile
--------------

### Overview

As **suanPan** uses new language features, please use compilers that support ***C++14*** standard. For example,

-   **MSVC 14** (Visual Studio 2015) and/or later version,
-   **GNU GCC 5.4** and/or later version.
-   **Clang 3.4** and/or later version.

On Ubuntu, the external libraries are compiled with **GCC 4.8** using **binutils 2.24**, which is the default configuration on **Ubuntu 14.04 LTS**. On Windows, all libraries are compiled as dynamic libraries with **GCC 5.4** in **MinGW-w64**. To avoid any potential linking error, please use later versions which are available from the **MinGW-w64** project.

1. [**GCC 5.4 x86_64-posix-seh**](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/5.4.0/threads-posix/seh/x86_64-5.4.0-release-posix-seh-rt_v5-rev0.7z)
4. [**GCC 8.1 x86_64-posix-seh**](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/8.1.0/threads-posix/seh/x86_64-8.1.0-release-posix-seh-rt_v6-rev0.7z)

The software is *deliberately designed to disable the backward compatibility*. **suanPan** uses [**CMake**](https://cmake.org/) to manage builds.

### Windows

The package is tested under Windows with **MSVC++**, **GNU GCC**, **Intel C++ Compiler** and **Clang**. The libraries (only 64-bit, no 32-bit support anymore) currently used:

-   **ARPACK** version 0.96 (optional, need a valid Fortran compiler)
-   **SuperLU** version 5.2.1
-   **OpenBLAS** version 0.3.5
-   **TBB** Threading Building Blocks version 2018U3
-   **HDF5** version 1.10.1
-   **MUMPS** version 5.12

are bundled with the source code package.

#### Visual Studio

The solution file under `MSVC` can be compiled directly. The `Release` and `Debug` configurations use external **OpenBLAS**. The `Release-MKL` uses **Intel MKL** and **VTK**. The `Release-GPU` further uses **CUDA** and **magma** to offload sparse matrix solving on GPUs. Simply hit `Build (F7)`  to build the solution. You can change the linked library to other equivalent libraries such as **Intel MKL** manually, if those libraries are available. The bundled **OpenBLAS** is compiled with **GCC 5.4** and dynamic architecture. If it does not work on your machine, please compile your own version. To do so, you may need to download [**MSYS**](http://www.mingw.org/wiki/msys), or [**Cygwin**](https://www.cygwin.com/) with necessary packages such as **gcc** and **perl**.

The compiled program cannot run directly as it depends on other dynamic libraries. Please copy all dynamic link library into the path that can be found by the program.

```text
/Libs/gcc-win/*.dll
```

In addition, the **TBB** libraries shall be copied as well since by default multithreading is enabled.

```text
/Libs/vs/tbb.dll
/Libs/vs/tbbmalloc.dll
/Libs/vs/tbbmalloc_proxy.dll
```

The VS solution can also be generated via **CMake**.

#### GNU GCC

Use **CMake** to generate Makefiles, assume current folder is the root of the package and your are using **MinGW**, the commands should look like this.

``` bash
# current folder is /suanPan-source-code-path
mkdir cmake-build && cd cmake-build
cmake -G "MinGW Makefiles" ..
make
```

For **MSYS**, you may change the generator to "MSYS Makefiles". The GUI may be a better tool for beginners to configure the build. There are several options provided to build the software with different configurations.

After successful compilation, the executable file is under `/cmake-build` folder. The dynamic libraries should also be copied.

```bash
# current folder is cmake-build
cp ../Libs/gcc-win/*.dll .
# run the program
./suanPan.exe
```

### Ubuntu

Again, the shipped **OpenBLAS** may not work on your platform, please compile your own library if any error occurs.

Make sure the compilers are installed.

```bash
sudo apt install gcc g++ gfortran binutils cmake cmake-qt-gui
```

Please do check the versions of those tools. A default configuration is enough for most cases if the platform is not too old. Simply create a build folder next to the source code folder and configure/make the program, such as

``` bash
# current folder is /suanPan-source-code-path
mkdir cmake-build
cd cmake-build
cmake ../suanPan
make
```

The multi-threaded version uses **TBB**, the corresponding path can be added.

```bash
# current folder cmake-build
LD_LIBRAY_PATH=$LD_LIBRARY_PATH:../Libs/gcc-linux
export LD_LIBRAY_PATH
# run the program
./suanPan
```

Dependency
----------

Additional libraries that are used in **suanPan** are

- [Armadillo](http://arma.sourceforge.net/) --- Linear Algebra Library
- [MKL](https://software.intel.com/en-us/mkl) --- High Performance Linear Algebra Driver
- [VTK](https://www.vtk.org/) --- Visualization Toolkit
- [MUMPS](http://mumps.enseeiht.fr/) --- Parallel Sparse Direct Solver

Additional tools may be used by **suanPan**, they are

- [UPX](https://upx.github.io/) --- Executable File Packer

Those libraries may depend on the libraries that are used in the project, such as

- [zlib](https://zlib.net/)
- [Szip](https://support.hdfgroup.org/doc_resource/SZIP/)

