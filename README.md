ELM Functional Physics Model
================================
This project is in active development.
Beta release V0.1 is nearing completion. 

Description
-----------
This library contains land-surface physics functions derived from E3SM's land model, ELM.

https://github.com/E3SM-Project/E3SM

This library's functional design and cell-by-cell physics implementation is designed 
to expose a general interface for use with hydrology models and data/task-parallel 
programming models.


Installation
--------------
Building this library is simople and should only take a few minutes.

I recommend setting up the folling diretory structure:

export KOKKOS_DIR=/path/to/kokkos/install
export ELM_DIR=/path/to/ELMKernels
export ELM_INPUT_DIR=/path/to/ELM/inout/for/driver     - not included yet

Add these directory definitions to your environement and include them in a place where they
will be sourced by new shells.

Then clone Kokkos
git clone -b master git@github.com:kokkos/kokkos.git $KOKKOS_DIR

And this repository
git clone -b main https://github.com/CANGA/ELMKernels.git $ELM_DIR

Build Kokkos
cd $KOKKOS_DIR
mkdir "build"; cd "build"
cmake ..\
  -DCMAKE_CXX_COMPILER=$path_to_compiler \
  -DCMAKE_INSTALL_PREFIX=${KOKKOS_DIR} \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ARCH_HSW=On \
  -DKokkos_HWLOC_DIR=$path_to_hwloc
make -j6 install

And now we can build this library
cd $ELM_DIR
sh buildELMphys.sh

If you set up both the KOKKOS_DIR AND ELM_INPUT_DIR, that script shpuld build without issue.

Run some tests
./install/bin/test_*

This code has been built with various minor versipns of GCC 11.xx - 13.xx and the past year
of Clang versions. Device simulation is not currently supported, only CPU computation. The
Kokkos driver at driver/kokkos/kokkos_driver.cc demonstrates simple usage of these kernels.

    
