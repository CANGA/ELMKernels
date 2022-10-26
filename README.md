ELM Functional Physics Model
================================
This project is in active development.

Description
-----------
This library contains land-surface physics functions derived from E3SM's land
model, ELM.

https://github.com/E3SM-Project/E3SM

This library's cell-by-cell physics implementation is designed to expose a
general interface for use with hydrology models anddata/task-parallel 
programming models.



Installation
--------------
Building this library is simple and should only take a few minutes.

I recommend setting up the following directory structure:

```
export KOKKOS_DIR=/path/to/kokkos/install
export ELM_DIR=/path/to/ELMKernels
export ELM_INPUT_DIR=/path/to/ELM/inout/for/driver     ## input data not yet included - will update soon
```

Add these directory definitions to your environment and include them in a place where they
will be sourced by new shells.

Then clone Kokkos
```
git clone -b master "https://github.com/kokkos/kokkos.git" $KOKKOS_DIR
```

And this repository
```
git clone -b main "https://github.com/CANGA/ELMKernels.git" $ELM_DIR
```

Build Kokkos
Your build options may differ from those shown here
```
cd $KOKKOS_DIR
mkdir "build"; cd "build"
cmake ..\
  -DCMAKE_CXX_COMPILER=$path_to_compiler \
  -DCMAKE_INSTALL_PREFIX=${KOKKOS_DIR} \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ARCH_HSW=On \
  -DKokkos_HWLOC_DIR=$path_to_hwloc
make -j6 install
```

And now we can build this library
```
cd $ELM_DIR
sh buildELMphys.sh
```

If you set up both the KOKKOS_DIR and ELM_INPUT_DIR, the build script work without issue.

Run some tests
```
./install/bin/test_*
```

This code has been built with various minor versions of GCC 11.xx - 13.xx and the past year
of Clang versions. Device simulation is not currently supported, only CPU computation. The
Kokkos driver at driver/kokkos/kokkos_driver.cc demonstrates simple usage of these kernels.

    
