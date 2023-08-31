#!/usr/bin/env bash

if [ -d "build" ]; then 
  rm -rf "build"
fi
if [ -d "install" ]; then 
  rm -rf "install"
fi
ORIGIN_DIR=`pwd`
mkdir "install" ; mkdir "build" ; cd "build"
cmake ..\
    -DBUILD_SHARED_LIBS:BOOL=true \
    -DCMAKE_CXX_COMPILER:STRING=mpicxx \
    -DCMAKE_C_COMPILER:STRING=mpicc \
    -DSOURCE_PATH:FILEPATH=${ORIGIN_DIR} \
    -DCMAKE_INSTALL_PREFIX:FILEPATH="${ORIGIN_DIR}/install" \
    -DCMAKE_BUILD_TYPE:STRING=Debug \
    -DKokkos_ROOT:FILEPATH=${KOKKOS_DIR} \
    -DENABLE_KOKKOS:BOOL=ON \
    -DENABLE_TESTS:BOOL=ON \
    -DELM_INPUT_DATA_PREFIX:FILEPATH=${ELM_INPUT_DIR}
    ##-DCMAKE_CXX_FLAGS:STRING="-pedantic-errors -Wall -Wextra"
make -j6 VERBOSE=1
make install
cd $ORIGIN_DIR
