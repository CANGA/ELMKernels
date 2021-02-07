#!/bin/bash

if [ -d "build" ]; then 
  rm -rf "build"
fi
if [ -d "install" ]; then 
  rm -rf "install"
fi
ORIGIN_DIR=${pwd}
mkdir "install" ; mkdir "build" ; cd "build"
cmake ..\
    -DBUILD_SHARED_LIBS:BOOL=true \
    -DCMAKE_CXX_COMPILER:STRING=g++ \
    -DCMAKE_INSTALL_PREFIX:FILEPATH=`pwd`/../install \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DKokkos_ROOT:FILEPATH=${KOKKOS_DIR} \
    -DENABLE_Kokkos:BOOL=ON \
    -DENABLE_CC:BOOL=OFF
make VERBOSE=1
make install
cd $ORIGIN_DIR
