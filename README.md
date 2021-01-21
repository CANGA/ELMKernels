ELM Functional Physics Model
================================



Install
-------

sample configure/build/test looks like:

    git clone -b ecoon/physics https://code.ornl.gov/uec/ELM_Kernels.git ELM_Kernels
    mkdir build-dir
    mkdir install-dir
    
    cd build-dir
    cmake \
        -DNetCDF_ROOT:FILEPATH=PATH/TO/NETCDF \
        -DBUILD_SHARED_LIBS:BOOL=true \
        -DCMAKE_INSTALL_PREFIX:FILEPATH=`pwd`/../install-dir
        -DCMAKE_BUILD_TYPE:STRING=Debug \
        -DENABLE_KOKKOS:BOOL=true \
      `pwd`/../ELM_Kernels
      
    make
    make test
    make install

    
