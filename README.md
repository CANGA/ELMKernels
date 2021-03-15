ELM Functional Physics Model
================================
!!work-in-progress!! this project is in active development and is subject to all caveats

Description
-----------
This library contains LSM physics functions derived from E3SM's land model, ELM. 
This library's functional design and cell-by-cell physics implementation is designed 
to expose a general interface for use with task-parallel programming models.

Current Status
--------------
Currently, the conversion of ELM's fixed-phenology vegetation biophysics modules is complete,
as well as a majority of the supporting data processing infrastructure.

A driver for these physics kernels should call this library's functions in the following manner:
    
    Initialization:
        serial region
            ELM::ReadPFTConstants() - read data from params NetCDF file
            ELM::ReadAtmForcing() - read forcing NetCDF
            ELM::ReadLandData() - read time-invariant land data from NetCDF
        parallel region
            ELM::InitSnowLayers() - initial set of snow layers
            ELM::InitTopoSlope() - set minimum slope
            ELM::InitMicroTopo() - set topography parameters
    Time-stepping
        serial region
            ELM::InterpMonthlyVeg() - read from phenology NetCDF if needed
        parallel region
            ELM::InitTimestep() - initialize a few variables
            ELM::GetAtmTimestep() - process forcing data for this time-step
            ELM::SatellitePhenology() - process prescribed vegetation data

            CanopyHydrology calls
            ELM::Interception()
            ELM::GroundFlux()
            ELM::FracWet()
            ELM::SnowInit()
            ELM::FracH2OSfc()

            SurfaceRadiation calls
            ELM::CanopySunShadeFractions()
            ELM::SurfRadZeroFluxes()
            ELM::SurfRadAbsorbed()
            ELM::SurfRadLayers()
            ELM::SurfRadReflected()

            CanopyTemperature calls
            ELM::SaveGroundTemp()
            ELM::CalculateGroundTemp()
            ELM::CalculateSoilAlpha()
            ELM::CalculateSoilBeta()
            ELM::CalculateHumidities()
            ELM::GroundProperties()
            ELM::CalculateForcingHeight()
            ELM::InitializeEnergyFluxes()

            BareGroundFluxes calls
            ELM::InitializeFlux_BG()
            ELM::StabilityIteration_BG()
            ELM::ComputeFlux_BG()

            CanopyFluxes calls
            ELM::InitializeFlux_Can()
            ELM::StabilityIteration_Can()
            ELM::ComputeFlux_Can()

Documentation
-------------
Documentation for the physics library can be built:
    
    cd docs
    doxygen Doxyfile.in 

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

or you can edit buildELMphys.sh with appropriate paths and options and invoke
    
    sh buildELMphys.sh

    
