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
            ELM::read_atm_forcing() - read forcing NetCDF
            ELM::ReadLandData() - read time-invariant land data from NetCDF
        parallel region
            ELM::InitSnowLayers() - initial set of snow layers
            ELM::InitTopoSlope() - set minimum slope
            ELM::InitMicroTopo() - set topography parameters
    Time-stepping
        serial region
            ELM::interp_monthly_veg() - read from phenology NetCDF if needed
        parallel region
            ELM::InitTimestep() - initialize a few variables
            ELM::get_atm_timestep() - process forcing data for this time-step
            ELM::satellite_phenology() - process prescribed vegetation data

            CanopyHydrology calls
            ELM::canopy_hydrology::interception()
            ELM::canopy_hydrology::ground_flux()
            ELM::canopy_hydrology::fraction_wet()
            ELM::canopy_hydrology::snow_init()
            ELM::canopy_hydrology::fraction_h2osfc()

            SurfaceRadiation calls
            ELM::surface_radiation::canopy_sunshade_fractions()
            ELM::surface_radiation::initialize_flux()
            ELM::surface_radiation::total_absorbed_radiation()
            ELM::surface_radiation::layer_absorbed_radiation()
            ELM::surface_radiation::reflected_radiation()

            CanopyTemperature calls
            ELM::canopy_temperature::old_ground_temp()
            ELM::canopy_temperature::ground_temp()
            ELM::canopy_temperature::calc_soilalpha()
            ELM::canopy_temperature::calc_soilbeta()
            ELM::canopy_temperature::humidities()
            ELM::canopy_temperature::ground_properties()
            ELM::canopy_temperature::forcing_height()
            ELM::canopy_temperature::init_energy_fluxes()

            Bareground_fluxes calls
            ELM::bareground_fluxes::initialize_flux()
            ELM::bareground_fluxes::stability_iteration()
            ELM::bareground_fluxes::compute_flux()

            CanopyFluxes calls
            ELM::canopy_fluxes::initialize_flux()
            ELM::canopy_fluxes::stability_iteration()
            ELM::canopy_fluxes::compute_flux()

Documentation
-------------
Documentation for the physics library can be built:
    
    cd docs
    doxygen Doxyfile.in 

Install
-------
sample configure/build/test looks like:

    git clone -b physics https://code.ornl.gov/uec/ELM_Kernels.git ELM_Kernels
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

    
