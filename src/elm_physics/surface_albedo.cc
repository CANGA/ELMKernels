#include <cmath>
#include <stdexcept>

#include "elm_constants.h"
#include "landtype.h"
#include "surface_albedo.h"


// for build testing
//#include "array.hh"
//using ArrayB1 = ELM::Array<bool, 1>;
//using ArrayI1 = ELM::Array<int, 1>;
//using ArrayS1 = ELM::Array<std::string, 1>;
//using ArrayD1 = ELM::Array<double, 1>;
//using ArrayD2 = ELM::Array<double, 2>;

namespace ELM {
namespace SurfaceAlbedo {

// need to figure out vcmaxcintsha vcmaxcintsun - probably do need to calc twice  - check! 
// call sequence -> SurfAlbInitTimestep() -> SoilAlbedo() -> SNICAR_AD_RT() (once for both wavebands) -> GroundAlbedo() -> SnowAbsorptionFactor()
// CanopyLayers() -> TwoStream()

// will need to finish zenith angle ATS code, Aerosol functions
// need to initialize h2osoi_vol[nlevgrnd]
// need to read in soil color NetCDF, initialize isoicol, albsat, albdry
// need to figure out doalb - during nstep == 0 SurfaceAlbedo doesn't get called and values for the outputs are provided by SurfaceAlbedoType::InitCold
//  why does SurfAlb calc at nstep+1? there must be a reason
//  from CLM 4.5 tech note:
/*

The land model calculations are implemented in two steps. The land
model proceeds with the calculation of surface energy, constituent, momentum, and
radiative fluxes using the snow and soil hydrologic states from the previous time step.
The land model then updates the soil and snow hydrology calculations based on these
fluxes. These fields are passed to the atmosphere (Table 2.4). The albedos sent to the
atmosphere are for the solar zenith angle at the next time step but with surface conditions
from the current time step.

So it's for the ATM coupling
so we can probably calc at beginning of nstep with current solar zenith angle


these don't need to be persistent at the driver level, but will need to be passed from SNICAR_AD_RT() to SnowAbsorptionFactor()
flx_absd_snw
flx_absi_snw
mss_cnc_aer_in_fdb - from SurfAlbInitTimestep() to SNICAR_AD_RT()

*/  

// not sure if necessary
// maybe for nstep == 0 this should be run, don't call surfalb, but call evrything else?
//  will need to investigate initialization of coupled system - ELM calls dummy nstep == 0 for a reason
//void SurfaceAlbedo_InitCold() {
//albgrd[numrad] = {0.2};
//albgri[numrad] = {0.2};
//albsod[numrad] = {0.2};
//albsoi[numrad] = {0.2};
//albd[numrad] = {0.2};
//albi[numrad] = {0.2};
//fabi[numrad] = {0.0};
//fabd[numrad] = {0.0};
//fabi_sun[numrad] = {0.0};
//fabd_sun[numrad] = {0.0};
//fabd_sha[numrad] = {0.0};
//fabi_sha[numrad] = {0.0};
//ftdd[numrad] = {1.0};
//ftid[numrad] = {0.0};
//ftii[numrad] = {1.0};
//}



// FUNCTION to return the cosine of the solar zenith angle.
// Assumes 365.0 days/year
//jday    Julian cal day (1.xx to 365.xx)
//lat     Centered latitude (radians)
//lon     Centered longitude (radians)
//declin  Solar declination (radians)
double calc_cosz(const double jday, const double lat, const double lon, const double declin) {

  return sin(lat) * sin(declin) - cos(lat) * cos(declin) * cos((jday - floor(jday)) * 2.0 * ELM_PI + lon);
}

/*
returns true if !urbpoi && coszen > 0 && landtype is vegetated

inputs:
Land    [LandType] struct containing information about landtype
coszen  [double]   solar zenith angle factor
elai    [double] one-sided leaf area index with burying by snow
esai    [double] one-sided stem area index with burying by snow
*/
//inline bool vegsol(const LandType &Land, const double &coszen, const double &elai, const double &esai) {
//  if (!Land.urbpoi && coszen > 0.0 && (Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0) {
//    return true;
//  } else {
//    return false;
//  }
//}


/*
returns true if !urbpoi && coszen > 0 && landtype is not vegetated

inputs:
inputs:
Land    [LandType] struct containing information about landtype
coszen  [double]   solar zenith angle factor
elai    [double] one-sided leaf area index with burying by snow
esai    [double] one-sided stem area index with burying by snow
*/
//inline bool novegsol(const LandType &Land, const double &coszen, const double &elai, const double &esai) {
//  if (!Land.urbpoi && coszen > 0.0) {
//    if (!((Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0)) {
//      return true;
//    }
//  }
//  return false;
//}



//subroutine SurfaceAlbedoInitTimeConst(bounds)
//    !
//    ! !DESCRIPTION:
//    ! Initialize module time constant variables
//    !
//    ! !USES:
//    use shr_log_mod, only : errMsg => shr_log_errMsg
//    use fileutils  , only : getfil
//    use abortutils , only : endrun
//    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_pio_openfile, ncd_pio_closefile
//    use spmdMod    , only : masterproc
//    !
//    ! !ARGUMENTS:
//    type(bounds_type), intent(in) :: bounds  
//    !
//    ! !LOCAL VARIABLES:
//    integer            :: c,g          ! indices
//    integer            :: mxsoil_color ! maximum number of soil color classes
//    type(file_desc_t)  :: ncid         ! netcdf id
//    character(len=256) :: locfn        ! local filename
//    integer            :: ier          ! error status
//    logical            :: readvar 
//    integer  ,pointer  :: soic2d (:)   ! read in - soil color 
//    !---------------------------------------------------------------------
//
//    ! Allocate module variable for soil color
//
//    allocate(isoicol(bounds%begc:bounds%endc)) 
//
//    ! Determine soil color and number of soil color classes 
//    ! if number of soil color classes is not on input dataset set it to 8
//
//    call getfil (fsurdat, locfn, 0)
//    call ncd_pio_openfile (ncid, locfn, 0)
//
//    call ncd_io(ncid=ncid, varname='mxsoil_color', flag='read', data=mxsoil_color, readvar=readvar)
//    if ( .not. readvar ) mxsoil_color = 8  
//
//    allocate(soic2d(bounds%begg:bounds%endg)) 
//    call ncd_io(ncid=ncid, varname='SOIL_COLOR', flag='read', data=soic2d, dim1name=grlnd, readvar=readvar)
//    if (.not. readvar) then
//       call endrun(msg=' ERROR: SOIL_COLOR NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
//    end if
//    do c = bounds%begc, bounds%endc
//       g = col_pp%gridcell(c)
//       isoicol(c) = soic2d(g)
//    end do
//    deallocate(soic2d)
//
//    call ncd_pio_closefile(ncid)
//
//    ! Determine saturated and dry soil albedos for n color classes and 
//    ! numrad wavebands (1=vis, 2=nir)
//
//    allocate(albsat(mxsoil_color,numrad), albdry(mxsoil_color,numrad), stat=ier)
//    if (ier /= 0) then
//       write(iulog,*)'allocation error for albsat, albdry'
//       call endrun(msg=errMsg(__FILE__, __LINE__)) 
//    end if
//
//    if (masterproc) then
//       write(iulog,*) 'Attempting to read soil colo data .....'
//    end if
//    
//    if (mxsoil_color == 8) then
//       albsat(1:8,1) = (/0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8/)
//       albsat(1:8,2) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
//       albdry(1:8,1) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
//       albdry(1:8,2) = (/0.48_r8,0.44_r8,0.40_r8,0.36_r8,0.32_r8,0.28_r8,0.24_r8,0.20_r8/)
//    else if (mxsoil_color == 20) then
//       albsat(1:20,1) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
//            0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
//       albsat(1:20,2) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
//            0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
//       albdry(1:20,1) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
//            0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
//       albdry(1:20,2) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
//            0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)
//    else
//       write(iulog,*)'maximum color class = ',mxsoil_color,' is not supported'
//       call endrun(msg=errMsg(__FILE__, __LINE__)) 
//    end if
//
//    ! Set alblakwi
//    alblakwi(:) = lake_melt_icealb(:)
//
//  end subroutine SurfaceAlbedoInitTimeConst

} // namespace SurfaceAlbedo
} // namespace ELM
