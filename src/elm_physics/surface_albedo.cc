
#include "surface_albedo.h"

namespace ELM {
namespace surface_albedo {

// need to figure out vcmaxcintsha vcmaxcintsun - probably do need to calc twice  - check!
// call sequence -> SurfAlbInitTimestep() -> soil_albedo() -> SNICAR_AD_RT() (once for both wavebands) ->
// ground_albedo() -> flux_absorption_factor() CanopyLayers() -> two_stream_solver()

// will need to finish zenith angle ATS code, Aerosol functions
// need to initialize h2osoi_vol[nlevgrnd]
// need to read in soil color NetCDF, initialize isoicol, albsat, albdry
// need to figure out doalb - during nstep == 0 SurfaceAlbedo doesn't get called and values for the outputs are provided
// by SurfaceAlbedoType::InitCold
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


these don't need to be persistent at the driver level, but will need to be passed from SNICAR_AD_RT() to
flux_absorption_factor() flx_absd_snw flx_absi_snw mss_cnc_aer_in_fdb - from SurfAlbInitTimestep() to SNICAR_AD_RT()

*/

// not sure if necessary
// maybe for nstep == 0 this should be run, don't call surfalb, but call evrything else?
//  will need to investigate initialization of coupled system - ELM calls dummy nstep == 0 for a reason
// void SurfaceAlbedo_InitCold() {
// albgrd[numrad] = {0.2};
// albgri[numrad] = {0.2};
// albsod[numrad] = {0.2};
// albsoi[numrad] = {0.2};
// albd[numrad] = {0.2};
// albi[numrad] = {0.2};
// fabi[numrad] = {0.0};
// fabd[numrad] = {0.0};
// fabi_sun[numrad] = {0.0};
// fabd_sun[numrad] = {0.0};
// fabd_sha[numrad] = {0.0};
// fabi_sha[numrad] = {0.0};
// ftdd[numrad] = {1.0};
// ftid[numrad] = {0.0};
// ftii[numrad] = {1.0};
//}

// FUNCTION to return the cosine of the solar zenith angle.
// Assumes 365.0 days/year
// jday    Julian cal day (1.xx to 365.xx)
// lat     Centered latitude (radians)
// lon     Centered longitude (radians)
// declin  Solar declination (radians)
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
// inline bool vegsol(const LandType &Land, const double &coszen, const double &elai, const double &esai) {
//   if (!Land.urbpoi && coszen > 0.0 && (Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0) {
//     return true;
//   } else {
//     return false;
//   }
// }

/*
returns true if !urbpoi && coszen > 0 && landtype is not vegetated

inputs:
inputs:
Land    [LandType] struct containing information about landtype
coszen  [double]   solar zenith angle factor
elai    [double] one-sided leaf area index with burying by snow
esai    [double] one-sided stem area index with burying by snow
*/
// inline bool novegsol(const LandType &Land, const double &coszen, const double &elai, const double &esai) {
//   if (!Land.urbpoi && coszen > 0.0) {
//     if (!((Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0)) {
//       return true;
//     }
//   }
//   return false;
// }

} // namespace surface_albedo
} // namespace ELM
