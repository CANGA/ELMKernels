/*! \file surface_albedo.h
\brief Functions derived from SurfaceAlbedoMod.F90

// call sequence -> init_timestep() -> soil_albedo() -> SNICAR_AD_RT() (once for both wavebands) -> ground_albedo() ->
flux_absorption_factor()
// CanopyLayers() -> two_stream_solver()

// will need to finish zenith angle ATS code, Aerosol functions
// need to write coszen function
// need to initialize h2osoi_vol[nlevgrnd]
// need to read in soil color NetCDF, initialize isoicol, albsat, albdry

// need to figure out doalb - during nstep == 0 SurfaceAlbedo doesn't get called and values for the outputs are provided
by SurfaceAlbedoType::InitCold
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
flux_absorption_factor() flx_absd_snw flx_absi_snw mss_cnc_aer_in_fdb - from InitTimestep() to SNICAR_AD_RT()

*/

// not sure if necessary
// maybe for nstep == 0 this should be run, don't call surfalb, but call evrything else?
//  will need to investigate initialization of coupled system - ELM calls dummy nstep == 0 for a reason

#pragma once

#include "elm_constants.h"
#include "land_data.h"
#include "pft_data.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ELM {
namespace surface_albedo {

constexpr double dincmax = 0.25; // maximum lai+sai increment for canopy layer
constexpr double mpe = 1.e-06;   // prevents overflow for division by zero
constexpr double extkn = 0.30;
constexpr double albice[numrad] = {0.8, 0.55};    // albedo land ice by waveband (0=vis, 1=nir)
constexpr double alblak[numrad] = {0.60, 0.40};   // albedo frozen lakes by waveband (0=vis, 1=nir)
constexpr double alblakwi[numrad] = {0.10, 0.10}; // albedo of melting lakes due to puddling, open water, or white ice
                                                  // From D. Mironov (2010) Boreal Env. Research
constexpr double calb = 95.6; // Coefficient for calculating ice "fraction" for lake surface albedo From D. Mironov
                              // (2010) Boreal Env. Research
constexpr bool lakepuddling = false;          // puddling (not extensively tested and currently hardwired off)
constexpr double omegas[numrad] = {0.8, 0.4}; // two-stream parameter omega for snow by band
constexpr double betads = 0.5;                // two-stream parameter betad for snow
constexpr double betais = 0.5;                // two-stream parameter betai for snow

/*
inputs:
urbpoi                                   [bool]   true if urban point, false otherwise
elai                                     [double] one-sided leaf area index with burying by snow
mss_cnc_bcphi[nlevsno]                   [double] mass concentration of hydrophilic BC in snow [kg/kg]
mss_cnc_bcpho[nlevsno]                   [double] mass concentration of hydrophilic BC in snow [kg/kg]
mss_cnc_dst1[nlevsno]                    [double] mass concentration of dust species 1 in snow [kg/kg]
mss_cnc_dst2[nlevsno]                    [double] mass concentration of dust species 2 in snow [kg/kg]
mss_cnc_dst3[nlevsno]                    [double] mass concentration of dust species 3 in snow [kg/kg]
mss_cnc_dst4[nlevsno]                    [double] mass concentration of dust species 4 in snow [kg/kg]

outputs:
vcmaxcintsun                             [double] leaf to canopy scaling coefficient, sunlit leaf vcmax
vcmaxcintsha                             [double] leaf to canopy scaling coefficient, shaded leaf vcmax
albsod[numrad]                           [double] ground albedo (direct)
albsoi[numrad]                           [double] ground albedo (diffuse)
albgrd[numrad]                           [double] direct-beam soil albedo [frc]
albgri[numrad]                           [double] diffuse soil albedo [frc]
albd[numrad]                             [double] surface albedo (direct)
albi[numrad]                             [double] surface albedo (diffuse)
fabd[numrad]                             [double] flux absorbed by canopy per unit direct flux
fabd_sun[numrad]                         [double] flux absorbed by sunlit canopy per unit direct flux
fabd_sha[numrad]                         [double] flux absorbed by shaded canopy per unit direct flux
fabi[numrad]                             [double] flux absorbed by canopy per unit diffuse flux
fabi_sun[numrad]                         [double] flux absorbed by sunlit canopy per unit diffuse flux
fabi_sha[numrad]                         [double] flux absorbed by shaded canopy per unit diffuse flux
ftdd[numrad]                             [double] down direct flux below canopy per unit direct flux
ftid[numrad]                             [double] down diffuse flux below canopy per unit direct flux
ftii[numrad]                             [double] down diffuse flux below canopy per unit diffuse flux
flx_absdv[nlevsno]                       [double] direct flux absorption factor : VIS [frc]
flx_absdn[nlevsno]                       [double] direct flux absorption factor : NIR [frc]
flx_absiv[nlevsno]                       [double] diffuse flux absorption factor : VIS [frc]
flx_absin[nlevsno]                       [double] diffuse flux absorption factor : NIR [frc]
mss_cnc_aer_in_fdb[nlevsno][sno_nbr_aer] [double] mass concentration of all aerosol species for feedback calculation [kg
kg-1]
*/
template <class ArrayD1, class ArrayD2>
void init_timestep(const bool& urbpoi, const double& elai, const ArrayD1 mss_cnc_bcphi, const ArrayD1 mss_cnc_bcpho,
                   const ArrayD1 mss_cnc_dst1, const ArrayD1 mss_cnc_dst2, const ArrayD1 mss_cnc_dst3,
                   const ArrayD1 mss_cnc_dst4, double& vcmaxcintsun, double& vcmaxcintsha, ArrayD1 albsod,
                   ArrayD1 albsoi, ArrayD1 albgrd, ArrayD1 albgri, ArrayD1 albd, ArrayD1 albi, ArrayD1 fabd,
                   ArrayD1 fabd_sun, ArrayD1 fabd_sha, ArrayD1 fabi, ArrayD1 fabi_sun, ArrayD1 fabi_sha, ArrayD1 ftdd,
                   ArrayD1 ftid, ArrayD1 ftii, ArrayD1 flx_absdv, ArrayD1 flx_absdn, ArrayD1 flx_absiv,
                   ArrayD1 flx_absin, ArrayD2 mss_cnc_aer_in_fdb);

/*
Compute ground albedo from weighted snow and soil albedos

inputs:
urbpoi          [bool]     true if urban point, false otherwise
coszen          [double]   solar zenith angle factor
frac_sno        [double]   fraction of ground covered by snow (0 to 1)
albsod[numrad]  [double]   direct-beam soil albedo [frc]
albsoi[numrad]  [double]   diffuse soil albedo [frc]
albsnd[numrad]  [double]   direct-beam snow albedo [frc]
albsni[numrad]  [double]   diffuse snow albedo [frc]

outputs:
albgrd[numrad]  [double] direct-beam ground albedo [frc]
albgri[numrad]  [double] diffuse ground albedo [frc]
*/
template <class ArrayD1>
void ground_albedo(const bool& urbpoi, const double& coszen, const double& frac_sno, const ArrayD1 albsod,
                   const ArrayD1 albsoi, const ArrayD1 albsnd, const ArrayD1 albsni, ArrayD1 albgrd, ArrayD1 albgri);

/*
weight snow layer radiative absorption factors based on snow fraction and soil albedo

inputs:
Land            [LandType] struct containing information about landtype
coszen          [double]   solar zenith angle factor
frac_sno        [double]   fraction of ground covered by snow (0 to 1)
albsod[numrad]  [double]   direct-beam soil albedo [frc]
albsoi[numrad]  [double]   diffuse soil albedo [frc]
albsnd[numrad]  [double]   direct-beam snow albedo [frc]
albsni[numrad]  [double]   diffuse snow albedo [frc]
flx_absd_snw[nlevsno+1][numrad] [double] flux absorption factor for just snow (direct) [frc]
flx_absi_snw[nlevsno+1][numrad] [double] flux absorption factor for just snow (diffuse) [frc]

outputs:
flx_absdv[nlevsno]                       [double] direct flux absorption factor : VIS [frc]
flx_absdn[nlevsno]                       [double] direct flux absorption factor : NIR [frc]
flx_absiv[nlevsno]                       [double] diffuse flux absorption factor : VIS [frc]
flx_absin[nlevsno]                       [double] diffuse flux absorption factor : NIR [frc]
*/
template <class ArrayD1, class ArrayD2>
void flux_absorption_factor(const LandType& Land, const double& coszen, const double& frac_sno, const ArrayD1 albsod,
                            const ArrayD1 albsoi, const ArrayD1 albsnd, const ArrayD1 albsni,
                            const ArrayD2 flx_absd_snw, const ArrayD2 flx_absi_snw, ArrayD1 flx_absdv,
                            ArrayD1 flx_absdn, ArrayD1 flx_absiv, ArrayD1 flx_absin);

/*
Diagnose number of canopy layers for radiative transfer, in increments of dincmax.
Add to number of layers so long as cumulative leaf+stem area does not exceed total
leaf+stem area. Then add any remaining leaf+stem area to next layer and exit the loop.
Do this first for elai and esai (not buried by snow) and then for the part of the
canopy that is buried by snow. Sun/shade big leaf code uses only one layer
(nrad = ncan = 1), triggered by nlevcan == 1.
------------------

tlai_z summed from 1 to nrad = elai
tlai_z summed from 1 to ncan = tlai
tsai_z summed from 1 to nrad = esai
tsai_z summed from 1 to ncan = tsai
------------------


inputs:
urbpoi               [bool]   true if urban point, false otherwise
elai                 [double] one-sided leaf area index with burying by snow
esai                 [double] one-sided stem area index with burying by snow
tlai                 [double] one-sided leaf area index, no burying by snow
tsai                 [double] one-sided stem area index, no burying by snow

outputs:
nrad                 [double] number of canopy layers above snow
ncan                 [double] total number of canopy layers
tlai_z[nlevcan]      [double] leaf area increment for a layer
tsai_z[nlevcan]      [double] stem area increment for a layer
fsun_z[nlevcan]      [double] sunlit fraction of canopy layer
fabd_sun_z[nlevcan]  [double] absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
fabd_sha_z[nlevcan]  [double] absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
fabi_sun_z[nlevcan]  [double] absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
fabi_sha_z[nlevcan]  [double] absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
*/
template <class ArrayD1>
void canopy_layer_lai(const int& urbpoi, const double& elai, const double& esai, const double& tlai, const double& tsai,
                      int& nrad, int& ncan, ArrayD1 tlai_z, ArrayD1 tsai_z, ArrayD1 fsun_z, ArrayD1 fabd_sun_z,
                      ArrayD1 fabd_sha_z, ArrayD1 fabi_sun_z, ArrayD1 fabi_sha_z);

/*
returns true if !urbpoi && coszen > 0 && landtype is vegetated

Land    [LandType] struct containing information about landtype
coszen  [double] solar zenith angle factor
elai    [double] one-sided leaf area index with burying by snow
esai    [double] one-sided stem area index with burying by snow
*/
inline bool vegsol(const LandType& Land, const double& coszen, const double& elai, const double& esai);

/*
returns true if !urbpoi && coszen > 0 && landtype is not vegetated

Land    [LandType] struct containing information about landtype
coszen  [double] solar zenith angle factor
elai    [double] one-sided leaf area index with burying by snow
esai    [double] one-sided stem area index with burying by snow
*/
inline bool novegsol(const LandType& Land, const double& coszen, const double& elai, const double& esai);

/*
DESCRIPTION:
Two-stream fluxes for canopy radiative transfer
Use two-stream approximation of Dickinson (1983) Adv Geophysics
25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
to calculate fluxes absorbed by vegetation, reflected by vegetation,
and transmitted through vegetation for unit incoming direct or diffuse
flux given an underlying surface with known albedo.
Calculate sunlit and shaded fluxes as described by
Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
a multi-layer canopy to calculate APAR profile

inputs:
Land                     [LandType] struct containing information about landtype
nrad                     [int]  number of canopy layers above snow for radiative transfer
coszen                   [double] cos solar zenith angle next time step
t_veg                    [double] vegetation temperature (Kelvin)
fwet                     [double] fraction of canopy that is wet (0 to 1)
elai                     [double] one-sided leaf area index with burying by snow
esai                     [double] one-sided stem area index with burying by snow
tlai_z[nlevcan]          [double] tlai increment for canopy layer
tsai_z[nlevcan]          [double] tsai increment for canopy layer
albgrd[numrad]           [double] ground albedo (direct) (column-level)
albgri[numrad]           [double] ground albedo (diffuse)(column-level)

outputs:
vcmaxcintsun             [double] leaf to canopy scaling coefficient, sunlit leaf vcmax
vcmaxcintsha             [double] leaf to canopy scaling coefficient, shaded leaf vcmax
albd[numrad]             [double] Upward scattered flux above canopy (per unit direct beam flux)
ftid[numrad]             [double] Downward scattered flux below canopy (per unit direct beam flux)
ftdd[numrad]             [double] Transmitted direct beam flux below canopy (per unit direct beam flux)
fabd[numrad]             [double] Flux absorbed by canopy (per unit direct beam flux)
fabd_sun[numrad]         [double] Sunlit portion of fabd
fabd_sha[numrad]         [double] Shaded portion of fabd
albi[numrad]             [double] Upward scattered flux above canopy (per unit diffuse flux)
ftii[numrad]             [double] Downward scattered flux below canopy (per unit diffuse flux)
fabi[numrad]             [double] Flux absorbed by canopy (per unit diffuse flux)
fabi_sun[numrad]         [double] Sunlit portion of fabi
fabi_sha[numrad]         [double] Shaded portion of fabi
fsun_z[nlevcan]          [double] sunlit fraction of canopy layer
fabd_sun_z[nlevcan]      [double] absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
fabd_sha_z[nlevcan]      [double] absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
fabi_sun_z[nlevcan]      [double] absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
fabi_sha_z[nlevcan]      [double] absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
*/
template <class ArrayD1>
void two_stream_solver(const LandType& Land, const int& nrad, const double& coszen, const double& t_veg,
                       const double& fwet, const double& elai, const double& esai, const ArrayD1 tlai_z,
                       const ArrayD1 tsai_z, const ArrayD1 albgrd, const ArrayD1 albgri, const AlbedoVegData& albveg,
                       double& vcmaxcintsun, double& vcmaxcintsha, ArrayD1 albd, ArrayD1 ftid, ArrayD1 ftdd,
                       ArrayD1 fabd, ArrayD1 fabd_sun, ArrayD1 fabd_sha, ArrayD1 albi, ArrayD1 ftii, ArrayD1 fabi,
                       ArrayD1 fabi_sun, ArrayD1 fabi_sha, ArrayD1 fsun_z, ArrayD1 fabd_sun_z, ArrayD1 fabd_sha_z,
                       ArrayD1 fabi_sun_z, ArrayD1 fabi_sha_z);

/*
Soil albedos
Note that soil albedo routine will only compute nonzero soil albedos where coszen > 0

inputs:
Land                       [LandType] struct containing information about landtype
snl                        [int]      number of snow layers
t_grnd                     [double]   ground temperature (Kelvin)
coszen                     [double]   solar zenith angle factor
//lake_icefrac[nlevlak]      [double]   mass fraction of lake layer that is frozen -- removed for now
h2osoi_vol[nlevgrnd]       [double]   volumetric soil water [m3/m3]
albsat[numrad]             [double]   wet soil albedo by color class and waveband (color class designated in
SurfaceAlbedoInitTimeConst) albdry[numrad]             [double]   dry soil albedo by color class and waveband (color
class designated in SurfaceAlbedoInitTimeConst)

outputs:
albsod[numrad]             [double]   direct-beam soil albedo [frc]
albsoi[numrad]             [double]   diffuse soil albedo [frc]
*/
template <class ArrayD1>
void soil_albedo(const LandType& Land, const int& snl, const double& t_grnd, const double& coszen,
                 const ArrayD1 h2osoi_vol, const ArrayD1 albsat, const ArrayD1 albdry, ArrayD1 albsod, ArrayD1 albsoi);

} // namespace surface_albedo
} // namespace ELM

#include "surface_albedo_impl.hh"
