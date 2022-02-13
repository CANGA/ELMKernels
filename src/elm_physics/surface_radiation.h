/*! \file surface_radiation.h
\brief Functions derived from SurfaceRadiationMod.F90

Compute canopy sun/shade fractions, radiation absorbed by surface/layers, and radiation reflected.

Call sequence:
CanopySunShadeFrac() -> initialize_flux() -> total_absorbed_radiation() -> layer_absorbed_radiation() -> reflected_radiation()
*/
#pragma once

#include "elm_constants.h"
#include "landtype.h"

#include <cassert>
#include <cmath>

namespace ELM::surface_radiation {

/*! Zero out fluxes before surface radiation calculations.

\param[in]  Land                [LandType] struct containing information about landtype
\param[out] sabg_soil           [double] solar radiation absorbed by soil (W/m**2)
\param[out] sabg_snow           [double] solar radiation absorbed by snow (W/m**2)
\param[out] sabg                [double] solar radiation absorbed by ground (W/m**2)
\param[out] sabv                [double] solar radiation absorbed by vegetation (W/m**2)
\param[out] fsa                 [double] solar radiation absorbed (total) (W/m**2)
\param[out] sabg_lyr[nlevsno+1] [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
template <class ArrayD1>
void initialize_flux(const LandType &Land, double &sabg_soil, double &sabg_snow, double &sabg, double &sabv,
                       double &fsa, ArrayD1 sabg_lyr);

/*! Calculate solar flux absorbed by canopy, soil, snow, and ground.

\param[in] Land               [LandType] struct containing information about landtype
\param[in] snl                [int] number of snow layers
\param[in] ftdd[numrad]       [double] down direct flux below canopy per unit direct flux
\param[in] ftid[numrad]       [double] down diffuse flux below canopy per unit direct flux
\param[in] ftii[numrad]       [double] down diffuse flux below canopy per unit diffuse flux
\param[in] forc_solad[numrad] [double] direct beam radiation (W/m**2)
\param[in] forc_solai[numrad] [double] diffuse radiation (W/m**2)
\param[in] fabd[numrad]       [double] flux absorbed by canopy per unit direct flux
\param[in] fabi[numrad]       [double] flux absorbed by canopy per unit diffuse flux
\param[in] albsod[numrad]     [double] soil albedo: direct
\param[in] albsoi[numrad]     [double] soil albedo: diffuse
\param[in] albsnd_hst[numrad] [double] snow albedo, direct , for history files
\param[in] albsni_hst[numrad] [double] snow albedo, diffuse, for history files
\param[in] albgrd[numrad]     [double] ground albedo (direct)
\param[in] albgri[numrad]     [double] ground albedo (diffuse)
\param[out] sabv              [double] solar radiation absorbed by vegetation (W/m**2)
\param[out] fsa               [double] solar radiation absorbed (total) (W/m**2)
\param[out] sabg              [double] solar radiation absorbed by ground (W/m**2)
\param[out] sabg_soil         [double] solar radiation absorbed by soil (W/m**2)
\param[out] sabg_snow         [double] solar radiation absorbed by snow (W/m**2))
\param[out] trd[numrad]       [double] transmitted solar radiation: direct (W/m**2)
\param[out] tri[numrad]       [double] transmitted solar radiation: diffuse (W/m**2)
*/
template <class ArrayD1>
void total_absorbed_radiation(const LandType &Land, const int &snl, const ArrayD1 ftdd, const ArrayD1 ftid,
                     const ArrayD1 ftii, const ArrayD1 forc_solad, const ArrayD1 forc_solai,
                     const ArrayD1 fabd, const ArrayD1 fabi, const ArrayD1 albsod, const ArrayD1 albsoi,
                     const ArrayD1 albsnd_hst, const ArrayD1 albsni_hst, const ArrayD1 albgrd,
                     const ArrayD1 albgri, double &sabv, double &fsa, double &sabg, double &sabg_soil,
                     double &sabg_snow, double trd[numrad], double tri[numrad]);

/*! Compute absorbed flux in each snow layer and top soil layer.

\param[in]  Land                 [LandType] struct containing information about landtype
\param[in]  snl                  [int] number of snow layers
\param[in]  sabg                 [double] solar radiation absorbed by ground (W/m**2)
\param[in]  sabg_snow            [double] solar radiation absorbed by snow (W/m**2)
\param[in]  snow_depth           [double] snow height (m)
\param[in]  flx_absdv[nlevsno+1] [double] direct flux absorption factor: VIS
\param[in]  flx_absdn[nlevsno+1] [double] direct flux absorption factor: NIR
\param[in]  flx_absiv[nlevsno+1] [double] diffuse flux absorption factor: VIS
\param[in]  flx_absin[nlevsno+1] [double] diffuse flux absorption factor: NIR
\param[in]  trd[numrad]          [double] transmitted solar radiation: direct (W/m**2)
\param[in]  tri[numrad]          [double] transmitted solar radiation: diffuse (W/m**2)
\param[out] sabg_lyr[nlevsno+1]  [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
template <class ArrayD1>
void layer_absorbed_radiation(const LandType &Land, const int &snl, const double &sabg, const double &sabg_snow,
                   const double &snow_depth, const ArrayD1 flx_absdv, const ArrayD1 flx_absdn,
                   const ArrayD1 flx_absiv, const ArrayD1 flx_absin, const double trd[numrad],
                   const double tri[numrad], ArrayD1 sabg_lyr);

/*! Calculate reflected solar radiation.

\param[in]  Land               [LandType] struct containing information about landtype
\param[in]  albd[numrad]       [double] surface albedo (direct)
\param[in]  albi[numrad]       [double] surface albedo (diffuse)
\param[in]  forc_solad[numrad] [double] direct beam radiation (W/m**2)
\param[in]  forc_solai[numrad] [double] diffuse radiation (W/m**2)
\param[out] fsr                [double] solar radiation reflected (W/m**2)
*/
template <class ArrayD1>
void reflected_radiation(const LandType &Land, const ArrayD1 albd, const ArrayD1 albi,
                      const ArrayD1 forc_solad, const ArrayD1 forc_solai, double &fsr);

/*!
This subroutine calculates and returns:
1) absorbed PAR for sunlit leaves in canopy layer
2) absorbed PAR for shaded leaves in canopy layer
3) sunlit leaf area
4) shaded  leaf area
5) sunlit leaf area for canopy layer
6) shaded leaf area for canopy layer

\param[in]  Land                [LandType] struct containing information about landtype
\param[in]  nrad                [int] number of canopy layers above snow for radiative transfer
\param[in]  elai                [double] one-sided leaf area index
\param[in]  tlai_z[nlevcan]     [double] tlai increment for canopy layer
\param[in]  fsun_z[nlevcan]     [double] sunlit fraction of canopy layer
\param[in]  forc_solad[numrad]  [double] direct beam radiation (W/m**2)
\param[in]  forc_solai[numrad]  [double] diffuse radiation (W/m**2))
\param[in]  fabd_sun_z[nlevcan] [double] absorbed sunlit leaf direct PAR
\param[in]  fabd_sha_z[nlevcan] [double] absorbed shaded leaf direct PAR
\param[in]  fabi_sun_z[nlevcan] [double] absorbed sunlit leaf diffuse PAR
\param[in]  fabi_sha_z[nlevcan] [double] absorbed shaded leaf diffuse PAR
\param[out] parsun_z[nlevcan]   [double] absorbed PAR for sunlit leaves
\param[out] parsha_z[nlevcan]   [double] absorbed PAR for shaded leaves
\param[out] laisun_z[nlevcan]   [double] sunlit leaf area for canopy layer
\param[out] laisha_z[nlevcan]   [double] shaded leaf area for canopy layer
\param[out] laisun              [double] sunlit leaf area
\param[out] laisha              [double] shaded  leaf area
*/
template <class ArrayD1>
void canopy_sunshade_fractions(const LandType &Land, const int &nrad, const double &elai, const ArrayD1 tlai_z,
                             const ArrayD1 fsun_z, const ArrayD1 forc_solad, const ArrayD1 forc_solai,
                             const ArrayD1 fabd_sun_z, const ArrayD1 fabd_sha_z, const ArrayD1 fabi_sun_z,
                             const ArrayD1 fabi_sha_z, ArrayD1 parsun_z, ArrayD1 parsha_z,
                             ArrayD1 laisun_z, ArrayD1 laisha_z, double &laisun, double &laisha);

} // namespace ELM::surface_radiation

#include "surface_radiation_impl.hh"
