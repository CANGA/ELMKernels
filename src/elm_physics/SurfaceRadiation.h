/*! \file SurfaceRadiation.h
\brief Functions derived from SurfaceRadiationMod.F90

Compute canopy sun/shade fractions, radiation absorbed by surface/layers, and radiation reflected.

Call sequence:
CanopySunShadeFrac() -> SurfRadZeroFluxes() -> SurfRadAbsorbed() -> SurfRadLayers() -> SurfRadReflected()
*/
#pragma once

#include "clm_constants.h"
#include "landtype.h"
#include <assert.h>
#include <cmath>

namespace ELM {

/*! Zero out fluxes before surface radiation calculations.

\param[in]  Land                [LandType] struct containing information about landtype
\param[out] sabg_soil           [double] solar radiation absorbed by soil (W/m**2)
\param[out] sabg_snow           [double] solar radiation absorbed by snow (W/m**2)
\param[out] sabg                [double] solar radiation absorbed by ground (W/m**2)
\param[out] sabv                [double] solar radiation absorbed by vegetation (W/m**2)
\param[out] fsa                 [double] solar radiation absorbed (total) (W/m**2)
\param[out] sabg_lyr[nlevsno+1] [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
template <class dArray_type>
void SurfRadZeroFluxes(const LandType &Land, double &sabg_soil, double &sabg_snow, double &sabg, double &sabv,
                       double &fsa, dArray_type sabg_lyr);

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
template <class dArray_type>
void SurfRadAbsorbed(const LandType &Land, const int &snl, const dArray_type ftdd, const dArray_type ftid,
                     const dArray_type ftii, const dArray_type forc_solad, const dArray_type forc_solai,
                     const dArray_type fabd, const dArray_type fabi, const dArray_type albsod, const dArray_type albsoi,
                     const dArray_type albsnd_hst, const dArray_type albsni_hst, const dArray_type albgrd,
                     const dArray_type albgri, double &sabv, double &fsa, double &sabg, double &sabg_soil,
                     double &sabg_snow, double *trd, double *tri);

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
template <class dArray_type>
void SurfRadLayers(const LandType &Land, const int &snl, const double &sabg, const double &sabg_snow,
                   const double &snow_depth, const dArray_type flx_absdv, const dArray_type flx_absdn,
                   const dArray_type flx_absiv, const dArray_type flx_absin, const double *trd, const double *tri,
                   dArray_type sabg_lyr);

/*! Calculate reflected solar radiation.

\param[in]  Land               [LandType] struct containing information about landtype
\param[in]  albd[numrad]       [double] surface albedo (direct)
\param[in]  albi[numrad]       [double] surface albedo (diffuse)
\param[in]  forc_solad[numrad] [double] direct beam radiation (W/m**2)
\param[in]  forc_solai[numrad] [double] diffuse radiation (W/m**2)
\param[out] fsr                [double] solar radiation reflected (W/m**2)
*/
template <class dArray_type>
void SurfRadReflected(const LandType &Land, const dArray_type albd, const dArray_type albi,
                      const dArray_type forc_solad, const dArray_type forc_solai, double &fsr);

/*!
This subroutine calculates and returns:
1) absorbed PAR for sunlit leaves in canopy layer
2) absorbed PAR for shaded leaves in canopy layer
3) sunlit leaf area
4) shaded  leaf area
5) sunlit leaf area for canopy layer
6) shaded leaf area for canopy layer
7) sunlit fraction of canopy

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
template <class dArray_type>
void CanopySunShadeFractions(const LandType &Land, const int &nrad, const double &elai, const dArray_type tlai_z,
                             const dArray_type fsun_z, const dArray_type forc_solad, const dArray_type forc_solai,
                             const dArray_type fabd_sun_z, const dArray_type fabd_sha_z, const dArray_type fabi_sun_z,
                             const dArray_type fabi_sha_z, dArray_type parsun_z, dArray_type parsha_z,
                             dArray_type laisun_z, dArray_type laisha_z, double &laisun, double &laisha);

} // namespace ELM

#include "SurfaceRadiation_impl.hh"
