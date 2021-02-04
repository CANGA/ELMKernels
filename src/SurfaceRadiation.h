/* functions derived from SurfaceRadiationMod.F90 

Call sequence:
CanopySunShadeFrac()
SurfRadZeroFluxes()
SurfRadAbsorbed()
SurfRadLayers()
SurfRadReflected
*/
#pragma once 

#include <cmath>
#include <assert.h>
#include "clm_constants.h"
#include "landtype.h"

namespace ELM {

/* SurfaceRadiation::SurfRadZeroFluxes()
DESCRIPTION: zero out fluxes before surface radiation calculations

INPUTS:
Land                [LandType] struct containing information about landtype 

OUTPUTS:
sabg_soil           [double] solar radiation absorbed by soil (W/m**2) 
sabg_snow           [double] solar radiation absorbed by snow (W/m**2)
sabg                [double] solar radiation absorbed by ground (W/m**2)
sabv                [double] solar radiation absorbed by vegetation (W/m**2)
fsa                 [double] solar radiation absorbed (total) (W/m**2)
sabg_lyr[nlevsno+1] [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
template<class dArray_type>
void SurfRadZeroFluxes (
  const LandType& Land,
  double& sabg_soil,
  double& sabg_snow,
  double& sabg,
  double& sabv,
  double& fsa,
  dArray_type sabg_lyr);

/* SurfaceRadiation::SurfRadAbsorbed()
DESCRIPTION: calculate solar flux absorbed by canopy, soil, snow, and ground

INPUTS:
Land               [LandType] struct containing information about landtype 
snl                [int] number of snow layers
ftdd[numrad]       [double] down direct flux below canopy per unit direct flux
ftid[numrad]       [double] down diffuse flux below canopy per unit direct flux
ftii[numrad]       [double] down diffuse flux below canopy per unit diffuse flux
forc_solad[numrad] [double] direct beam radiation (W/m**2)
forc_solai[numrad] [double] diffuse radiation (W/m**2)
fabd[numrad]       [double] flux absorbed by canopy per unit direct flux
fabi[numrad]       [double] flux absorbed by canopy per unit diffuse flux
albsod[numrad]     [double] soil albedo: direct 
albsoi[numrad]     [double] soil albedo: diffuse
albsnd_hst[numrad] [double] snow albedo, direct , for history files
albsni_hst[numrad] [double] snow albedo, diffuse, for history files
albgrd[numrad]     [double] ground albedo (direct) 
albgri[numrad]     [double] ground albedo (diffuse)

OUTPUTS:
sabv               [double] solar radiation absorbed by vegetation (W/m**2)
fsa                [double] solar radiation absorbed (total) (W/m**2)
sabg               [double] solar radiation absorbed by ground (W/m**2)
sabg_soil          [double] solar radiation absorbed by soil (W/m**2)
sabg_snow          [double] solar radiation absorbed by snow (W/m**2))
trd[numrad]        [double] transmitted solar radiation: direct (W/m**2)
tri[numrad]        [double] transmitted solar radiation: diffuse (W/m**2)

*/
template<class dArray_type>
void SurfRadAbsorbed(
  const LandType& Land,
  const int& snl,
  const dArray_type ftdd,
  const dArray_type ftid,
  const dArray_type ftii,
  const dArray_type forc_solad,
  const dArray_type forc_solai,
  const dArray_type fabd,
  const dArray_type fabi,
  const dArray_type albsod,
  const dArray_type albsoi,
  const dArray_type albsnd_hst,
  const dArray_type albsni_hst,
  const dArray_type albgrd,
  const dArray_type albgri,
  double& sabv,
  double& fsa,
  double& sabg,
  double& sabg_soil,
  double& sabg_snow,
  double *trd,
  double *tri);


/* SurfaceRadiation::SurfRadLayers()
DESCRIPTION: compute absorbed flux in each snow layer and top soil layer

INPUTS:
Land                 [LandType] struct containing information about landtype
snl                  [int] number of snow layers
sabg                 [double] solar radiation absorbed by ground (W/m**2)
sabg_snow            [double] solar radiation absorbed by snow (W/m**2)
snow_depth           [double] snow height (m) 
flx_absdv[nlevsno+1] [double] direct flux absorption factor: VIS 
flx_absdn[nlevsno+1] [double] direct flux absorption factor: NIR 
flx_absiv[nlevsno+1] [double] diffuse flux absorption factor: VIS
flx_absin[nlevsno+1] [double] diffuse flux absorption factor: NIR
trd[numrad]          [double] transmitted solar radiation: direct (W/m**2)
tri[numrad]          [double] transmitted solar radiation: diffuse (W/m**2)

OUTPUTS:
sabg_lyr[nlevsno+1]  [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
template<class dArray_type>
void SurfRadLayers(
  const LandType& Land,
  const int& snl,
  const double& sabg,
  const double& sabg_snow,
  const double& snow_depth,
  const dArray_type flx_absdv,
  const dArray_type flx_absdn,
  const dArray_type flx_absiv,
  const dArray_type flx_absin,
  const double *trd,
  const double *tri,
  dArray_type sabg_lyr);


/* SurfaceRadiation::SurfRadReflected()
DESCRIPTION: calculate reflected solar radiation

INPUTS:
Land               [LandType] struct containing information about landtype 
albd[numrad]       [double] surface albedo (direct)
albi[numrad]       [double] surface albedo (diffuse)
forc_solad[numrad] [double] direct beam radiation (W/m**2)
forc_solai[numrad] [double] diffuse radiation (W/m**2)

OUTPUTS:
fsr                [double] solar radiation reflected (W/m**2)
*/
template<class dArray_type>
void SurfRadReflected(
  const LandType& Land,
  const dArray_type albd,
  const dArray_type albi,
  const dArray_type forc_solad,
  const dArray_type forc_solai,
  double& fsr);


/* SurfaceRadiation::CanopySunShadeFractions()
DESCRIPTION:
This subroutine calculates and returns:
1) absorbed PAR for sunlit leaves in canopy layer
2) absorbed PAR for shaded leaves in canopy layer
3) sunlit leaf area
4) shaded  leaf area
5) sunlit leaf area for canopy layer
6) shaded leaf area for canopy layer
7) sunlit fraction of canopy

INPUTS:
Land                [LandType] struct containing information about landtype 
nrad                [int] number of canopy layers above snow for radiative transfer
elai                [double] one-sided leaf area index
tlai_z[nlevcan]     [double] tlai increment for canopy layer
fsun_z[nlevcan]     [double] sunlit fraction of canopy layer
forc_solad[numrad]  [double] direct beam radiation (W/m**2)
forc_solai[numrad]  [double] diffuse radiation (W/m**2))
fabd_sun_z[nlevcan] [double] absorbed sunlit leaf direct PAR
fabd_sha_z[nlevcan] [double] absorbed shaded leaf direct PAR
fabi_sun_z[nlevcan] [double] absorbed sunlit leaf diffuse PAR
fabi_sha_z[nlevcan] [double] absorbed shaded leaf diffuse PAR

OUTPUTS:
parsun_z[nlevcan]   [double] absorbed PAR for sunlit leaves
parsha_z[nlevcan]   [double] absorbed PAR for shaded leaves
laisun_z[nlevcan]   [double] sunlit leaf area for canopy layer
laisha_z[nlevcan]   [double] shaded leaf area for canopy layer
laisun              [double] sunlit leaf area
laisha              [double] shaded  leaf area
*/
template<class dArray_type>
void CanopySunShadeFractions(
  const LandType& Land,
  const int& nrad,
  const double& elai,
  const dArray_type tlai_z,
  const dArray_type fsun_z,
  const dArray_type forc_solad,
  const dArray_type forc_solai,
  const dArray_type fabd_sun_z,
  const dArray_type fabd_sha_z,
  const dArray_type fabi_sun_z,
  const dArray_type fabi_sha_z,
  dArray_type parsun_z,
  dArray_type parsha_z,
  dArray_type laisun_z,
  dArray_type laisha_z,
  double& laisun,
  double& laisha);

}