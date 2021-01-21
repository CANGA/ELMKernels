/* functions derived from SurfaceRadiationMod.F90 

Call sequence:
CanopySunShadeFrac()
SurfRadZeroFluxes()
SurfRadAbsorbed()
SurfRadLayers()
SurfRadReflected
*/

#include <cmath>
#include <assert.h>
#include "clm_constants.h"
#include "landtype.h"

namespace ELM {

class SurfaceRadiation {

private:
  double trd[numrad]; // transmitted solar radiation: direct (W/m**2)
  double tri[numrad]; // transmitted solar radiation: diffuse (W/m**2)
  double parveg; // absorbed par by vegetation (W/m**2)

public:

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
fsa_r               [double] rural solar radiation absorbed (total) (W/m**2)
fsun                [double] sunlit fraction of canopy 
sabg_lyr[nlevsno+1] [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
  void SurfRadZeroFluxes (
    const LandType& Land,
    double& sabg_soil,
    double& sabg_snow,
    double& sabg,
    double& sabv,
    double& fsa,
    double& fsa_r,
    double& fsun,
    double* sabg_lyr) {

    // Initialize fluxes
    if (!Land.urbpoi) {
      sabg_soil = 0.0;
      sabg_snow = 0.0;
      sabg      = 0.0;
      sabv      = 0.0;
      fsa       = 0.0;
      if (Land.ltype == istsoil || Land.ltype == istcrop) { fsa_r = 0.0; }
      for (int j = 0; j < nlevsno + 1; j++) { sabg_lyr[j] = 0.0; }
    }
  
    // zero-out fsun for the urban patches
    // the non-urban patches were set prior to this call
    // and split into ed and non-ed specific functions
    if (Land.urbpoi) { fsun = 0.0; }
  }

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
fsa_r              [double] rural solar radiation absorbed (total) (W/m**2)
sabv               [double] solar radiation absorbed by vegetation (W/m**2)
fsa                [double] solar radiation absorbed (total) (W/m**2)
sabg               [double] solar radiation absorbed by ground (W/m**2)
sabg_soil          [double] solar radiation absorbed by soil (W/m**2)
sabg_snow          [double] solar radiation absorbed by snow (W/m**2))

*/
  void SurfRadAbsorbed(
    const LandType& Land,
    const int& snl,
    const double* ftdd,
    const double* ftid,
    const double* ftii,
    const double* forc_solad,
    const double* forc_solai,
    const double* fabd,
    const double* fabi,
    const double* albsod,
    const double* albsoi,
    const double* albsnd_hst,
    const double* albsni_hst,
    const double* albgrd,
    const double* albgri,
  
    double& fsa_r,
    double& sabv,
    double& fsa,
    double& sabg,
    double& sabg_soil,
    double& sabg_snow) {

    double absrad, cad[numrad], cai[numrad];
    if (!Land.urbpoi) {
      for (int ib = 0; ib < numrad; ib++) {
  
        cad[ib] = forc_solad[ib] * fabd[ib];
        cai[ib] = forc_solai[ib] * fabi[ib];
        sabv += cad[ib] + cai[ib];
        fsa  += cad[ib] + cai[ib];
  
        if (ib == 0) { parveg = cad[ib] + cai[ib]; }
        if (Land.ltype == istsoil || Land.ltype == istcrop) { fsa_r  += cad[ib] + cai[ib]; }
  
        // Transmitted = solar fluxes incident on ground
        trd[ib] = forc_solad[ib] * ftdd[ib];
        tri[ib] = forc_solad[ib] * ftid[ib] + forc_solai[ib]*ftii[ib];
        // Solar radiation absorbed by ground surface
        // calculate absorbed solar by soil/snow separately
        absrad  = trd[ib] * (1.0 - albsod[ib]) + tri[ib] * (1.0 - albsoi[ib]);
        sabg_soil += absrad;
        absrad  = trd[ib] * (1.0 - albsnd_hst[ib]) + tri[ib] * (1.0 - albsni_hst[ib]);
        sabg_snow += absrad;
        absrad  = trd[ib]*(1.0 - albgrd[ib]) + tri[ib] * (1.0 - albgri[ib]);
        sabg += absrad;
        fsa  += absrad;
  
        if (Land.ltype == istsoil || Land.ltype == istcrop) { fsa_r += absrad; }
        if (snl == 0) {
          sabg_snow = sabg;
          sabg_soil = sabg;
        }
  
        if (subgridflag == 0) {
          sabg_snow = sabg;
          sabg_soil = sabg;
        }
      } // end of numrad
    } // end if not urban
  }

/* SurfaceRadiation::SurfRadLayers()
DESCRIPTION: compute absorbed flux in each snow layer and top soil layer

INPUTS:
Land                 [LandType] struct containing information about landtype
snl                  [int] number of snow layers
subgridflag          [int] use subgrid fluxes
sabg                 [double] solar radiation absorbed by ground (W/m**2)
sabg_snow            [double] solar radiation absorbed by snow (W/m**2)
snow_depth           [double] snow height (m) 
flx_absdv[nlevsno+1] [double] direct flux absorption factor: VIS 
flx_absdn[nlevsno+1] [double] direct flux absorption factor: NIR 
flx_absiv[nlevsno+1] [double] diffuse flux absorption factor: VIS
flx_absin[nlevsno+1] [double] diffuse flux absorption factor: NIR

OUTPUTS:
sub_surf_abs_SW      [double] percent of solar radiation absorbed below first snow layer (W/M**2)
sabg_pen             [double] solar (rural) radiation penetrating top soisno layer (W/m**2)
sabg_lyr[nlevsno+1]  [double] absorbed radiative flux (pft,lyr) [W/m2]
*/
  void SurfRadLayers(
    const LandType& Land,
    const int& snl,
    const int& subgridflag,
    const double& sabg,
    const double& sabg_snow,
    const double& snow_depth,
    const double* flx_absdv,
    const double* flx_absdn,
    const double* flx_absiv,
    const double* flx_absin,
    double& sub_surf_abs_SW,
    double& sabg_pen,
    double* sabg_lyr) {

    double err_sum = 0.0;
    double sabg_snl_sum;
    // compute absorbed flux in each snow layer and top soil layer,
    // based on flux factors computed in the radiative transfer portion of SNICAR.
    if (!Land.urbpoi) {
  
      sabg_snl_sum = 0.0;
      sub_surf_abs_SW = 0.0;
  
      // CASE1: No snow layers: all energy is absorbed in top soil layer
      if (snl == 0) {
        for (int i = 0; i <= nlevsno; i++) { sabg_lyr[i] = 0.0; }
        sabg_lyr[nlevsno] = sabg; 
        sabg_snl_sum  = sabg_lyr[nlevsno];
      } else { // CASE 2: Snow layers present: absorbed radiation is scaled according to flux factors computed by SNICAR
  
        for (int i = 0; i < nlevsno + 1; i++) {
          sabg_lyr[i] = flx_absdv[i] * trd[0] + flx_absdn[i] * trd[1] + flx_absiv[i] * tri[0] + flx_absin[i] * tri[1];
          // summed radiation in active snow layers:
          // if snow layer is at or below snow surface 
          if (i >= nlevsno - snl) { sabg_snl_sum += sabg_lyr[i]; }
          // if snow layer is below surface snow layer accumulate subsurface flux as a diagnostic
          if (i > nlevsno - snl) { sub_surf_abs_SW += + sabg_lyr[i]; }
        }
  
        // Divide absorbed by total, to get % absorbed in subsurface
        if (sabg_snl_sum != 0.0) { // inequality with double?
          sub_surf_abs_SW = sub_surf_abs_SW / sabg_snl_sum;
        } else {
          sub_surf_abs_SW = 0.0;
        }
        // Error handling: The situation below can occur when solar radiation is 
        // NOT computed every timestep.
        // When the number of snow layers has changed in between computations of the 
        // absorbed solar energy in each layer, we must redistribute the absorbed energy
        // to avoid physically unrealistic conditions. The assumptions made below are 
        // somewhat arbitrary, but this situation does not arise very frequently. 
        // This error handling is implemented to accomodate any value of the
        // radiation frequency.
        // change condition to match sabg_snow isntead of sabg
  
        if (std::abs(sabg_snl_sum - sabg_snow) > 0.00001) {
          if (snl == 0) {
            for (int j = 0; j < nlevsno; j++) { sabg_lyr[j] = 0.0; }
            sabg_lyr[nlevsno] = sabg;
          } else if (snl == 1) {
            for (int j = 0; j < nlevsno - 1; j++) { sabg_lyr[j] = 0.0; }
            sabg_lyr[nlevsno - 1] = sabg_snow * 0.6;
            sabg_lyr[nlevsno] = sabg_snow * 0.4;
          } else {
            for (int j = 0; j <= nlevsno; j++) { sabg_lyr[j] = 0.0; }
            sabg_lyr[nlevsno - snl] = sabg_snow * 0.75;
            sabg_lyr[nlevsno - snl + 1] = sabg_snow * 0.25;
          }
        }
  
        // If shallow snow depth, all solar radiation absorbed in top or top two snow layers
        // to prevent unrealistic timestep soil warming 
  
        if (subgridflag == 0) {
          if (snow_depth < 0.1) {
            if (snl == 0) {
              for (int j = 0; j < nlevsno; j++) { sabg_lyr[j] = 0.0; }
              sabg_lyr[nlevsno] = sabg;
            } else if (snl == 1) {
              for (int j = 0; j < nlevsno - 1; j++) { sabg_lyr[j] = 0.0; }
              sabg_lyr[nlevsno - 1] = sabg;
              sabg_lyr[nlevsno] = 0.0;
            } else {
              for (int j = 0; j <= nlevsno; j++) { sabg_lyr[j] = 0.0; }
              sabg_lyr[nlevsno - snl] = sabg_snow * 0.75;
              sabg_lyr[nlevsno - snl + 1] = sabg_snow * 0.25;
            }
          }
        }
      }
  
      // Error check - This situation should not happen:
      for (int j = 0; j <= nlevsno; j++) { err_sum += sabg_lyr[j]; }
      assert(!(std::abs(err_sum - sabg_snow) > 0.00001)); 
  
      // Diagnostic: shortwave penetrating ground (e.g. top layer)
      if (Land.ltype == istsoil || Land.ltype == istcrop) { sabg_pen = sabg - sabg_lyr[snl]; }
    }
  }


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
  void SurfRadReflected(
    const LandType& Land,
    const double* albd,
    const double* albi,
    const double* forc_solad,
    const double* forc_solai,
    double& fsr) {

    double fsr_vis_d, fsr_nir_d, fsr_vis_i, fsr_nir_i, rvis, rnir;
    // Radiation diagnostics
    if (!Land.urbpoi) {
      // NDVI and reflected solar radiation
      rvis = albd[0] * forc_solad[0] + albi[0] * forc_solai[0];
      rnir = albd[1] * forc_solad[1] + albi[1] * forc_solai[1];
      fsr = rvis + rnir;
    } else {
      // Solar reflected per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
      fsr_vis_d = albd[0] * forc_solad[0];
      fsr_nir_d = albd[1] * forc_solad[1];
      fsr_vis_i = albi[0] * forc_solai[0];
      fsr_nir_i = albi[1] * forc_solai[1];
  
      fsr = fsr_vis_d + fsr_nir_d + fsr_vis_i + fsr_nir_i;
    }
  }

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
nrad                [int] number of canopy layers
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
fsun                [double] sunlit fraction of canopy
*/
  void CanopySunShadeFractions(
    const LandType& Land,
    const int& nrad,
    const double& elai,
    const double* tlai_z,
    const double* fsun_z,
    const double* forc_solad,
    const double* forc_solai,
    const double* fabd_sun_z,
    const double* fabd_sha_z,
    const double* fabi_sun_z,
    const double* fabi_sha_z,
    double* parsun_z,
    double* parsha_z,
    double* laisun_z,
    double* laisha_z,
    double& laisun,
    double& laisha,
    double& fsun) {

    if (!Land.urbpoi) {
      int ipar = 0; // The band index for PAR
      for (int iv = 0; iv < nrad; iv++) {
        parsun_z[iv] = 0.0;
        parsha_z[iv] = 0.0;
        laisun_z[iv] = 0.0;
        laisha_z[iv] = 0.0;
      }
      // Loop over patches to calculate laisun_z and laisha_z for each layer.
      // Derive canopy laisun, laisha, and fsun from layer sums.
      // If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
      // SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
      laisun = 0.0;
      laisha = 0.0;
  
      for (int iv = 0; iv < nrad; iv++) {
        laisun_z[iv] = tlai_z[iv] * fsun_z[iv];
        laisha_z[iv] = tlai_z[iv] * (1.0 - fsun_z[iv]);
        laisun += laisun_z[iv];
        laisha += laisha_z[iv];
      }
      if (elai > 0.0) {
        fsun = laisun / elai;
      } else {
        fsun = 0.0;
      }
      // Absorbed PAR profile through canopy
      // If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
      // are canopy integrated so that layer values equal big leaf values.
      for (int iv = 0; iv < nrad; iv++) {
        parsun_z[iv] = forc_solad[ipar] * fabd_sun_z[iv] + forc_solai[ipar] * fabi_sun_z[iv];
        parsha_z[iv] = forc_solad[ipar] * fabd_sha_z[iv] + forc_solai[ipar] * fabi_sha_z[iv];
      }
    }
  }


};

} // namespace ELM
