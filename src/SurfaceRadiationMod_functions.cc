/*
Layer based radiative flux variables are indexed in CLM as var(-nlevsno+1:1). 
We use a simple array declared as var[6]. The index mapping is described below.
This mapping will be consistently used in all similar situations throughout the code.

CLM variable indices      This code's variable indices
|-4|   top snow layer       |5|
 --                         ---
|-3|   snow layer           |4|
 --                         ---
|-2|   snow layer           |3|
 --                         ---
|-1|   snow layer           |2|
 --                         ---
| 0|   bottom snow layer    |1| 
-----  ground interface  -------
| 1|   top soil layer       |0|

Call sequence: 
SurfRadZeroFluxes()
SurfRadAbsorbed()
SurfRadLayers()
SurfRadReflected
*/

#include <cmath>
#include <assert.h>
#include "clm_constants.hh"

/* 
DESCRIPTION: zero out fluxes before surface radiation calculations

INPUTS:
ltype [int] landunit type

OUTPUTS:
sabg_soil [double] solar radiation absorbed by soil (W/m**2) 
sabg_snow [double] solar radiation absorbed by snow (W/m**2)
sabg      [double] solar radiation absorbed by ground (W/m**2)
sabv      [double] solar radiation absorbed by vegetation (W/m**2)
fsa       [double] solar radiation absorbed (total) (W/m**2)
fsa_r     [double] rural solar radiation absorbed (total) (W/m**2)
fsun      [double] sunlit fraction of canopy 
sabg_lyr  [double] absorbed radiative flux (pft,lyr) [W/m2] */
void SurfRadZeroFluxes (
  const int& ltype,

  double& sabg_soil,
  double& sabg_snow,
  double& sabg,
  double& sabv,
  double& fsa,
  double& fsa_r,
  double& fsun,
  double* sabg_lyr) 
{
  // Initialize fluxes
  if (ltype < isturb_MIN) {
    sabg_soil = 0.0;
    sabg_snow = 0.0;
    sabg      = 0.0;
    sabv      = 0.0;
    fsa       = 0.0;
    if (ltype == istsoil || ltype == istcrop) { fsa_r = 0.0; }
    for (int j = nlevsno; j >= 0; j--) { sabg_lyr[j] = 0.0; }
  }

  // zero-out fsun for the urban patches
  // the non-urban patches were set prior to this call
  // and split into ed and non-ed specific functions
  if (ltype >= isturb_MIN && ltype <= isturb_MAX) { fsun = 0.0; }
}

/*
DESCRIPTION: calculate solat dlux absorbed by canopy, soil, snow, and ground

INPUTS:
ltype [int] landunit type
nband [int] number of solar radiation waveband classes (same as numrad)
snl   [int] number of snow layers

ftdd       [double] down direct flux below canopy per unit direct flux
ftid       [double] down diffuse flux below canopy per unit direct flux
ftii       [double] down diffuse flux below canopy per unit diffuse flux
forc_solad [double] direct beam radiation (W/m**2)
forc_solai [double] diffuse radiation (W/m**2)    


OUTPUTS:


*/
void SurfRadAbsorbed(
  const int& ltype,
  const int& nband,
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

  double& parveg,
  double& fsa_r,
  double& sabv,
  double& fsa,
  double& sabg,
  double& sabg_soil,
  double& sabg_snow,
  double* cad,
  double* cai,
  double* trd,
  double* tri)
{
  double absrad;

  if (ltype < isturb_MIN) {
    for (int ib = 0; ib < nband; ib++) {

      cad[ib] = forc_solad[ib] * fabd[ib];
      cai[ib] = forc_solai[ib] * fabi[ib];
      sabv += cad[ib] + cai[ib];
      fsa  += cad[ib] + cai[ib];

      if (ib == 0) { parveg = cad[ib] + cai[ib]; }
      if (ltype == istsoil || ltype == istcrop) { fsa_r  += cad[ib] + cai[ib]; }

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

      if (ltype == istsoil || ltype == istcrop) { fsa_r += absrad; }
      if (snl == 0) {
        sabg_snow = sabg;
        sabg_soil = sabg;
      }

      if (subgridflag == 0) {
        sabg_snow = sabg;
        sabg_soil = sabg;
      }
    } // end of nbands
  } // end if not urban
}


void SurfRadLayers(
  const int& ltype,
  const int& snl,
  const int& subgridflag,
  const double& sabg,
  const double& sabg_snow,
  const double& snow_depth,
  const double* flx_absdv,
  const double* flx_absdn,
  const double* flx_absiv,
  const double* flx_absin,
  const double* tri,
  const double* trd,

  double& sabg_snl_sum,
  double& sub_surf_abs_SW,
  double& sabg_pen,
  double* sabg_lyr)
{
  double err_sum = 0.0;

  // compute absorbed flux in each snow layer and top soil layer,
  // based on flux factors computed in the radiative transfer portion of SNICAR.
  if (ltype < isturb_MIN) {

    sabg_snl_sum = 0.0;
    sub_surf_abs_SW = 0.0;

    // CASE1: No snow layers: all energy is absorbed in top soil layer
    if (snl == 0) {
      for (int i = 0; i < nlevsno + 1; i++) { sabg_lyr[i] = 0.0; }
      sabg_lyr[0] = sabg; 
      sabg_snl_sum  = sabg_lyr[0];
    } else { // CASE 2: Snow layers present: absorbed radiation is scaled according to flux factors computed by SNICAR

      for (int i = nlevsno; i >= 0; i--) {
        sabg_lyr[i] = flx_absdv[i] * trd[0] + flx_absdn[i] * trd[1] + flx_absiv[i] * tri[0] + flx_absin[i] * tri[1];
        // summed radiation in active snow layers:
        // if snow layer is at or below snow surface 
        if (i <= snl) { sabg_snl_sum += sabg_lyr[i]; }
        // if snow layer is below surface snow layer accumulate subsurface flux as a diagnostic
        if (i < snl) { sub_surf_abs_SW += + sabg_lyr[i]; }
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
          for (int j = nlevsno; j > 0; j--) { sabg_lyr[j] = 0.0; }
          sabg_lyr[0] = sabg;
        } else if (snl == 1) {
          for (int j = nlevsno; j > 1; j--) { sabg_lyr[j] = 0.0; }
          sabg_lyr[1] = sabg_snow * 0.6;
          sabg_lyr[0] = sabg_snow * 0.4;
        } else {
          for (int j = nlevsno; j >= 0; j--) { sabg_lyr[j] = 0.0; }
          sabg_lyr[snl] = sabg_snow * 0.75;
          sabg_lyr[snl-1] = sabg_snow * 0.25;
        }
      }

      // If shallow snow depth, all solar radiation absorbed in top or top two snow layers
      // to prevent unrealistic timestep soil warming 

      if (subgridflag == 0) {
        if (snow_depth < 0.1) {
          if (snl == 0) {
            for (int j = nlevsno; j > 0; j--) { sabg_lyr[j] = 0.0; }
            sabg_lyr[0] = sabg;
          } else if (snl == 1) {
            for (int j = nlevsno; j > 1; j--) { sabg_lyr[j] = 0.0; }
            sabg_lyr[1] = sabg;
            sabg_lyr[0] = 0.0;
          } else {
            for (int j = nlevsno; j >= 0; j--) { sabg_lyr[j] = 0.0; }
            sabg_lyr[snl] = sabg_snow * 0.75;
            sabg_lyr[snl-1] = sabg_snow * 0.25;
          }
        }
      }
    }

    // Error check - This situation should not happen:
    for (int j = nlevsno; j >= 0; j--) { err_sum += sabg_lyr[j]; }
    assert(!(std::abs(err_sum - sabg_snow) > 0.00001)); 

    // Diagnostic: shortwave penetrating ground (e.g. top layer)
    if (ltype == istsoil || ltype == istcrop) { sabg_pen = sabg - sabg_lyr[snl]; }
  }
}



void SurfRadReflected(
  const int& ltype,
  const double* albd,
  const double* albi,
  const double* forc_solad,
  const double* forc_solai,

  double& rvis,
  double& rnir,
  double& fsr)
{
  double fsr_vis_d, fsr_nir_d, fsr_vis_i, fsr_nir_i;

  // Radiation diagnostics
  if (ltype < isturb_MIN) {
    // NDVI and reflected solar radiation
    rvis = albd[0] * forc_solad[0] + albi[0] * forc_solai[0];
    rnir = albd[1] * forc_solad[1] + albi[1] * forc_solai[1];
    fsr = rvis + rnir;
  }

  if (ltype >= isturb_MIN && ltype <= isturb_MAX) {
    // Solar reflected per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
    fsr_vis_d = albd[0] * forc_solad[0];
    fsr_nir_d = albd[1] * forc_solad[1];
    fsr_vis_i = albi[0] * forc_solai[0];
    fsr_nir_i = albi[1] * forc_solai[1];

    fsr = fsr_vis_d + fsr_nir_d + fsr_vis_i + fsr_nir_i;
  }
}
