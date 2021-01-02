/*
functions derived from SoilMoistStressMod.F90

These functions get called from CanopyFluxes. 
CLM calc_root_moist_stress() call tree:
calc_root_moist_stress() -> normalize_unfrozen_rootfr() -> array_normalization()
                            > 
                              calc_root_moist_stress_clm45default() -> soil_water_retention_curve%soil_suction()
I remove the original calc_root_moist_stress() and instead call 
normalize_unfrozen_rootfr(), which calls array_normalization() -- note: rootfr_unf[nlevgrnd] needs to be initialized to 0.0
calc_root_moist_stress_clm45default(), which calls soil_water_retention_curve%soil_suction()

*/
#include <algorithm>
#include <cmath>
#include "clm_constants.h"


void array_normalization (double *arr_inout)
{
  double arr_sum = 0.0;

  for (int i = 0; i < nlevgrnd; i++) { arr_sum += arr_inout[i]; }
  for (int i = 0; i < nlevgrnd; i++) {
    if (arr_sum > 0.0) { arr_inout[i] = arr_inout[i] / arr_sum; }
  }
}


/*
DESCRIPTION: compute soil suction potential, negative
INPUTS:
smpsat   [double] minimum soil suction, positive [mm]
s        [double] relative saturation [0-1]
bsw      [double] shape parameter

OUTPUTS:
smp      [double] soil suction, negative, [mm]
*/
void soil_suction(const double& smpsat, const double& s, const double& bsw, double& smp) 
{ smp = -smpsat * pow(s, (-bsw)); }


/*
DESCRIPTION: normalize root fraction for total unfrozen depth

INPUTS:
t_soisno[nlevgrnd+nlevsno] [double] soil temperature (Kelvin)
rootfr[nlevgrnd]           [double] fraction of roots in each soil layer
altmax_indx                [int] index corresponding to maximum active layer depth from current year
altmax_lastyear_indx       [int] index corresponding to maximum active layer depth from prior year

OUTPUTS:
rootfr_unf[nlevgrnd]       [double] root fraction defined for unfrozen layers only
*/
void normalize_unfrozen_rootfr(
  const double *t_soisno,
  const double *rootfr,
  const int& altmax_indx,
  const int& altmax_lastyear_indx,
  double *rootfr_unf) {

  if (perchroot || perchroot_alt) { // Define rootfraction for unfrozen soil only
    if (perchroot_alt) { // use total active layer (defined ass max thaw depth for current and prior year)
      for (int i = 0; i < nlevgrnd; i++) {
        if (i <= std::max(0,std::max(altmax_lastyear_indx, altmax_indx))) {
          rootfr_unf[i] = rootfr[i];
        } else {
          rootfr_unf[i] = 0.0;
        }
      }
    } else { // use instantaneous temperature
      for (int i = 0; i < nlevgrnd; i++) {
        if (t_soisno[nlevsno + i] > tfrz) {
          rootfr_unf[i] = rootfr[i];
        } else {
          rootfr_unf[i] = 0.0;
        }
      }
    }
  }
  array_normalization(rootfr_unf); // normalize the root fraction
}


/*
DESCRIPTION: compute the effective soil porosity

INPUTS:
watsat[nlevgrnd]             [double] volumetric soil water at saturation (porosity)
h2osoi_ice[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
eff_porosity[nlevgrnd]       [double] effective soil porosity
*/
void calc_effective_soilporosity(
  const double *watsat,
  const double *h2osoi_ice,
  const double *dz,
  double *eff_por) {

  double vol_ice;
  for (int i = 0; i < nlevgrnd; i++) {
    //compute the volumetric ice content
    vol_ice = std::min(watsat[i], (h2osoi_ice[nlevsno + i] / (denice * dz[nlevsno + i])));
    // compute the maximum soil space to fill liquid water and air
    eff_por[i] = watsat[i] - vol_ice;
  }
}


/*
DESCRIPTION: compute the volumetric liquid water content

INPUTS:
eff_porosity[nlevgrnd]       [double] effective soil porosity
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
vol_liq[nlevgrnd+nlevsno]    [double] volumetric liquid water content  
*/
void calc_volumetric_h2oliq(
  const double* eff_por,
  const double* h2osoi_liq,
  const double* dz,
  double* vol_liq) {

  for (int i = 0; i < nlevgrnd; i++) {
    // volume of liquid is no greater than effective void space
    vol_liq[nlevsno + i] = std::min(eff_por[i], (h2osoi_liq[nlevsno + i] / (dz[nlevsno + i] * denh2o)));
  }
}


/*
INPUTS:
h2osoi_liqvol[nlevgrnd+nlevsno] [double] liquid volumetric moisture, will be used for BeTR
rootfr[nlevgrnd]                [double] fraction of roots in each soil layer
rootfr_unf[nlevgrnd]            [double] root fraction defined for unfrozen layers only
t_soisno[nlevgrnd+nlevsno]      [double] col soil temperature (Kelvin)
tc_stress                       [double] critical soil temperature for soil water stress (C)
sucsat[nlevgrnd]                [double] minimum soil suction (mm)
watsat[nlevgrnd]                [double] volumetric soil water at saturation (porosity)
h2osoi_vol[nlevgrnd]            [double] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
bsw[nlevgrnd]                   [double] Clapp and Hornberger "b
smpso[numpft]                   [double] soil water potential at full stomatal opening (mm)                    
smpsc[numpft]                   [double] soil water potential at full stomatal closure (mm)
eff_porosity[nlevgrnd]          [double] effective soil porosity


OUTPUTS:
rootr[nlevgrnd]                 [double] effective fraction of roots in each soil layer
rresis[nlevgrnd]                [double] root soil water stress (resistance) by layer (0-1)
btran                           [double] transpiration wetness factor (0 to 1) (integrated soil water stress)
btran2                          [double] integrated soil water stress square
*/
void calc_root_moist_stress_clm45default(
  const int& vtype,
  const double *h2osoi_liqvol,
  const double *rootfr,
  const double *rootfr_unf,
  const double *t_soisno,
  const double& tc_stress,
  const double *sucsat,
  const double *watsat,
  const double *h2osoi_vol,
  const double *bsw,
  const double *smpso,
  const double *smpsc,
  const double *eff_porosity,

  double *rootr,
  double *rresis,
  double& btran,
  double& btran2)
{
  double s_node, smp_node, smp_node_lf;
  const double btran0 = 0.0;

  for (int i = 0; i < nlevgrnd; i++) {
    // Root resistance factors
    // rootr effectively defines the active root fraction in each layer      
    if (h2osoi_liqvol[nlevsno + i] <= 0.0 || t_soisno[nlevsno + i] <= tfrz + tc_stress) {
      rootr[i] = 0.0;
    } else {
      s_node = std::max(h2osoi_liqvol[nlevsno + i] / eff_porosity[i], 0.01);
      soil_suction(sucsat[i], s_node, bsw[i], smp_node);
      smp_node = std::max(smpsc[vtype], smp_node);
      rresis[i] = std::min( (eff_porosity[i]/watsat[i]) * (smp_node - smpsc[vtype]) / (smpso[vtype] - smpsc[vtype]), 1.0);

      if (!perchroot && !perchroot_alt) {
        rootr[i] = rootfr[i] * rresis[i];
      } else {
        rootr[i] = rootfr_unf[i] * rresis[i];
      }

      btran += std::max(rootr[i], 0.0);
      s_node = h2osoi_vol[i] / watsat[i];
      soil_suction(sucsat[i], s_node, bsw[i], smp_node_lf);
      smp_node_lf = std::max(smpsc[vtype], smp_node_lf);
      btran2 += rootfr[i] * std::min((smp_node_lf - smpsc[vtype]) / (smpso[vtype] - smpsc[vtype]), 1.0);
    }
  }
  // Normalize root resistances to get layer contribution to ET
  for (int i = 0; i < nlevgrnd; i++) {
    if (btran > btran0) {
      rootr[i] /= btran;
    } else {
      rootr[i] = 0.0;
    }
  }
}
