/*
functions derived from SoilMoistStressMod.F90

These functions get called from CanopyFluxes.
CLM calc_root_moist_stress() call tree:
calc_root_moist_stress() -> normalize_unfrozen_rootfr() -> array_normalization()
                            >
                              calc_root_moist_stress_clm45default() -> soil_water_retention_curve%soil_suction()
I remove the original calc_root_moist_stress() and instead call
normalize_unfrozen_rootfr(), which calls array_normalization() -- note: rootfr_unf[nlevgrnd] needs to be initialized to
0.0 calc_root_moist_stress_clm45default(), which calls soil_water_retention_curve%soil_suction()

*/
#pragma once

#include "clm_constants.h"
#include <algorithm>
#include <cmath>

namespace ELM {

void array_normalization(double *arr_inout) {
  double arr_sum = 0.0;

  for (int i = 0; i < nlevgrnd; i++) {
    arr_sum += arr_inout[i];
  }
  for (int i = 0; i < nlevgrnd; i++) {
    if (arr_sum > 0.0) {
      arr_inout[i] /= arr_sum;
    }
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
void soil_suction(const double &smpsat, const double &s, const double &bsw, double &smp) {
  smp = -smpsat * pow(s, (-bsw));
}

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
template <class dArray_type>
void normalize_unfrozen_rootfr(const dArray_type t_soisno, const dArray_type rootfr, const int &altmax_indx,
                               const int &altmax_lastyear_indx, double *rootfr_unf) {

  if (perchroot || perchroot_alt) { // Define rootfraction for unfrozen soil only
    if (perchroot_alt) {            // use total active layer (defined ass max thaw depth for current and prior year)
      for (int i = 0; i < nlevgrnd; i++) {
        if (i <= std::max(0, std::max(altmax_lastyear_indx, altmax_indx))) {
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
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
eff_porosity[nlevgrnd]       [double] effective soil porosity
*/
template <class dArray_type>
void calc_effective_soilporosity(const dArray_type watsat, const dArray_type h2osoi_ice, const dArray_type dz,
                                 dArray_type eff_por) {

  double vol_ice;
  for (int i = 0; i < nlevgrnd; i++) {
    // compute the volumetric ice content
    vol_ice = std::min(watsat[i], (h2osoi_ice[nlevsno + i] / (denice * dz[nlevsno + i])));
    // compute the maximum soil space to fill liquid water and air
    eff_por[i] = watsat[i] - vol_ice;
  }
}

/*
DESCRIPTION: compute the volumetric liquid water content

INPUTS:
eff_porosity[nlevgrnd]       [double] effective soil porosity
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
vol_liq[nlevgrnd+nlevsno]    [double] volumetric liquid water content
*/
template <class dArray_type>
void calc_volumetric_h2oliq(const dArray_type eff_por, const dArray_type h2osoi_liq, const dArray_type dz,
                            double *vol_liq) {

  for (int i = 0; i < nlevgrnd; i++) {
    // volume of liquid is no greater than effective void space
    vol_liq[nlevsno + i] = std::min(eff_por[i], (h2osoi_liq[nlevsno + i] / (dz[nlevsno + i] * denh2o)));
  }
}

/*
INPUTS:
h2osoi_liqvol[nlevgrnd+nlevsno] [double] liquid volumetric moisture
rootfr[nlevgrnd]                [double] fraction of roots in each soil layer
t_soisno[nlevgrnd+nlevsno]      [double] col soil temperature (Kelvin)
tc_stress                       [double] critical soil temperature for soil water stress (C)
sucsat[nlevgrnd]                [double] minimum soil suction (mm)
watsat[nlevgrnd]                [double] volumetric soil water at saturation (porosity)
bsw[nlevgrnd]                   [double] Clapp and Hornberger "b
smpso[numpft]                   [double] soil water potential at full stomatal opening (mm)
smpsc[numpft]                   [double] soil water potential at full stomatal closure (mm)
eff_porosity[nlevgrnd]          [double] effective soil porosity
altmax_indx                     [int] index corresponding to maximum active layer depth from current year
altmax_lastyear_indx            [int] index corresponding to maximum active layer depth from prior year

OUTPUTS:
rootr[nlevgrnd]                 [double] effective fraction of roots in each soil layer
btran                           [double] transpiration wetness factor (0 to 1) (integrated soil water stress)
*/
template <class dArray_type>
void calc_root_moist_stress(const int &vtype, const double *h2osoi_liqvol, const dArray_type rootfr,
                            const dArray_type t_soisno, const double &tc_stress, const dArray_type sucsat,
                            const dArray_type watsat, const dArray_type bsw, const dArray_type smpso,
                            const dArray_type smpsc, const dArray_type eff_porosity, const int &altmax_indx,
                            const int &altmax_lastyear_indx, dArray_type rootr, double &btran) {
  double s_node, smp_node;
  const double btran0 = 0.0;
  double rootfr_unf[nlevgrnd] = {0.0}; // unfrozen root fraction
  double rresis[nlevgrnd] = {0.0};     // root soil water stress (resistance) by layer (0-1)

  normalize_unfrozen_rootfr(t_soisno, rootfr, altmax_indx, altmax_lastyear_indx, rootfr_unf);

  for (int i = 0; i < nlevgrnd; i++) {
    // Root resistance factors
    // rootr effectively defines the active root fraction in each layer
    if (h2osoi_liqvol[nlevsno + i] <= 0.0 || t_soisno[nlevsno + i] <= tfrz + tc_stress) {
      rootr[i] = 0.0;
    } else {
      s_node = std::max(h2osoi_liqvol[nlevsno + i] / eff_porosity[i], 0.01);
      soil_suction(sucsat[i], s_node, bsw[i], smp_node);
      smp_node = std::max(smpsc[vtype], smp_node);
      rresis[i] =
          std::min((eff_porosity[i] / watsat[i]) * (smp_node - smpsc[vtype]) / (smpso[vtype] - smpsc[vtype]), 1.0);

      if (!perchroot && !perchroot_alt) {
        rootr[i] = rootfr[i] * rresis[i];
      } else {
        rootr[i] = rootfr_unf[i] * rresis[i];
      }

      btran += std::max(rootr[i], 0.0);
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

} // namespace ELM
