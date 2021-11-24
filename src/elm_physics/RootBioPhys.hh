#include "elm_constants.h"
#include <assert.h>
#include <cmath>

/*
DESCRIPTION:
compute root profile for soil water uptake using equation from Zeng 2001, J. Hydrometeorology

INPUTS:
vtype [int] vegetation type
roota_par [double] CLM rooting distribution parameter [1/m]
rootb_par [double] CLM rooting distribution parameter [1/m]


OUTPUT:
rootfr_road_perv[nlevgrnd]   [double] fraction of roots in each soil layer
*/
void init_vegrootfr(const int &vtype, double rootfr[nlevgrnd]) {
  // (computing from surface, d is depth in meter): Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
  // Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs given in Zeng et al. (1998).
  if (vtype != noveg) {
    double totrootfr = 0.0;
    for (int i = 0; i < nlevsoi - 1; i++) {
      rootfr[i] = 0.5 * (exp(-roota_par[vtype] * zi[i + nlevsno]) + exp(-rootb_par[vtype] * zi[i + nlevsno]) -
                         exp(-roota_par[vtype] * zi[i + 1 + nlevsno]) - exp(-rootb_par[vtype] * zi[i + 1 + nlevsno]));
      if (i < nlev2bed) {
        totrootfr += rootfr[i];
      }
    }
    rootfr[nlevsoi - 1] =
        0.5 * (exp(-roota_par[vtype] * zi[nlevsoi - 1 + nlevsno]) + exp(-rootb_par[vtype] * zi[nlevsoi - 1 + nlevsno]));
    // Adjust layer root fractions if nlev2bed < nlevsoi
    if (use_var_soil_thick && nlev2bed < nlevsoi) {
      for (int i = 0; i < nlev2bed; i++) {
        rootfr[i] /= totrootfr;
      }
      for (int i = nlev2bed; i < nlevsoi; i++) {
        rootfr[i] = 0.0;
      }
    }
  } else {
    for (int i = 0; i < nlevsoi; i++) {
      rootfr[i] = 0.0;
    }
  }
  for (int i = nlevsoi; i < nlevgrnd; i++) {
    rootfr[i] = 0.0;
  }
}
