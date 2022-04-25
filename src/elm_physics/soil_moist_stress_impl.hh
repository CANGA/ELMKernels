

#pragma once

namespace ELM::soil_moist_stress {

ACCELERATED
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

ACCELERATED
double soil_suction(const double& smpsat, const double& s, const double& bsw) { return -smpsat * pow(s, (-bsw)); }

ACCELERATED
double dsuction_dsat(const double& bsw, const double& smp, const double& s) { return -bsw * smp / s; }

template <class ArrayD1>
ACCELERATED
void normalize_unfrozen_rootfr(const ArrayD1 t_soisno, const ArrayD1 rootfr, const int& altmax_indx,
                               const int& altmax_lastyear_indx, double *rootfr_unf) {

  if (perchroot || perchroot_alt) { // Define rootfraction for unfrozen soil only
    if (perchroot_alt) {            // use total active layer (defined as max thaw depth for current and prior year)
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

template <class ArrayD1>
ACCELERATED
void calc_effective_soilporosity(const ArrayD1 watsat, const ArrayD1 h2osoi_ice, const ArrayD1 dz, ArrayD1 eff_por) {

  double vol_ice;
  for (int i = 0; i < nlevgrnd; i++) {
    // compute the volumetric ice content
    vol_ice = std::min(watsat[i], (h2osoi_ice[nlevsno + i] / (denice * dz[nlevsno + i])));
    // compute the maximum soil space to fill liquid water and air
    eff_por[i] = watsat[i] - vol_ice;
  }
}

template <class ArrayD1>
ACCELERATED
void calc_volumetric_h2oliq(const ArrayD1 eff_por, const ArrayD1 h2osoi_liq, const ArrayD1 dz, double *vol_liq) {

  for (int i = 0; i < nlevgrnd; i++) {
    // volume of liquid is no greater than effective void space
    vol_liq[nlevsno + i] = std::min(eff_por[i], (h2osoi_liq[nlevsno + i] / (dz[nlevsno + i] * denh2o)));
  }
}

template <class ArrayD1>
ACCELERATED
void calc_root_moist_stress(const double *h2osoi_liqvol, const ArrayD1 rootfr, const ArrayD1 t_soisno,
                            const double& tc_stress, const ArrayD1 sucsat, const ArrayD1 watsat, const ArrayD1 bsw,
                            const double& smpso, const double& smpsc, const ArrayD1 eff_porosity,
                            const int& altmax_indx, const int& altmax_lastyear_indx, ArrayD1 rootr, double& btran) {
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
      smp_node = soil_suction(sucsat[i], s_node, bsw[i]);
      smp_node = std::max(smpsc, smp_node);
      rresis[i] = std::min((eff_porosity[i] / watsat[i]) * (smp_node - smpsc) / (smpso - smpsc), 1.0);

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

} // namespace ELM::soil_moist_stress
