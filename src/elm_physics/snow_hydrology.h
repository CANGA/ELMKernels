
#pragma once

#include "elm_constants.h"

#include "kokkos_includes.hh"

namespace ELM::snow {

  template<typename ArrayD1>
  ACCELERATE
  void snow_water(const bool& do_capsnow,
                  const int& snl,
                  const double& dtime,
                  const double& frac_sno_eff,
                  const double& h2osno,
                  const double& qflx_sub_snow,
                  const double& qflx_evap_grnd,
                  const double& qflx_dew_snow,
                  const double& qflx_dew_grnd,
                  const double& qflx_rain_grnd,
                  const double& qflx_snomelt,
                  double& qflx_snow_melt,
                  double& qflx_top_soil,
                  double& int_snow,
                  double& frac_sno,
                  double& mflx_neg_snow,
                  ArrayD1 h2osoi_liq,
                  ArrayD1 h2osoi_ice,
                  ArrayD1 mss_bcphi,
                  ArrayD1 mss_bcpho,
                  ArrayD1 mss_dst1,
                  ArrayD1 mss_dst2,
                  ArrayD1 mss_dst3,
                  ArrayD1 mss_dst4,
                  ArrayD1 dz);

  template <typename ArrayD1>
  void aerosol_phase_change(const int& snl,
                            const double& dtime,
                            const double& qflx_sub_snow,
                            const ArrayD1 h2osoi_liq,
                            const ArrayD1 h2osoi_ice,
                            ArrayD1 mss_bcphi,
                            ArrayD1 mss_bcpho);

} // namespace ELM::snow

#include "snow_hydrology_impl.hh"
