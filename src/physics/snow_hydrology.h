
#pragma once

#include "elm_constants.h"
#include "snicar_data.h"

#include "compile_options.hh"

namespace ELM::snow {

  template <typename ArrayD1, typename ArrayD3>
  ACCELERATE
  void snow_aging(const bool& do_capsnow,
                 const int& snl,
                 const double& frac_sno,
                 const double& dtime,
                 const double& qflx_snwcp_ice,
                 const double& qflx_snow_grnd,
                 const double& h2osno,
                 const ArrayD1 dz,
                 const ArrayD1 h2osoi_liq,
                 const ArrayD1 h2osoi_ice,
                 const ArrayD1 t_soisno,
                 const ArrayD1 qflx_snofrz_lyr,
                 const SnwRdsTable<ArrayD3>& snw_table,
                 double& snot_top,
                 double& dTdz_top,
                 double& snw_rds_top,
                 double& sno_liq_top,
                 ArrayD1 snw_rds);


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
  ACCELERATE
  void aerosol_phase_change(const int& snl,
                            const double& dtime,
                            const double& qflx_sub_snow,
                            const ArrayD1 h2osoi_liq,
                            const ArrayD1 h2osoi_ice,
                            ArrayD1 mss_bcphi,
                            ArrayD1 mss_bcpho);


  template <typename ArrayI1, typename ArrayD1>
  ACCELERATE
  void snow_compaction(const int& snl,
                       const int& ltype,
                       const double& dtime,
                       const double& int_snow,
                       const double& n_melt,
                       const double frac_sno,
                       const ArrayI1 imelt,
                       const ArrayD1 swe_old,
                       const ArrayD1 h2osoi_liq,
                       const ArrayD1 h2osoi_ice,
                       const ArrayD1 t_soisno,
                       const ArrayD1 frac_iceold,
                       ArrayD1 dz);


  template <typename ArrayD1>
  ACCELERATE
  void combine_layers(const bool& urbpoi,
                      const int& ltype,
                      const double& dtime,
                      int& snl,
                      double& h2osno,
                      double& snow_depth,
                      double& frac_sno_eff,
                      double& frac_sno,
                      double& int_snow,
                      double& qflx_sl_top_soil,
                      double& qflx_snow2topsoi,
                      double& mflx_snowlyr_col,
                      ArrayD1 t_soisno,
                      ArrayD1 h2osoi_ice,
                      ArrayD1 h2osoi_liq,
                      ArrayD1 snw_rds,
                      ArrayD1 mss_bcphi,
                      ArrayD1 mss_bcpho,
                      ArrayD1 mss_dst1,
                      ArrayD1 mss_dst2,
                      ArrayD1 mss_dst3,
                      ArrayD1 mss_dst4,
                      ArrayD1 dz,
                      ArrayD1 z,
                      ArrayD1 zi);


  template <typename ArrayD1>
  ACCELERATE
  void divide_layers(const double& frac_sno,
                     int& snl,
                     ArrayD1 h2osoi_ice,
                     ArrayD1 h2osoi_liq,
                     ArrayD1 t_soisno,
                     ArrayD1 snw_rds,
                     ArrayD1 mss_bcphi,
                     ArrayD1 mss_bcpho,
                     ArrayD1 mss_dst1,
                     ArrayD1 mss_dst2,
                     ArrayD1 mss_dst3,
                     ArrayD1 mss_dst4,
                     ArrayD1 dz,
                     ArrayD1 z,
                     ArrayD1 zi);


  ACCELERATE
  void combine(const double& dz2,
               const double& wliq2,
               const double& wice2,
               const double& t2,
               double& dz,
               double& wliq,
               double& wice,
               double& t);


  template <typename ArrayD1>
  ACCELERATE
  void prune_snow_layers(const int& snl,
                         ArrayD1 h2osoi_ice,
                         ArrayD1 h2osoi_liq,
                         ArrayD1 t_soisno,
                         ArrayD1 dz,
                         ArrayD1 z,
                         ArrayD1 zi);



} // namespace ELM::snow

#include "snow_hydrology_impl.hh"
