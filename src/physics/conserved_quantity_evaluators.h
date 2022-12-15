
#pragma once

#include "elm_constants.h"
#include "compile_options.hh"

namespace ELM::conservation_eval {

  // total integrated column water (kg/m^2)
  // snow + surface water + canopy water + soil water + soil ice
  template <typename ArrayD1>
  ACCELERATE
  double column_water_mass(const double& h2ocan, const double& h2osno, const double& h2osfc,
                           const ArrayD1 h2osoi_ice, const ArrayD1 h2osoi_liq);

  // d water_mass / d time (kg/m^2/s)
  ACCELERATE
  double dh2o_dt(const double& begwb, const double& endwb, const double& dtime);

  // need to get hydrology model source_sink
  // should be column integrated in units of mm/s
  ACCELERATE
  double column_water_balance_error(const double& begwb, const double& endwb, const double& hydrology_source_sink, 
                                    const double& forc_rain, const double& forc_snow, const double& qflx_evap_tot,
                                    const double& qflx_snwcp_ice, const double& dtime);

  ACCELERATE
  double snow_water_balance_error(const int& snl, const double& qflx_dew_snow,
                                  const double& qflx_dew_grnd, const double& qflx_sub_snow,
                                  const double& qflx_evap_grnd, const double& qflx_snow_melt,
                                  const double& qflx_snwcp_ice, const double& qflx_snwcp_liq,
                                  const double& qflx_sl_top_soil, const double& frac_sno_eff,
                                  const double& qflx_rain_grnd, const double& qflx_snow_grnd,
                                  const double& qflx_h2osfc_ice, const double& h2osno,
                                  const double& h2osno_old, const double& dtime, const bool do_capsnow);

  // shortwave solar radiation energy balance error
  // (W/m^2)
  template <typename ArrayD1>
  ACCELERATE
  double solar_shortwave_balance_error(const double& fsa, const double& fsr,
                                       const ArrayD1 forc_solad, const ArrayD1 forc_solai);

  // longwave solar radiation energy balance error
  // (W/m^2)
  ACCELERATE
  double solar_longwave_balance_error(const double& eflx_lwrad_out, const double& eflx_lwrad_net,
                                      const double& forc_lwrad);

  // integrated SEB error
  ACCELERATE
  double surface_energy_balance_error(const double& sabv, const double& sabg_chk,
                                      const double& forc_lwrad, const double& eflx_lwrad_out,
                                      const double& eflx_sh_tot, const double& eflx_lh_tot,
                                      const double& eflx_soil_grnd);
  // net radiation
  // positive downward) (W/m^2)
  ACCELERATE
  double net_radiation(const double& fsa, const double& eflx_lwrad_net);

} // namespace ELM::conservation_eval

#include "conserved_quantity_evaluators_impl.hh"
