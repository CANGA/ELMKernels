#pragma once

namespace ELM::conservation_eval {

template <typename ArrayD1>
ACCELERATE
double column_water_mass(const double& h2ocan, const double& h2osno, const double& h2osfc,
                         const ArrayD1 h2osoi_ice, const ArrayD1 h2osoi_liq)
{
  double water = h2ocan + h2osno + h2osfc;
  for (int i = 0; i < ELMdims::nlevgrnd + ELMdims::nlevsno; ++i) {
    water += h2osoi_ice(i) + h2osoi_liq(i);
  }
  return water;
}


ACCELERATE
double dh2o_dt(const double& begwb, const double& endwb, const double& dtime)
{
  return (endwb - begwb) / dtime;
}


ACCELERATE
double column_water_balance_error(const double& begwb, const double& endwb, const double& hydrology_source_sink, 
                                  const double& forc_rain, const double& forc_snow, const double& qflx_evap_tot,
                                  const double& qflx_snwcp_ice, const double& dtime)
{

  return (endwb - begwb) - 
         (forc_rain + forc_snow - hydrology_source_sink - qflx_evap_tot - qflx_snwcp_ice) * dtime;
}


ACCELERATE
double snow_water_balance_error(const int& snl, const double& qflx_dew_snow,
                                const double& qflx_dew_grnd, const double& qflx_sub_snow,
                                const double& qflx_evap_grnd, const double& qflx_snow_melt,
                                const double& qflx_snwcp_ice, const double& qflx_snwcp_liq,
                                const double& qflx_sl_top_soil, const double& frac_sno_eff,
                                const double& qflx_rain_grnd, const double& qflx_snow_grnd,
                                const double& qflx_h2osfc_ice, const double& h2osno,
                                const double& h2osno_old, const double& dtime, const bool do_capsnow)
{
  double snow_sources{0.0}, snow_sinks{0.0};
  if (snl > 0) {
    double snow_sources = qflx_snow_grnd + qflx_rain_grnd + qflx_dew_snow + qflx_dew_grnd;
    double snow_sinks = qflx_sub_snow + qflx_evap_grnd + qflx_snow_melt +
      qflx_snwcp_ice + qflx_snwcp_liq + qflx_sl_top_soil;

    if (do_capsnow) {
      snow_sources = frac_sno_eff *
        (qflx_dew_snow + qflx_dew_grnd ) + qflx_h2osfc_ice + qflx_snow_grnd + qflx_rain_grnd;
      snow_sinks = frac_sno_eff * (qflx_sub_snow + qflx_evap_grnd) +
        qflx_snwcp_ice + qflx_snwcp_liq + qflx_snow_melt + qflx_sl_top_soil;
    } else {
      const double qflx_snow_h2osfc{0.0}; // snow falling onto water surface - add later
      snow_sources = (qflx_snow_grnd - qflx_snow_h2osfc ) + frac_sno_eff *
        (qflx_rain_grnd + qflx_dew_snow + qflx_dew_grnd ) + qflx_h2osfc_ice;
      snow_sinks = frac_sno_eff *
        (qflx_sub_snow + qflx_evap_grnd) + qflx_snow_melt + qflx_sl_top_soil;
    }
    return (h2osno - h2osno_old) - (snow_sources - snow_sinks) * dtime;
  } else {
    return 0.0;
  }
}


template <typename ArrayD1>
ACCELERATE
double solar_shortwave_balance_error(const double& fsa, const double& fsr,
                                     const ArrayD1 forc_solad, const ArrayD1 forc_solai)
{
  //if (!urbpoi)
    return fsa + fsr - (forc_solad(0) + forc_solad(1) + forc_solai(0) + forc_solai(1));
  //else
  //  return 0.0;
}


ACCELERATE
double solar_longwave_balance_error(const double& eflx_lwrad_out, const double& eflx_lwrad_net,
                                    const double& forc_lwrad)
{
  //if (!urbpoi)
    return eflx_lwrad_out - eflx_lwrad_net - forc_lwrad;
  //else
  //  return 0.0;
}


ACCELERATE
double surface_energy_balance_error(const double& sabv, const double& sabg_chk,
                                    const double& forc_lwrad, const double& eflx_lwrad_out,
                                    const double& eflx_sh_tot, const double& eflx_lh_tot,
                                    const double& eflx_soil_grnd)
{
  return sabv + sabg_chk + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot - eflx_soil_grnd;
}


ACCELERATE
double net_radiation(const double& fsa, const double& eflx_lwrad_net)
{
  return fsa - eflx_lwrad_net;
}

} // namespace ELM::conservation_eval
