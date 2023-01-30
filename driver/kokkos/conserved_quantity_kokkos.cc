
#include "invoke_kernel.hh"
#include "conserved_quantity_evaluators.h"
#include "conserved_quantity_kokkos.hh"

#include <iostream>

void ELM::kokkos_evaluate_conservation(ELMStateType& S,
                                       const double& dtime)
{

  size_t ncols =S.snl.extent(0);
  ViewD1 dtend_column_h2o("dtend_column_h2o", ncols);
  ViewD1 errh2o("errh2o", ncols);
  ViewD1 errh2osno("errh2osno", ncols);
  ViewD1 dwb("dwb", ncols);
  ViewD1 errsol("errsol", ncols);
  ViewD1 errlon("errlon", ncols);
  ViewD1 errseb("errseb", ncols);
  ViewD1 netrad("netrad", ncols);

  double hydrology_source_sink = 0.0; // hardwired for now
  auto conservation_evaluator_kernels = ELM_LAMBDA (const int& idx) {

    dtend_column_h2o(idx) = 
      ELM::conservation_eval::column_water_mass(
        S.h2ocan(idx), S.h2osno(idx), S.h2osfc(idx),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL));

    dwb(idx) = ELM::conservation_eval::dh2o_dt(
      S.dtbegin_column_h2o(idx), dtend_column_h2o(idx), dtime);

  // need to get hydrology model source_sink
  // should be column integrated in units of mm/s
  errh2o(idx) =
    ELM::conservation_eval::column_water_balance_error(S.dtbegin_column_h2o(idx), dtend_column_h2o(idx), hydrology_source_sink, 
                                   S.forc_rain(idx), S.forc_snow(idx), S.qflx_evap_tot(idx),
                                   S.qflx_snwcp_ice(idx), dtime);

  errh2osno(idx) =
    ELM::conservation_eval::snow_water_balance_error(S.snl(idx), S.qflx_dew_snow(idx),
                                  S.qflx_dew_grnd(idx), S.qflx_sub_snow(idx),
                                  S.qflx_evap_grnd(idx), S.qflx_snow_melt(idx),
                                  S.qflx_snwcp_ice(idx), S.qflx_snwcp_liq(idx),
                                  S.qflx_sl_top_soil(idx), S.frac_sno_eff(idx),
                                  S.qflx_rain_grnd(idx), S.qflx_snow_grnd(idx),
                                  S.qflx_h2osfc_ice(idx),
                                  S.h2osno(idx), S.h2osno_old(idx), dtime, S.do_capsnow(idx));

  errsol(idx) = ELM::conservation_eval::solar_shortwave_balance_error(S.fsa(idx), S.fsr(idx),
                                        Kokkos::subview(S.forc_solad, idx, Kokkos::ALL),
                                        Kokkos::subview(S.forc_solai, idx, Kokkos::ALL));

  errlon(idx) =
    ELM::conservation_eval::solar_longwave_balance_error(S.eflx_lwrad_out(idx), S.eflx_lwrad_net(idx), S.forc_lwrad(idx));

  // integrated SEB error
  errseb(idx) =
    ELM::conservation_eval::surface_energy_balance_error(S.sabv(idx), S.sabg_chk(idx),
                                      S.forc_lwrad(idx), S.eflx_lwrad_out(idx),
                                      S.eflx_sh_tot(idx), S.eflx_lh_tot(idx),
                                      S.eflx_soil_grnd(idx));

  netrad(idx) =
    ELM::conservation_eval::net_radiation(S.fsa(idx), S.eflx_lwrad_net(idx));

  };

  apply_parallel_for(conservation_evaluator_kernels, "conservation_evaluator_kernels", S.snl.extent(0));

  std::cout << "dtend_column_h2o::  " << 
dtend_column_h2o(0) << "\nerrh2o  " <<
errh2o(0) << "\nerrh2osno  " <<
errh2osno(0) << "\ndwb  " <<
dwb(0) << "\nerrsol  " <<
errsol(0) << "\nerrlon  " <<
errlon(0) << "\nerrseb  " <<
errseb(0) << " \nnetrad " <<
netrad(0) << "  " << std::endl;
}

