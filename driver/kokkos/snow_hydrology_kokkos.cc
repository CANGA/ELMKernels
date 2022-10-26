
#include "invoke_kernel.hh"
#include "snow_hydrology.h"

#include "snow_hydrology_kokkos.hh"
#include "transpiration_impl.hh"

void ELM::kokkos_snow_hydrology(ELMStateType& S,
                                ELM::AerosolMasses<ViewD2>& aerosol_masses,
                                ELM::AerosolDataManager<ViewD1>& aerosol_data,
                                const ELM::SnwRdsTable<ViewD3>& snw_rds_table,
                                const ELM::Utils::Date& time_plus_half_dt,
                                const double& dtime)
{

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call snow hydrology kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

  // evaluate change in snow mass due to water movement
  // and snow->soil fluxes
  auto snow_water_kernel = ELM_LAMBDA (const int& idx) {
    ELM::snow::snow_water(
        S.do_capsnow(idx),
        S.snl(idx),
        dtime,
        S.frac_sno_eff(idx),
        S.h2osno(idx),
        S.qflx_sub_snow(idx),
        S.qflx_evap_grnd(idx),
        S.qflx_dew_snow(idx),
        S.qflx_dew_grnd(idx),
        S.qflx_rain_grnd(idx),
        S.qflx_snomelt(idx),
        S.qflx_snow_melt(idx),
        S.qflx_top_soil(idx),
        S.int_snow(idx),
        S.frac_sno(idx),
        S.mflx_neg_snow(idx),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst1, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst2, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst3, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst4, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL));
  }; // end snow_water lambda
  invoke_kernel(snow_water_kernel, std::make_tuple(S.snl.extent(0)), "kokkos_snow_hydrology::snow_water");
    
  // aerosol deposition must be called between these two snow hydrology functions
  // invokes it's own parallel loop
  ELM::aerosols::invoke_aerosol_source(time_plus_half_dt, dtime, S.snl, aerosol_data, aerosol_masses);

  // evolve snowpack - adds/combines/divides/removes snow layers
  // updates snow mesh
  // updates snow water, radius, density
  // calculates transpiration -- this should be moved
  auto snow_update_kernels = ELM_LAMBDA (const int& idx) {
    ELM::snow::aerosol_phase_change(
        S.snl(idx),
        dtime,
        S.qflx_sub_snow(idx),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL));

    // move this!!
    ELM::trans::transpiration(S.veg_active(idx), S.qflx_tran_veg(idx),
        Kokkos::subview(S.rootr, idx, Kokkos::ALL),
        Kokkos::subview(S.qflx_rootsoi, idx, Kokkos::ALL));

    ELM::snow::snow_compaction(S.snl(idx),
        S.Land.ltype,
        dtime,
        S.int_snow(idx),
        S.n_melt(idx),
        S.frac_sno(idx),
        Kokkos::subview(S.imelt, idx, Kokkos::ALL),
        Kokkos::subview(S.swe_old, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.frac_iceold, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL));


    ELM::snow::combine_layers(
        S.Land.urbpoi,
        S.Land.ltype,
        dtime,
        S.snl(idx),
        S.h2osno(idx),
        S.snow_depth(idx),
        S.frac_sno_eff(idx),
        S.frac_sno(idx),
        S.int_snow(idx),
        S.qflx_sl_top_soil(idx),
        S.qflx_snow2topsoi(idx),
        S.mflx_snowlyr_col(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.snw_rds, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst1, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst2, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst3, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst4, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.zsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.zisoi, idx, Kokkos::ALL));


    ELM::snow::divide_layers(
        S.frac_sno(idx),
        S.snl(idx),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.snw_rds, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst1, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst2, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst3, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_masses.mss_dst4, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.zsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.zisoi, idx, Kokkos::ALL));


    ELM::snow::prune_snow_layers(
        S.snl(idx),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.zsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.zisoi, idx, Kokkos::ALL));


    ELM::snow::snow_aging(
        S.do_capsnow(idx),
        S.snl(idx),
        S.frac_sno(idx),
        dtime,
        S.qflx_snwcp_ice(idx),
        S.qflx_snow_grnd(idx),
        S.h2osno(idx),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.qflx_snofrz_lyr, idx, Kokkos::ALL),
        snw_rds_table,
        Kokkos::subview(S.snw_rds, idx, Kokkos::ALL));

  }; // end snow_update_kernels kernels
  invoke_kernel(snow_update_kernels, std::make_tuple(S.snl.extent(0)), "kokkos_snow_hydrology::snow_update");
}
