
#include "invoke_kernel.hh"
#include "surface_fluxes.h"
#include "surface_fluxes_kokkos.hh"

void ELM::kokkos_surface_fluxes(ELMStateType& S,
                                const double& dtime)
{

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call surface_fluxes kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  const auto soitop = ELMdims::nlevsno;

  auto surf_flux_kernels = ELM_LAMBDA (const int& idx) {
    const auto snotop = soitop - S.snl(idx);
    ELM::surface_fluxes::initial_flux_calc(
        S.Land.urbpoi,
        S.snl(idx),
        S.frac_sno_eff(idx),
        S.frac_h2osfc(idx),
        S.t_h2osfc_bef(idx),
        S.tssbef(idx, snotop),
        S.tssbef(idx, soitop),
        S.t_grnd(idx),
        S.cgrnds(idx),
        S.cgrndl(idx),
        S.eflx_sh_grnd(idx),
        S.qflx_evap_soi(idx),
        S.qflx_ev_snow(idx),
        S.qflx_ev_soil(idx),
        S.qflx_ev_h2osfc(idx));

    ELM::surface_fluxes::update_surface_fluxes(
        S.Land.urbpoi,
        S.do_capsnow(idx),
        S.snl(idx),
        dtime,
        S.t_grnd(idx),
        S.htvp(idx),
        S.frac_sno_eff(idx),
        S.frac_h2osfc(idx),
        S.t_h2osfc_bef(idx),
        S.sabg_soil(idx),
        S.sabg_snow(idx),
        S.dlrad(idx),
        S.frac_veg_nosno(idx),
        S.emg(idx),
        S.forc_lwrad(idx),
        S.tssbef(idx, snotop),
        S.tssbef(idx, soitop),
        S.h2osoi_ice(idx, snotop),
        S.h2osoi_liq(idx, soitop),
        S.eflx_sh_veg(idx),
        S.qflx_evap_veg(idx),
        S.qflx_evap_soi(idx),
        S.eflx_sh_grnd(idx),
        S.qflx_ev_snow(idx),
        S.qflx_ev_soil(idx),
        S.qflx_ev_h2osfc(idx),
        S.eflx_soil_grnd(idx),
        S.eflx_sh_tot(idx),
        S.qflx_evap_tot(idx),
        S.eflx_lh_tot(idx),
        S.qflx_evap_grnd(idx),
        S.qflx_sub_snow(idx),
        S.qflx_dew_snow(idx),
        S.qflx_dew_grnd(idx),
        S.qflx_snwcp_liq(idx),
        S.qflx_snwcp_ice(idx));

    ELM::surface_fluxes::lwrad_outgoing(
        S.Land.urbpoi,
        S.snl(idx),
        S.frac_veg_nosno(idx),
        S.forc_lwrad(idx),
        S.frac_sno_eff(idx),
        S.tssbef(idx, snotop),
        S.tssbef(idx, soitop),
        S.frac_h2osfc(idx),
        S.t_h2osfc_bef(idx),
        S.t_grnd(idx),
        S.ulrad(idx),
        S.emg(idx),
        S.eflx_lwrad_out(idx),
        S.eflx_lwrad_net(idx));

    S.soil_e_balance(idx) = ELM::surface_fluxes::soil_energy_balance(
        S.Land.ctype,
        S.snl(idx),
        S.eflx_soil_grnd(idx),
        S.xmf(idx),
        S.xmf_h2osfc(idx),
        S.frac_h2osfc(idx),
        S.t_h2osfc(idx),
        S.t_h2osfc_bef(idx),
        dtime,
        S.eflx_h2osfc_snow(idx),
        S.frac_sno_eff(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.tssbef, idx, Kokkos::ALL),
        Kokkos::subview(S.fact, idx, Kokkos::ALL));
  }; // end surf_flux lambda
  apply_parallel_for(surf_flux_kernels, "kokkos_surface_fluxes", S.snl.extent(0));
}
