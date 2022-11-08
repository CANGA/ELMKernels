
#include "invoke_kernel.hh"
#include "bareground_fluxes.h"
#include "atm_physics.h"
#include "bareground_fluxes_kokkos.hh"

void ELM::kokkos_bareground_fluxes(ELMStateType& S)
{
  size_t ncells = S.snl.extent(0);

  ViewD1 zldis("zldis", ncells);     // reference height "minus" zero displacement height [m]
  ViewD1 displa("displa", ncells);   // displacement height [m]
  ViewD1 dth("dth", ncells);         // diff of virtual temp. between ref. height and surface
  ViewD1 dqh("dqh", ncells);         // diff of humidity between ref. height and surface
  ViewD1 obu("obu", ncells);         // Monin-Obukhov length (m)
  ViewD1 ur("ur", ncells);           // wind speed at reference height [m/s]
  ViewD1 um("um", ncells);           // wind speed including the stablity effect [m/s]
  ViewD1 temp1("temp1", ncells);     // relation for potential temperature profile
  ViewD1 temp2("temp2", ncells);     // relation for specific humidity profile
  ViewD1 temp12m("temp12m", ncells); // relation for potential temperature profile applied at 2-m
  ViewD1 temp22m("temp22m", ncells); // relation for specific humidity profile applied at 2-m
  ViewD1 ustar("ustar", ncells);     // friction velocity [m/s]

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call bareground_fluxes kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto bgflux_kernels = ELM_LAMBDA (const int& idx) {

    double forc_rho = ELM::atm_forcing_physics::derive_forc_rho(S.forc_pbot(idx), S.forc_qbot(idx), S.forc_tbot(idx));

    ELM::bareground_fluxes::initialize_flux(
        S.Land,
        S.frac_veg_nosno(idx),
        S.forc_u(idx),
        S.forc_v(idx),
        S.forc_qbot(idx),
        forc_rho,
        S.forc_hgt_u_patch(idx),
        S.thm(idx),
        S.thv(idx),
        S.t_grnd(idx),
        S.qg(idx),
        S.z0mg(idx),
        S.dlrad(idx),
        S.ulrad(idx),
        zldis(idx),
        displa(idx),
        dth(idx),
        dqh(idx),
        obu(idx),
        ur(idx),
        um(idx));

    ELM::bareground_fluxes::stability_iteration(
        S.Land,
        S.frac_veg_nosno(idx),
        S.forc_hgt_t_patch(idx),
        S.forc_hgt_u_patch(idx),
        S.forc_hgt_q_patch(idx),
        S.z0mg(idx),
        zldis(idx),
        displa(idx),
        dth(idx),
        dqh(idx),
        ur(idx),
        S.forc_qbot(idx),
        S.forc_thbot(idx),
        S.thv(idx),
        S.z0hg(idx),
        S.z0qg(idx),
        obu(idx),
        um(idx),
        temp1(idx),
        temp2(idx),
        temp12m(idx),
        temp22m(idx),
        ustar(idx));

    ELM::bareground_fluxes::compute_flux(
        S.Land,
        S.frac_veg_nosno(idx),
        S.snl(idx),
        forc_rho,
        S.soilbeta(idx),
        S.dqgdT(idx),
        S.htvp(idx),
        S.t_h2osfc(idx),
        S.qg_snow(idx),
        S.qg_soil(idx),
        S.qg_h2osfc(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        S.forc_pbot(idx),
        dth(idx),
        dqh(idx),
        temp1(idx),
        temp2(idx),
        temp12m(idx),
        temp22m(idx),
        ustar(idx),
        S.forc_qbot(idx),
        S.thm(idx),
        S.cgrnds(idx),
        S.cgrndl(idx),
        S.cgrnd(idx),
        S.eflx_sh_grnd(idx),
        S.eflx_sh_tot(idx),
        S.eflx_sh_snow(idx),
        S.eflx_sh_soil(idx),
        S.eflx_sh_h2osfc(idx),
        S.qflx_evap_soi(idx),
        S.qflx_evap_tot(idx),
        S.qflx_ev_snow(idx),
        S.qflx_ev_soil(idx),
        S.qflx_ev_h2osfc(idx),
        S.t_ref2m(idx),
        S.q_ref2m(idx),
        S.rh_ref2m(idx));
  }; // end bgflux lambda
  invoke_kernel(bgflux_kernels, std::make_tuple(S.snl.extent(0)), "kokkos_bareground_fluxes");
}
