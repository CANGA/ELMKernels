
#include "invoke_kernel.hh"
#include "bareground_fluxes.h"

#include "bareground_fluxes_kokkos.hh"

void ELM::kokkos_bareground_fluxes(ELMStateType& S)
{

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call bareground_fluxes kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto bgflux_kernels = ELM_LAMBDA (const int& idx) {
    // thread local
    double zldis;   // reference height "minus" zero displacement height [m]
    double displa;  // displacement height [m]
    double dth;     // diff of virtual temp. between ref. height and surface
    double dqh;     // diff of humidity between ref. height and surface
    double obu;     // Monin-Obukhov length (m)
    double ur;      // wind speed at reference height [m/s]
    double um;      // wind speed including the stablity effect [m/s]
    double temp1;   // relation for potential temperature profile
    double temp2;   // relation for specific humidity profile
    double temp12m; // relation for potential temperature profile applied at 2-m
    double temp22m; // relation for specific humidity profile applied at 2-m
    double ustar;   // friction velocity [m/s]

    ELM::bareground_fluxes::initialize_flux(
        S.Land,
        S.frac_veg_nosno(idx),
        S.forc_u(idx),
        S.forc_v(idx),
        S.forc_qbot(idx),
        S.forc_thbot(idx),
        S.forc_hgt_u_patch(idx),
        S.thm(idx),
        S.thv(idx),
        S.t_grnd(idx),
        S.qg(idx),
        S.z0mg(idx),
        S.dlrad(idx),
        S.ulrad(idx),
        zldis,
        displa,
        dth,
        dqh,
        obu,
        ur,
        um);

    ELM::bareground_fluxes::stability_iteration(
        S.Land,
        S.frac_veg_nosno(idx),
        S.forc_hgt_t_patch(idx),
        S.forc_hgt_u_patch(idx),
        S.forc_hgt_q_patch(idx),
        S.z0mg(idx),
        zldis,
        displa,
        dth,
        dqh,
        ur,
        S.forc_qbot(idx),
        S.forc_thbot(idx),
        S.thv(idx),
        S.z0hg(idx),
        S.z0qg(idx),
        obu,
        um,
        temp1,
        temp2,
        temp12m,
        temp22m,
        ustar);

    ELM::bareground_fluxes::compute_flux(
        S.Land,
        S.frac_veg_nosno(idx),
        S.snl(idx),
        S.forc_rho(idx),
        S.soilbeta(idx),
        S.dqgdT(idx),
        S.htvp(idx),
        S.t_h2osfc(idx),
        S.qg_snow(idx),
        S.qg_soil(idx),
        S.qg_h2osfc(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        S.forc_pbot(idx),
        dth,
        dqh,
        temp1,
        temp2,
        temp12m,
        temp22m,
        ustar,
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
