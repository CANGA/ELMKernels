
#include "invoke_kernel.hh"
#include "canopy_fluxes.h"
#include "canopy_fluxes_kokkos.hh"

void ELM::kokkos_canopy_fluxes(ELMStateType& S,
                               const double& dtime)
{
    size_t ncols =S.snl.extent(0);

    ViewD1 wtg("wtg", ncols);                  // heat conductance for ground [m/s]
    ViewD1 wtgq("wtgq", ncols);                // latent heat conductance for ground [m/s]
    ViewD1 wtalq("wtalq", ncols);              // normalized latent heat cond. for air and leaf [-]
    ViewD1 wtlq0("wtlq0", ncols);              // normalized latent heat conductance for leaf [-]
    ViewD1 wtaq0("wtaq0", ncols);              // normalized latent heat conductance for air [-]
    ViewD1 wtl0("wtl0", ncols);                // normalized heat conductance for leaf [-]
    ViewD1 wta0("wta0", ncols);                // normalized heat conductance for air [-]
    ViewD1 wtal("wtal", ncols);                // normalized heat conductance for air and leaf [-]
    ViewD1 dayl_factor("dayl_factor", ncols);  // scalar (0-1) for daylength effect on Vcmax
    ViewD1 air("air", ncols);                  // atmos. radiation temporay set
    ViewD1 bir("bir", ncols);                  // atmos. radiation temporay set
    ViewD1 cir("cir", ncols);                  // atmos. radiation temporay set
    ViewD1 el("el", ncols);                    // vapor pressure on leaf surface [pa]
    ViewD1 qsatl("qsatl", ncols);              // leaf specific humidity [kg/kg]
    ViewD1 qsatldT("qsatldT", ncols);          // derivative of "qsatl" on "t_veg"
    ViewD1 taf("taf", ncols);                  // air temperature within canopy space [K]
    ViewD1 qaf("qaf", ncols);                  // humidity of canopy air [kg/kg]
    ViewD1 um("um", ncols);                    // wind speed including the stablity effect [m/s]
    ViewD1 ur("ur", ncols);                    // wind speed at reference height [m/s]
    ViewD1 dth("dth", ncols);                  // diff of virtual temp. between ref. height and surface
    ViewD1 dqh("dqh", ncols);                  // diff of humidity between ref. height and surface
    ViewD1 obu("obu", ncols);                  // Monin-Obukhov length (m)
    ViewD1 zldis("zldis", ncols);              // reference height "minus" zero displacement height [m]
    ViewD1 temp1("temp1", ncols);              // relation for potential temperature profile
    ViewD1 temp2("temp2", ncols);              // relation for specific humidity profile
    ViewD1 temp12m("temp12m", ncols);          // relation for potential temperature profile applied at 2-m
    ViewD1 temp22m("temp22m", ncols);          // relation for specific humidity profile applied at 2-m
    ViewD1 tlbef("tlbef", ncols);              // leaf temperature from previous iteration [K]
    ViewD1 delq("delq", ncols);                // temporary
    ViewD1 dt_veg("dt_veg", ncols);            // change in t_veg, last iteration (Kelvin)

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call canopy_fluxes kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto canflux_kernels = ELM_LAMBDA (const int& idx) {

    double forc_po2 = ELM::atm_forcing_physics::derive_forc_po2(S.forc_pbot(idx));
    double forc_pco2 = ELM::atm_forcing_physics::derive_forc_pco2(S.forc_pbot(idx));
    double forc_rho = ELM::atm_forcing_physics::derive_forc_rho(S.forc_pbot(idx), S.forc_qbot(idx), S.forc_tbot(idx));

    ELM::canopy_fluxes::initialize_flux(
        S.Land,
        S.snl(idx),
        S.frac_veg_nosno(idx),
        S.frac_sno(idx),
        S.forc_hgt_u_patch(idx),
        S.thm(idx),
        S.thv(idx),
        S.max_dayl,
        S.dayl,
        S.altmax_indx(idx),
        S.altmax_lastyear_indx(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.rootfr, idx, Kokkos::ALL),
        S.psn_pft(idx).tc_stress,
        Kokkos::subview(S.sucsat, idx, Kokkos::ALL),
        Kokkos::subview(S.watsat, idx, Kokkos::ALL),
        Kokkos::subview(S.bsw, idx, Kokkos::ALL),
        S.psn_pft(idx).smpso,
        S.psn_pft(idx).smpsc,
        S.elai(idx),
        S.esai(idx),
        S.emv(idx),
        S.emg(idx),
        S.qg(idx),
        S.t_grnd(idx),
        S.forc_tbot(idx),
        S.forc_pbot(idx),
        S.forc_lwrad(idx),
        S.forc_u(idx),
        S.forc_v(idx),
        S.forc_qbot(idx),
        S.forc_thbot(idx),
        S.z0mg(idx),
        S.btran(idx),
        S.displa(idx),
        S.z0mv(idx),
        S.z0hv(idx),
        S.z0qv(idx),
        Kokkos::subview(S.rootr, idx, Kokkos::ALL),
        Kokkos::subview(S.eff_porosity, idx, Kokkos::ALL),
        dayl_factor(idx),
        air(idx),
        bir(idx),
        cir(idx),
        el(idx),
        qsatl(idx),
        qsatldT(idx),
        taf(idx),
        qaf(idx),
        um(idx),
        ur(idx),
        obu(idx),
        zldis(idx),
        delq(idx),
        S.t_veg(idx));

    ELM::canopy_fluxes::stability_iteration(
        S.Land,
        dtime,
        S.snl(idx),
        S.frac_veg_nosno(idx),
        S.frac_sno(idx),
        S.forc_hgt_u_patch(idx),
        S.forc_hgt_t_patch(idx),
        S.forc_hgt_q_patch(idx),
        S.fwet(idx),
        S.fdry(idx),
        S.laisun(idx),
        S.laisha(idx),
        forc_rho,
        S.snow_depth(idx),
        S.soilbeta(idx),
        S.frac_h2osfc(idx),
        S.t_h2osfc(idx),
        S.sabv(idx),
        S.h2ocan(idx),
        S.htop(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        air(idx),
        bir(idx),
        cir(idx),
        ur(idx),
        zldis(idx),
        S.displa(idx),
        S.elai(idx),
        S.esai(idx),
        S.t_grnd(idx),
        S.forc_pbot(idx),
        S.forc_qbot(idx),
        S.forc_thbot(idx),
        S.z0mg(idx),
        S.z0mv(idx),
        S.z0hv(idx),
        S.z0qv(idx),
        S.thm(idx),
        S.thv(idx),
        S.qg(idx),
        S.psn_pft(idx),
        S.nrad(idx),
        S.t10(idx),
        Kokkos::subview(S.tlai_z, idx, Kokkos::ALL),
        S.vcmaxcintsha(idx),
        S.vcmaxcintsun(idx),
        Kokkos::subview(S.parsha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.parsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.laisha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.laisun_z, idx, Kokkos::ALL),
        forc_pco2,
        forc_po2,
        dayl_factor(idx),
        S.btran(idx),
        S.qflx_tran_veg(idx),
        S.qflx_evap_veg(idx),
        S.eflx_sh_veg(idx),
        wtg(idx),
        wtl0(idx),
        wta0(idx),
        wtal(idx),
        el(idx),
        qsatl(idx),
        qsatldT(idx),
        taf(idx),
        qaf(idx),
        um(idx),
        dth(idx),
        dqh(idx),
        obu(idx),
        temp1(idx),
        temp2(idx),
        temp12m(idx),
        temp22m(idx),
        tlbef(idx),
        delq(idx),
        dt_veg(idx),
        S.t_veg(idx),
        wtgq(idx),
        wtalq(idx),
        wtlq0(idx),
        wtaq0(idx));

    ELM::canopy_fluxes::compute_flux(
        S.Land,
        dtime,
        S.snl(idx),
        S.frac_veg_nosno(idx),
        S.frac_sno(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        S.frac_h2osfc(idx),
        S.t_h2osfc(idx),
        S.sabv(idx),
        S.qg_snow(idx),
        S.qg_soil(idx),
        S.qg_h2osfc(idx),
        S.dqgdT(idx),
        S.htvp(idx),
        wtg(idx),
        wtl0(idx),
        wta0(idx),
        wtal(idx),
        air(idx),
        bir(idx),
        cir(idx),
        qsatl(idx),
        qsatldT(idx),
        dth(idx),
        dqh(idx),
        temp1(idx),
        temp2(idx),
        temp12m(idx),
        temp22m(idx),
        tlbef(idx),
        delq(idx),
        dt_veg(idx),
        S.t_veg(idx),
        S.t_grnd(idx),
        S.forc_pbot(idx),
        S.qflx_tran_veg(idx),
        S.qflx_evap_veg(idx),
        S.eflx_sh_veg(idx),
        S.forc_qbot(idx),
        forc_rho,
        S.thm(idx),
        S.emv(idx),
        S.emg(idx),
        S.forc_lwrad(idx),
        wtgq(idx),
        wtalq(idx),
        wtlq0(idx),
        wtaq0(idx),
        S.h2ocan(idx),
        S.eflx_sh_grnd(idx),
        S.eflx_sh_snow(idx),
        S.eflx_sh_soil(idx),
        S.eflx_sh_h2osfc(idx),
        S.qflx_evap_soi(idx),
        S.qflx_ev_snow(idx),
        S.qflx_ev_soil(idx),
        S.qflx_ev_h2osfc(idx),
        S.dlrad(idx),
        S.ulrad(idx),
        S.cgrnds(idx),
        S.cgrndl(idx),
        S.cgrnd(idx),
        S.t_ref2m(idx),
        S.q_ref2m(idx),
        S.rh_ref2m(idx));
  }; // end canflux lambda
  invoke_kernel(canflux_kernels, std::make_tuple(S.snl.extent(0)), "kokkos_canopy_fluxes");
}
