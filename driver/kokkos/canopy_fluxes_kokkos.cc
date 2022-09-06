
#include "invoke_kernel.hh"
#include "canopy_fluxes.h"

#include "canopy_fluxes_kokkos.hh"

void ELM::kokkos_canopy_fluxes(ELMStateType& S,
                               const double& dtime)
{

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call canopy_fluxes kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto canflux_kernels = ELM_LAMBDA (const int& idx) {
    // thread local
    double wtg = 0.0;         // heat conductance for ground [m/s]
    double wtgq = 0.0;        // latent heat conductance for ground [m/s]
    double wtalq = 0.0;       // normalized latent heat cond. for air and leaf [-]
    double wtlq0 = 0.0;       // normalized latent heat conductance for leaf [-]
    double wtaq0 = 0.0;       // normalized latent heat conductance for air [-]
    double wtl0 = 0.0;        // normalized heat conductance for leaf [-]
    double wta0 = 0.0;        // normalized heat conductance for air [-]
    double wtal = 0.0;        // normalized heat conductance for air and leaf [-]
    double dayl_factor = 0.0; // scalar (0-1) for daylength effect on Vcmax
    double air = 0.0;         // atmos. radiation temporay set
    double bir = 0.0;         // atmos. radiation temporay set
    double cir = 0.0;         // atmos. radiation temporay set
    double el = 0.0;          // vapor pressure on leaf surface [pa]
    double qsatl = 0.0;       // leaf specific humidity [kg/kg]
    double qsatldT = 0.0;     // derivative of "qsatl" on "t_veg"
    double taf = 0.0;         // air temperature within canopy space [K]
    double qaf = 0.0;         // humidity of canopy air [kg/kg]
    double um = 0.0;          // wind speed including the stablity effect [m/s]
    double ur = 0.0;          // wind speed at reference height [m/s]
    double dth = 0.0;         // diff of virtual temp. between ref. height and surface
    double dqh = 0.0;         // diff of humidity between ref. height and surface
    double obu = 0.0;         // Monin-Obukhov length (m)
    double zldis = 0.0;       // reference height "minus" zero displacement height [m]
    double temp1 = 0.0;       // relation for potential temperature profile
    double temp2 = 0.0;       // relation for specific humidity profile
    double temp12m = 0.0;     // relation for potential temperature profile applied at 2-m
    double temp22m = 0.0;     // relation for specific humidity profile applied at 2-m
    double tlbef = 0.0;       // leaf temperature from previous iteration [K]
    double delq = 0.0;        // temporary
    double dt_veg = 0.0;      // change in t_veg, last iteration (Kelvin)

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
        dayl_factor,
        air,
        bir,
        cir,
        el,
        qsatl,
        qsatldT,
        taf,
        qaf,
        um,
        ur,
        obu,
        zldis,
        delq,
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
        S.forc_rho(idx),
        S.snow_depth(idx),
        S.soilbeta(idx),
        S.frac_h2osfc(idx),
        S.t_h2osfc(idx),
        S.sabv(idx),
        S.h2ocan(idx),
        S.htop(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        air,
        bir,
        cir,
        ur,
        zldis,
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
        S.forc_pco2(idx),
        S.forc_po2(idx),
        dayl_factor,
        S.btran(idx),
        S.qflx_tran_veg(idx),
        S.qflx_evap_veg(idx),
        S.eflx_sh_veg(idx),
        wtg,
        wtl0,
        wta0,
        wtal,
        el,
        qsatl,
        qsatldT,
        taf,
        qaf,
        um,
        dth,
        dqh,
        obu,
        temp1,
        temp2,
        temp12m,
        temp22m,
        tlbef,
        delq,
        dt_veg,
        S.t_veg(idx),
        wtgq,
        wtalq,
        wtlq0,
        wtaq0);

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
        wtg,
        wtl0,
        wta0,
        wtal,
        air,
        bir,
        cir,
        qsatl,
        qsatldT,
        dth,
        dqh,
        temp1,
        temp2,
        temp12m,
        temp22m,
        tlbef,
        delq,
        dt_veg,
        S.t_veg(idx),
        S.t_grnd(idx),
        S.forc_pbot(idx),
        S.qflx_tran_veg(idx),
        S.qflx_evap_veg(idx),
        S.eflx_sh_veg(idx),
        S.forc_qbot(idx),
        S.forc_rho(idx),
        S.thm(idx),
        S.emv(idx),
        S.emg(idx),
        S.forc_lwrad(idx),
        wtgq,
        wtalq,
        wtlq0,
        wtaq0,
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
