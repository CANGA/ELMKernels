#include "BareGroundFluxes.h"
#include "ELMConstants.h"
#include "Kokkos_Core.hpp"
#include "LandType.h"

using ArrayD1 = Kokkos::View<double *>;
using ArrayI1 = Kokkos::View<int *>;
using ArrayD2 = Kokkos::View<double **>;

struct CallBareGroundFluxes {

  CallBareGroundFluxes(ELM::LandType &Land_, ArrayI1 &frac_veg_nosno_, ArrayD1 &forc_u_, ArrayD1 &forc_v_,
                       ArrayD1 &forc_q_, ArrayD1 &forc_th_, ArrayD1 &forc_hgt_u_patch_, ArrayD1 &thm_, ArrayD1 &thv_,
                       ArrayD1 &t_grnd_, ArrayD1 &qg_, ArrayD1 &dlrad_, ArrayD1 &ulrad_, ArrayD1 &forc_hgt_t_patch_,
                       ArrayD1 &forc_hgt_q_patch_, ArrayD1 &z0mg_, ArrayD1 &z0hg_,
                       ArrayD1 &z0qg_, ArrayI1 &snl_, ArrayD1 &forc_rho_, ArrayD1 &soilbeta_, ArrayD1 &dqgdT_,
                       ArrayD1 &htvp_, ArrayD1 &t_h2osfc_, ArrayD1 &qg_snow_, ArrayD1 &qg_soil_, ArrayD1 &qg_h2osfc_,
                       ArrayD2 &t_soisno_, ArrayD1 &forc_pbot_, ArrayD1 &cgrnds_, ArrayD1 &cgrndl_, ArrayD1 &cgrnd_,
                       ArrayD1 &eflx_sh_grnd_, ArrayD1 &eflx_sh_tot_, ArrayD1 &eflx_sh_snow_, ArrayD1 &eflx_sh_soil_,
                       ArrayD1 &eflx_sh_h2osfc_, ArrayD1 &qflx_evap_soi_, ArrayD1 &qflx_evap_tot_,
                       ArrayD1 &qflx_ev_snow_, ArrayD1 &qflx_ev_soil_, ArrayD1 &qflx_ev_h2osfc_, ArrayD1 &t_ref2m_,
                       ArrayD1 &t_ref2m_r_, ArrayD1 &q_ref2m_, ArrayD1 &rh_ref2m_, ArrayD1 &rh_ref2m_r_)
      : Land(Land_), frac_veg_nosno(frac_veg_nosno_), forc_u(forc_u_), forc_v(forc_v_), forc_q(forc_q_),
        forc_th(forc_th_), forc_hgt_u_patch(forc_hgt_u_patch_), thm(thm_), thv(thv_), t_grnd(t_grnd_), qg(qg_),
        z0mg(z0mg_), dlrad(dlrad_), ulrad(ulrad_), forc_hgt_t_patch(forc_hgt_t_patch_),
        forc_hgt_q_patch(forc_hgt_q_patch_), z0hg(z0hg_), z0qg(z0qg_), snl(snl_),
        forc_rho(forc_rho_), soilbeta(soilbeta_), dqgdT(dqgdT_), htvp(htvp_), t_h2osfc(t_h2osfc_), qg_snow(qg_snow_),
        qg_soil(qg_soil_), qg_h2osfc(qg_h2osfc_), t_soisno(t_soisno_), forc_pbot(forc_pbot_), cgrnds(cgrnds_),
        cgrndl(cgrndl_), cgrnd(cgrnd_), eflx_sh_grnd(eflx_sh_grnd_), eflx_sh_tot(eflx_sh_tot_),
        eflx_sh_snow(eflx_sh_snow_), eflx_sh_soil(eflx_sh_soil_), eflx_sh_h2osfc(eflx_sh_h2osfc_),
        qflx_evap_soi(qflx_evap_soi_), qflx_evap_tot(qflx_evap_tot_), qflx_ev_snow(qflx_ev_snow_),
        qflx_ev_soil(qflx_ev_soil_), qflx_ev_h2osfc(qflx_ev_h2osfc_), t_ref2m(t_ref2m_), t_ref2m_r(t_ref2m_r_),
        q_ref2m(q_ref2m_), rh_ref2m(rh_ref2m_), rh_ref2m_r(rh_ref2m_r_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
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

    InitializeFlux_BG(Land, frac_veg_nosno[i], forc_u[i], forc_v[i], forc_q[i], forc_th[i], forc_hgt_u_patch[i], thm[i],
                      thv[i], t_grnd[i], qg[i], z0mg[i], dlrad[i], ulrad[i], zldis, displa, dth, dqh, obu, ur, um);

    StabilityIteration_BG(Land, frac_veg_nosno[i], forc_hgt_t_patch[i], forc_hgt_u_patch[i], forc_hgt_q_patch[i],
                          z0mg[i], zldis, displa, dth, dqh, ur, forc_q[i], forc_th[i], thv[i], z0hg[i],
                          z0qg[i], obu, um, temp1, temp2, temp12m, temp22m, ustar);

    ComputeFlux_BG(Land, frac_veg_nosno[i], snl[i], forc_rho[i], soilbeta[i], dqgdT[i], htvp[i], t_h2osfc[i],
                   qg_snow[i], qg_soil[i], qg_h2osfc[i], Kokkos::subview(t_soisno, i, Kokkos::ALL), forc_pbot[i], dth,
                   dqh, temp1, temp2, temp12m, temp22m, ustar, forc_q[i], thm[i], cgrnds[i], cgrndl[i], cgrnd[i],
                   eflx_sh_grnd[i], eflx_sh_tot[i], eflx_sh_snow[i], eflx_sh_soil[i], eflx_sh_h2osfc[i],
                   qflx_evap_soi[i], qflx_evap_tot[i], qflx_ev_snow[i], qflx_ev_soil[i], qflx_ev_h2osfc[i], t_ref2m[i],
                   t_ref2m_r[i], q_ref2m[i], rh_ref2m[i], rh_ref2m_r[i]);
  }

private:
  ELM::LandType Land;
  ArrayI1 frac_veg_nosno, snl;

  ArrayD1 forc_u, forc_v, forc_q, forc_th, forc_hgt_u_patch, thm, thv, t_grnd, qg, dlrad, ulrad, forc_hgt_t_patch,
      forc_hgt_q_patch, z0mg, z0hg, z0qg, forc_rho, soilbeta, dqgdT, htvp, t_h2osfc, qg_snow, qg_soil,
      qg_h2osfc, forc_pbot, cgrnds, cgrndl, cgrnd, eflx_sh_grnd, eflx_sh_tot, eflx_sh_snow, eflx_sh_soil,
      eflx_sh_h2osfc, qflx_evap_soi, qflx_evap_tot, qflx_ev_snow, qflx_ev_soil, qflx_ev_h2osfc, t_ref2m, t_ref2m_r,
      q_ref2m, rh_ref2m, rh_ref2m_r;

  ArrayD2 t_soisno;
}; // CallBareGroundFluxes

void bareGroundFluxesInvoke(const int &ncells_, ELM::LandType &Land_, ArrayI1 &frac_veg_nosno_, ArrayD1 &forc_u_,
                            ArrayD1 &forc_v_, ArrayD1 &forc_q_, ArrayD1 &forc_th_, ArrayD1 &forc_hgt_u_patch_,
                            ArrayD1 &thm_, ArrayD1 &thv_, ArrayD1 &t_grnd_, ArrayD1 &qg_, ArrayD1 &z0mg_,
                            ArrayD1 &dlrad_, ArrayD1 &ulrad_, ArrayD1 &forc_hgt_t_patch_, ArrayD1 &forc_hgt_q_patch_,
                            ArrayD1 &z0hg_, ArrayD1 &z0qg_, ArrayI1 &snl_,
                            ArrayD1 &forc_rho_, ArrayD1 &soilbeta_, ArrayD1 &dqgdT_, ArrayD1 &htvp_, ArrayD1 &t_h2osfc_,
                            ArrayD1 &qg_snow_, ArrayD1 &qg_soil_, ArrayD1 &qg_h2osfc_, ArrayD2 &t_soisno_,
                            ArrayD1 &forc_pbot_, ArrayD1 &cgrnds_, ArrayD1 &cgrndl_, ArrayD1 &cgrnd_,
                            ArrayD1 &eflx_sh_grnd_, ArrayD1 &eflx_sh_tot_, ArrayD1 &eflx_sh_snow_,
                            ArrayD1 &eflx_sh_soil_, ArrayD1 &eflx_sh_h2osfc_, ArrayD1 &qflx_evap_soi_,
                            ArrayD1 &qflx_evap_tot_, ArrayD1 &qflx_ev_snow_, ArrayD1 &qflx_ev_soil_,
                            ArrayD1 &qflx_ev_h2osfc_, ArrayD1 &t_ref2m_, ArrayD1 &t_ref2m_r_, ArrayD1 &q_ref2m_,
                            ArrayD1 &rh_ref2m_, ArrayD1 &rh_ref2m_r_) {

  CallBareGroundFluxes call_bgf(
      Land_, frac_veg_nosno_, forc_u_, forc_v_, forc_q_, forc_th_, forc_hgt_u_patch_, thm_, thv_, t_grnd_, qg_, dlrad_,
      ulrad_, forc_hgt_t_patch_, forc_hgt_q_patch_, z0mg_, z0hg_, z0qg_, snl_, forc_rho_, soilbeta_,
      dqgdT_, htvp_, t_h2osfc_, qg_snow_, qg_soil_, qg_h2osfc_, t_soisno_, forc_pbot_, cgrnds_, cgrndl_, cgrnd_,
      eflx_sh_grnd_, eflx_sh_tot_, eflx_sh_snow_, eflx_sh_soil_, eflx_sh_h2osfc_, qflx_evap_soi_, qflx_evap_tot_,
      qflx_ev_snow_, qflx_ev_soil_, qflx_ev_h2osfc_, t_ref2m_, t_ref2m_r_, q_ref2m_, rh_ref2m_, rh_ref2m_r_);

  Kokkos::parallel_for(ncells_, call_bgf);
}
