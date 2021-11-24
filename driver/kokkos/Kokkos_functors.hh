#include "bareground_fluxes.h"
#include "elm_constants.h"
#include "Kokkos_Core.hpp"
#include "landtype.h"

using ArrayD1 = Kokkos::View<double *>;
using ArrayI1 = Kokkos::View<int *>;
using ArrayD2 = Kokkos::View<double **>;

struct BGFCaller {

  // make private?
  Kokkos::View<ELM::BareGroundFluxes *> bg_fluxes;
  ELM::LandType land;
  ArrayI1 frac_vec_nosno_in, snl;

  ArrayD1 forc_u, forc_v, forc_q_in, forc_th_in, forc_hgt_u_patch_in, thm_in, thv_in, t_grnd, qg, z0mg_in, dlrad, ulrad,
      forc_hgt_t_patch, forc_hgt_q_patch, z0mg, zii, beta, z0hg, z0qg, forc_rho, soilbeta, dqgdT, htvp, t_h2osfc,
      qg_snow, qg_soil, qg_h2osfc, forc_pbot, cgrnds, cgrndl, cgrnd, eflx_sh_grnd, eflx_sh_tot, eflx_sh_snow,
      eflx_sh_soil, eflx_sh_h2osfc, qflx_evap_soi, qflx_evap_tot, qflx_ev_snow, qflx_ev_soil, qflx_ev_h2osfc, t_ref2m,
      t_ref2m_r, q_ref2m, rh_ref2m, rh_ref2m_r;

  ArrayD2 t_soisno;

  BGFCaller(Kokkos::View<ELM::BareGroundFluxes *> &bg_fluxes_, ELM::LandType &land_, ArrayI1 &frac_vec_nosno_in_,
            ArrayD1 &forc_u_, ArrayD1 &forc_v_, ArrayD1 &forc_q_in_, ArrayD1 &forc_th_in_,
            ArrayD1 &forc_hgt_u_patch_in_, ArrayD1 &thm_in_, ArrayD1 &thv_in_, ArrayD1 &t_grnd_, ArrayD1 &qg_,
            ArrayD1 &z0mg_in_, ArrayD1 &dlrad_, ArrayD1 &ulrad_, ArrayD1 &forc_hgt_t_patch_, ArrayD1 &forc_hgt_q_patch_,
            ArrayD1 &z0mg_, ArrayD1 &zii_, ArrayD1 &beta_, ArrayD1 &z0hg_, ArrayD1 &z0qg_, ArrayI1 &snl_,
            ArrayD1 &forc_rho_, ArrayD1 &soilbeta_, ArrayD1 &dqgdT_, ArrayD1 &htvp_, ArrayD1 &t_h2osfc_,
            ArrayD1 &qg_snow_, ArrayD1 &qg_soil_, ArrayD1 &qg_h2osfc_, ArrayD2 &t_soisno_, ArrayD1 &forc_pbot_,
            ArrayD1 &cgrnds_, ArrayD1 &cgrndl_, ArrayD1 &cgrnd_, ArrayD1 &eflx_sh_grnd_, ArrayD1 &eflx_sh_tot_,
            ArrayD1 &eflx_sh_snow_, ArrayD1 &eflx_sh_soil_, ArrayD1 &eflx_sh_h2osfc_, ArrayD1 &qflx_evap_soi_,
            ArrayD1 &qflx_evap_tot_, ArrayD1 &qflx_ev_snow_, ArrayD1 &qflx_ev_soil_, ArrayD1 &qflx_ev_h2osfc_,
            ArrayD1 &t_ref2m_, ArrayD1 &t_ref2m_r_, ArrayD1 &q_ref2m_, ArrayD1 &rh_ref2m_, ArrayD1 &rh_ref2m_r_)
      : bg_fluxes(bg_fluxes_), land(land_), frac_vec_nosno_in(frac_vec_nosno_in_), forc_u(forc_u_), forc_v(forc_v_),
        forc_q_in(forc_q_in_), forc_th_in(forc_th_in_), forc_hgt_u_patch_in(forc_hgt_u_patch_in_), thm_in(thm_in_),
        thv_in(thv_in_), t_grnd(t_grnd_), qg(qg_), z0mg_in(z0mg_in_), dlrad(dlrad_), ulrad(ulrad_),
        forc_hgt_t_patch(forc_hgt_t_patch_), forc_hgt_q_patch(forc_hgt_q_patch_), z0mg(z0mg_), zii(zii_), beta(beta_),
        z0hg(z0hg_), z0qg(z0qg_), snl(snl_), forc_rho(forc_rho_), soilbeta(soilbeta_), dqgdT(dqgdT_), htvp(htvp_),
        t_h2osfc(t_h2osfc_), qg_snow(qg_snow_), qg_soil(qg_soil_), qg_h2osfc(qg_h2osfc_), t_soisno(t_soisno_),
        forc_pbot(forc_pbot_), cgrnds(cgrnds_), cgrndl(cgrndl_), cgrnd(cgrnd_), eflx_sh_grnd(eflx_sh_grnd_),
        eflx_sh_tot(eflx_sh_tot_), eflx_sh_snow(eflx_sh_snow_), eflx_sh_soil(eflx_sh_soil_),
        eflx_sh_h2osfc(eflx_sh_h2osfc_), qflx_evap_soi(qflx_evap_soi_), qflx_evap_tot(qflx_evap_tot_),
        qflx_ev_snow(qflx_ev_snow_), qflx_ev_soil(qflx_ev_soil_), qflx_ev_h2osfc(qflx_ev_h2osfc_), t_ref2m(t_ref2m_),
        t_ref2m_r(t_ref2m_r_), q_ref2m(q_ref2m_), rh_ref2m(rh_ref2m_), rh_ref2m_r(rh_ref2m_r_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int c) const {
    bg_fluxes[c].InitializeFlux(land, frac_vec_nosno_in[c], forc_u[c], forc_v[c], forc_q_in[c], forc_th_in[c],
                                forc_hgt_u_patch_in[c], thm_in[c], thv_in[c], t_grnd[c], qg[c], z0mg_in[c], dlrad[c],
                                ulrad[c]);

    bg_fluxes[c].StabilityIteration(land, forc_hgt_t_patch[c], forc_hgt_q_patch[c], z0mg[c], zii[c], beta[c], z0hg[c],
                                    z0qg[c]);

    bg_fluxes[c].ComputeFlux(land, snl[c], forc_rho[c], soilbeta[c], dqgdT[c], htvp[c], t_h2osfc[c], qg_snow[c],
                             qg_soil[c], qg_h2osfc[c], Kokkos::subview(t_soisno, c, Kokkos::ALL), forc_pbot[c],
                             cgrnds[c], cgrndl[c], cgrnd[c], eflx_sh_grnd[c], eflx_sh_tot[c], eflx_sh_snow[c],
                             eflx_sh_soil[c], eflx_sh_h2osfc[c], qflx_evap_soi[c], qflx_evap_tot[c], qflx_ev_snow[c],
                             qflx_ev_soil[c], qflx_ev_h2osfc[c], t_ref2m[c], t_ref2m_r[c], q_ref2m[c], rh_ref2m[c],
                             rh_ref2m_r[c]);
  }
}; // BGFCaller
