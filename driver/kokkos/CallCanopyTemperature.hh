// CallCanopyTemperature.hh
#pragma once
#include "CanopyTemperature.h"
#include "Kokkos_Core.hpp"
#include "clm_constants.h"
#include "landtype.h"

using ArrayI1 = Kokkos::View<int *>;
using ArrayB1 = Kokkos::View<bool *>;
using ArrayD1 = Kokkos::View<double *>;
using ArrayD2 = Kokkos::View<double **>;

struct CallCanTemp {

  CallCanTemp(ELM::LandType &Land_, ArrayB1 &veg_active_, ArrayI1 &snl_, ArrayI1 &frac_veg_nosno_, ArrayD1 &t_h2osfc_,
              ArrayD1 &t_h2osfc_bef_, ArrayD1 &frac_sno_eff_, ArrayD1 &frac_h2osfc_, ArrayD1 &frac_sno_,
              ArrayD1 &smpmin_, ArrayD1 &forc_q_, ArrayD1 &forc_pbot_, ArrayD1 &forc_th_, ArrayD1 &elai_,
              ArrayD1 &esai_, ArrayD1 &htop_, ArrayD1 &forc_hgt_u_, ArrayD1 &forc_hgt_t_, ArrayD1 &forc_hgt_q_,
              ArrayD1 &z_0_town_, ArrayD1 &z_d_town_, ArrayD1 &forc_t_, ArrayD1 &t_grnd_, ArrayD1 &soilalpha_,
              ArrayD1 &soilalpha_u_, ArrayD1 &soilbeta_, ArrayD1 &qg_snow_, ArrayD1 &qg_soil_, ArrayD1 &qg_,
              ArrayD1 &qg_h2osfc_, ArrayD1 &dqgdT_, ArrayD1 &emg_, ArrayD1 &emv_, ArrayD1 &htvp_, ArrayD1 &z0mg_,
              ArrayD1 &z0hg_, ArrayD1 &z0qg_, ArrayD1 &z0mv_, ArrayD1 &z0hv_, ArrayD1 &z0qv_, ArrayD1 &beta_,
              ArrayD1 &zii_, ArrayD1 &thv_, ArrayD1 &z0m_, ArrayD1 &displa_, ArrayD1 &cgrnd_, ArrayD1 &cgrnds_,
              ArrayD1 &cgrndl_, ArrayD1 &forc_hgt_u_patch_, ArrayD1 &forc_hgt_t_patch_, ArrayD1 &forc_hgt_q_patch_,
              ArrayD1 &thm_, ArrayD1 &eflx_sh_tot_, ArrayD1 &eflx_sh_tot_u_, ArrayD1 &eflx_sh_tot_r_,
              ArrayD1 &eflx_lh_tot_, ArrayD1 &eflx_lh_tot_u_, ArrayD1 &eflx_lh_tot_r_, ArrayD1 &eflx_sh_veg_,
              ArrayD1 &qflx_evap_tot_, ArrayD1 &qflx_evap_veg_, ArrayD1 &qflx_tran_veg_, ArrayD2 &tssbef_,
              ArrayD2 &t_soisno_, ArrayD2 &h2osoi_liq_, ArrayD2 &h2osoi_ice_, ArrayD2 &dz_, ArrayD2 &watsat_,
              ArrayD2 &sucsat_, ArrayD2 &bsw_, ArrayD2 &watdry_, ArrayD2 &watopt_, ArrayD2 &rootfr_road_perv_,
              ArrayD2 &rootr_road_perv_, ArrayD2 &watfc_, ArrayD2 &displar_, ArrayD2 &z0mr_)
      : Land(Land_), veg_active(veg_active_), snl(snl_), frac_veg_nosno(frac_veg_nosno_), t_h2osfc(t_h2osfc_),
        t_h2osfc_bef(t_h2osfc_bef_), frac_sno_eff(frac_sno_eff_), frac_h2osfc(frac_h2osfc_), frac_sno(frac_sno_),
        smpmin(smpmin_), forc_q(forc_q_), forc_pbot(forc_pbot_), forc_th(forc_th_), elai(elai_), esai(esai_),
        htop(htop_), forc_hgt_u(forc_hgt_u_), forc_hgt_t(forc_hgt_t_), forc_hgt_q(forc_hgt_q_), z_0_town(z_0_town_),
        z_d_town(z_d_town_), forc_t(forc_t_), t_grnd(t_grnd_), soilalpha(soilalpha_), soilalpha_u(soilalpha_u_),
        soilbeta(soilbeta_), qg_snow(qg_snow_), qg_soil(qg_soil_), qg(qg_), qg_h2osfc(qg_h2osfc_), dqgdT(dqgdT_),
        emg(emg_), emv(emv_), htvp(htvp_), z0mg(z0mg_), z0hg(z0hg_), z0qg(z0qg_), z0mv(z0mv_), z0hv(z0hv_), z0qv(z0qv_),
        beta(beta_), zii(zii_), thv(thv_), z0m(z0m_), displa(displa_), cgrnd(cgrnd_), cgrnds(cgrnds_), cgrndl(cgrndl_),
        forc_hgt_u_patch(forc_hgt_u_patch_), forc_hgt_t_patch(forc_hgt_t_patch_), forc_hgt_q_patch(forc_hgt_q_patch_),
        thm(thm_), eflx_sh_tot(eflx_sh_tot_), eflx_sh_tot_u(eflx_sh_tot_u_), eflx_sh_tot_r(eflx_sh_tot_r_),
        eflx_lh_tot(eflx_lh_tot_), eflx_lh_tot_u(eflx_lh_tot_u_), eflx_lh_tot_r(eflx_lh_tot_r_),
        eflx_sh_veg(eflx_sh_veg_), qflx_evap_tot(qflx_evap_tot_), qflx_evap_veg(qflx_evap_veg_),
        qflx_tran_veg(qflx_tran_veg_), tssbef(tssbef_), t_soisno(t_soisno_), h2osoi_liq(h2osoi_liq_),
        h2osoi_ice(h2osoi_ice_), dz(dz_), watsat(watsat_), sucsat(sucsat_), bsw(bsw_), watdry(watdry_), watopt(watopt_),
        rootfr_road_perv(rootfr_road_perv_), rootr_road_perv(rootr_road_perv_), watfc(watfc_), displar(displar_),
        z0mr(z0mr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    double qred; // soil surface relative humidity
    double hr;   // relative humidity

    ELM::SaveGroundTemp(Land, t_h2osfc[i], Kokkos::subview(t_soisno, i, Kokkos::ALL), t_h2osfc_bef[i],
                        Kokkos::subview(tssbef, i, Kokkos::ALL));

    ELM::CalculateGroundTemp(Land, snl[i], frac_sno_eff[i], frac_h2osfc[i], t_h2osfc[i],
                             Kokkos::subview(t_soisno, i, Kokkos::ALL), t_grnd[i]);

    ELM::CalculateSoilAlpha(Land, frac_sno[i], frac_h2osfc[i], smpmin[i], Kokkos::subview(h2osoi_liq, i, Kokkos::ALL),
                            Kokkos::subview(h2osoi_ice, i, Kokkos::ALL), Kokkos::subview(dz, i, Kokkos::ALL),
                            Kokkos::subview(t_soisno, i, Kokkos::ALL), Kokkos::subview(watsat, i, Kokkos::ALL),
                            Kokkos::subview(sucsat, i, Kokkos::ALL), Kokkos::subview(bsw, i, Kokkos::ALL),
                            Kokkos::subview(watdry, i, Kokkos::ALL), Kokkos::subview(watopt, i, Kokkos::ALL),
                            Kokkos::subview(rootfr_road_perv, i, Kokkos::ALL),
                            Kokkos::subview(rootr_road_perv, i, Kokkos::ALL), qred, hr, soilalpha[i], soilalpha_u[i]);

    ELM::CalculateSoilBeta(Land, frac_sno[i], frac_h2osfc[i], Kokkos::subview(watsat, i, Kokkos::ALL),
                           Kokkos::subview(watfc, i, Kokkos::ALL), Kokkos::subview(h2osoi_liq, i, Kokkos::ALL),
                           Kokkos::subview(h2osoi_ice, i, Kokkos::ALL), Kokkos::subview(dz, i, Kokkos::ALL),
                           soilbeta[i]);

    ELM::CalculateHumidities(Land, snl[i], forc_q[i], forc_pbot[i], t_h2osfc[i], t_grnd[i], frac_sno[i],
                             frac_sno_eff[i], frac_h2osfc[i], qred, hr, Kokkos::subview(t_soisno, i, Kokkos::ALL),
                             qg_snow[i], qg_soil[i], qg[i], qg_h2osfc[i], dqgdT[i]);

    ELM::GroundProperties(Land, snl[i], frac_sno[i], forc_th[i], forc_q[i], elai[i], esai[i], htop[i],
                          Kokkos::subview(displar, i, Kokkos::ALL), Kokkos::subview(z0mr, i, Kokkos::ALL),
                          Kokkos::subview(h2osoi_liq, i, Kokkos::ALL), Kokkos::subview(h2osoi_ice, i, Kokkos::ALL),
                          emg[i], emv[i], htvp[i], z0mg[i], z0hg[i], z0qg[i], z0mv[i], z0hv[i], z0qv[i], beta[i],
                          zii[i], thv[i], z0m[i], displa[i], cgrnd[i], cgrnds[i], cgrndl[i]);

    ELM::CalculateForcingHeight(Land, veg_active[i], frac_veg_nosno[i], forc_hgt_u[i], forc_hgt_t[i], forc_hgt_q[i],
                                z0m[i], z0mg[i], z_0_town[i], z_d_town[i], forc_t[i], displa[i], forc_hgt_u_patch[i],
                                forc_hgt_t_patch[i], forc_hgt_q_patch[i], thm[i]);

    ELM::InitializeEnergyFluxes(Land, eflx_sh_tot[i], eflx_sh_tot_u[i], eflx_sh_tot_r[i], eflx_lh_tot[i],
                                eflx_lh_tot_u[i], eflx_lh_tot_r[i], eflx_sh_veg[i], qflx_evap_tot[i], qflx_evap_veg[i],
                                qflx_tran_veg[i]);
  }

private:
  ELM::LandType Land;
  ArrayB1 veg_active;
  ArrayI1 snl, frac_veg_nosno;

  ArrayD1 t_h2osfc, t_h2osfc_bef, frac_sno_eff, frac_h2osfc, frac_sno, smpmin, forc_q, forc_pbot, forc_th, elai, esai,
      htop, forc_hgt_u, forc_hgt_t, forc_hgt_q, z_0_town, z_d_town, forc_t, t_grnd, soilalpha, soilalpha_u, soilbeta,
      qg_snow, qg_soil, qg, qg_h2osfc, dqgdT, emg, emv, htvp, z0mg, z0hg, z0qg, z0mv, z0hv, z0qv, beta, zii, thv, z0m,
      displa, cgrnd, cgrnds, cgrndl, forc_hgt_u_patch, forc_hgt_t_patch, forc_hgt_q_patch, thm, eflx_sh_tot,
      eflx_sh_tot_u, eflx_sh_tot_r, eflx_lh_tot, eflx_lh_tot_u, eflx_lh_tot_r, eflx_sh_veg, qflx_evap_tot,
      qflx_evap_veg, qflx_tran_veg;

  ArrayD2 tssbef, t_soisno, h2osoi_liq, h2osoi_ice, dz, watsat, sucsat, bsw, watdry, watopt, rootfr_road_perv,
      rootr_road_perv, watfc, displar, z0mr;
};

void canopyTemperatureInvoke(
    const int &ncells_, ELM::LandType &Land_, ArrayB1 &veg_active_, ArrayI1 &snl_, ArrayI1 &frac_veg_nosno_,
    ArrayD1 &t_h2osfc_, ArrayD1 &t_h2osfc_bef_, ArrayD1 &frac_sno_eff_, ArrayD1 &frac_h2osfc_, ArrayD1 &frac_sno_,
    ArrayD1 &smpmin_, ArrayD1 &forc_q_, ArrayD1 &forc_pbot_, ArrayD1 &forc_th_, ArrayD1 &elai_, ArrayD1 &esai_,
    ArrayD1 &htop_, ArrayD1 &forc_hgt_u_, ArrayD1 &forc_hgt_t_, ArrayD1 &forc_hgt_q_, ArrayD1 &z_0_town_,
    ArrayD1 &z_d_town_, ArrayD1 &forc_t_, ArrayD1 &t_grnd_, ArrayD1 &soilalpha_, ArrayD1 &soilalpha_u_,
    ArrayD1 &soilbeta_, ArrayD1 &qg_snow_, ArrayD1 &qg_soil_, ArrayD1 &qg_, ArrayD1 &qg_h2osfc_, ArrayD1 &dqgdT_,
    ArrayD1 &emg_, ArrayD1 &emv_, ArrayD1 &htvp_, ArrayD1 &z0mg_, ArrayD1 &z0hg_, ArrayD1 &z0qg_, ArrayD1 &z0mv_,
    ArrayD1 &z0hv_, ArrayD1 &z0qv_, ArrayD1 &beta_, ArrayD1 &zii_, ArrayD1 &thv_, ArrayD1 &z0m_, ArrayD1 &displa_,
    ArrayD1 &cgrnd_, ArrayD1 &cgrnds_, ArrayD1 &cgrndl_, ArrayD1 &forc_hgt_u_patch_, ArrayD1 &forc_hgt_t_patch_,
    ArrayD1 &forc_hgt_q_patch_, ArrayD1 &thm_, ArrayD1 &eflx_sh_tot_, ArrayD1 &eflx_sh_tot_u_, ArrayD1 &eflx_sh_tot_r_,
    ArrayD1 &eflx_lh_tot_, ArrayD1 &eflx_lh_tot_u_, ArrayD1 &eflx_lh_tot_r_, ArrayD1 &eflx_sh_veg_,
    ArrayD1 &qflx_evap_tot_, ArrayD1 &qflx_evap_veg_, ArrayD1 &qflx_tran_veg_, ArrayD2 &tssbef_, ArrayD2 &t_soisno_,
    ArrayD2 &h2osoi_liq_, ArrayD2 &h2osoi_ice_, ArrayD2 &dz_, ArrayD2 &watsat_, ArrayD2 &sucsat_, ArrayD2 &bsw_,
    ArrayD2 &watdry_, ArrayD2 &watopt_, ArrayD2 &rootfr_road_perv_, ArrayD2 &rootr_road_perv_, ArrayD2 &watfc_,
    ArrayD2 &displar_, ArrayD2 &z0mr_) {

  CallCanTemp call_cantemp(
      Land_, veg_active_, snl_, frac_veg_nosno_, t_h2osfc_, t_h2osfc_bef_, frac_sno_eff_, frac_h2osfc_, frac_sno_,
      smpmin_, forc_q_, forc_pbot_, forc_th_, elai_, esai_, htop_, forc_hgt_u_, forc_hgt_t_, forc_hgt_q_, z_0_town_,
      z_d_town_, forc_t_, t_grnd_, soilalpha_, soilalpha_u_, soilbeta_, qg_snow_, qg_soil_, qg_, qg_h2osfc_, dqgdT_,
      emg_, emv_, htvp_, z0mg_, z0hg_, z0qg_, z0mv_, z0hv_, z0qv_, beta_, zii_, thv_, z0m_, displa_, cgrnd_, cgrnds_,
      cgrndl_, forc_hgt_u_patch_, forc_hgt_t_patch_, forc_hgt_q_patch_, thm_, eflx_sh_tot_, eflx_sh_tot_u_,
      eflx_sh_tot_r_, eflx_lh_tot_, eflx_lh_tot_u_, eflx_lh_tot_r_, eflx_sh_veg_, qflx_evap_tot_, qflx_evap_veg_,
      qflx_tran_veg_, tssbef_, t_soisno_, h2osoi_liq_, h2osoi_ice_, dz_, watsat_, sucsat_, bsw_, watdry_, watopt_,
      rootfr_road_perv_, rootr_road_perv_, watfc_, displar_, z0mr_);

  Kokkos::parallel_for(ncells_, call_cantemp);
}
