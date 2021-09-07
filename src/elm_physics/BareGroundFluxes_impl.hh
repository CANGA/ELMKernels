#pragma once

#include "QSat.hh"

namespace ELM {

template <class ArrayD1>
void ComputeFlux_BG(const LandType &Land, const int &frac_veg_nosno, const int &snl, const double &forc_rho,
                    const double &soilbeta, const double &dqgdT, const double &htvp, const double &t_h2osfc,
                    const double &qg_snow, const double &qg_soil, const double &qg_h2osfc, const ArrayD1 t_soisno,
                    const double &forc_pbot, const double &dth, const double &dqh, const double &temp1,
                    const double &temp2, const double &temp12m, const double &temp22m, const double &ustar,
                    const double &forc_q, const double &thm, double &cgrnds, double &cgrndl, double &cgrnd,
                    double &eflx_sh_grnd, double &eflx_sh_tot, double &eflx_sh_snow, double &eflx_sh_soil,
                    double &eflx_sh_h2osfc, double &qflx_evap_soi, double &qflx_evap_tot, double &qflx_ev_snow,
                    double &qflx_ev_soil, double &qflx_ev_h2osfc, double &t_ref2m, double &t_ref2m_r, double &q_ref2m,
                    double &rh_ref2m, double &rh_ref2m_r) {

  if (!Land.lakpoi) {
    // Initial set for calculation
    cgrnd = 0.0;
    cgrnds = 0.0;
    cgrndl = 0.0;
  }

  if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {

    double rah, raw, raih, raiw;
    double e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT;
    // Determine aerodynamic resistances
    rah = 1.0 / (temp1 * ustar);
    raw = 1.0 / (temp2 * ustar);
    raih = forc_rho * cpair / rah;

    // Soil evaporation resistance - changed by K.Sakaguchi. Soilbeta is used for evaporation
    if (dqh > 0.0) { // dew  (beta is not applied, just like rsoil used to be)
      raiw = forc_rho / raw;
    } else {
      // Lee and Pielke 1992 beta is applied
      raiw = soilbeta * forc_rho / raw;
    }

    // Output to pft-level data structures
    // Derivative of fluxes with respect to ground temperature
    cgrnds = raih;
    cgrndl = raiw * dqgdT;
    cgrnd = cgrnds + htvp * cgrndl;

    // Surface fluxes of momentum, sensible and latent heat
    // using ground temperatures from previous time step
    eflx_sh_grnd = -raih * dth;
    eflx_sh_tot = eflx_sh_grnd;

    // compute sensible heat fluxes individually
    eflx_sh_snow = -raih * (thm - t_soisno[nlevsno - snl]);
    eflx_sh_soil = -raih * (thm - t_soisno[nlevsno]);
    eflx_sh_h2osfc = -raih * (thm - t_h2osfc);

    // water fluxes from soil
    qflx_evap_soi = -raiw * dqh;
    qflx_evap_tot = qflx_evap_soi;

    // compute latent heat fluxes individually
    qflx_ev_snow = -raiw * (forc_q - qg_snow);
    qflx_ev_soil = -raiw * (forc_q - qg_soil);
    qflx_ev_h2osfc = -raiw * (forc_q - qg_h2osfc);

    // 2 m height air temperature
    t_ref2m = thm + temp1 * dth * (1.0 / temp12m - 1.0 / temp1);

    // 2 m height specific humidity
    q_ref2m = forc_q + temp2 * dqh * (1.0 / temp22m - 1.0 / temp2);

    // 2 m height relative humidity
    QSat(t_ref2m, forc_pbot, e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT);

    rh_ref2m = std::min(100.0, (q_ref2m / qsat_ref2m * 100.0));

    if (Land.ltype == istsoil || Land.ltype == istcrop) {
      rh_ref2m_r = rh_ref2m;
      t_ref2m_r = t_ref2m;
    }
  }
} // ComputeFlux_BG

} // namespace ELM
