
#pragma once

namespace ELM::bareground_fluxes {

ACCELERATE
void initialize_flux(const LandType& Land, const int& frac_veg_nosno, const double& forc_u,
                     const double& forc_v, const double& forc_q, const double& forc_th,
                     const double& forc_hgt_u_patch, const double& thm, const double& thv,
                     const double& t_grnd, const double& qg, const double& z0mg, double& dlrad,
                     double& ulrad, double& zldis, double& displa, double& dth, double& dqh, double& obu,
                     double& ur, double& um)
{
  if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {
    ur = std::max(1.0, std::sqrt(forc_u * forc_u + forc_v * forc_v));
    dth = thm - t_grnd;
    dqh = forc_q - qg;
    zldis = forc_hgt_u_patch;

    double dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
    displa = 0.0;
    dlrad = 0.0;
    ulrad = 0.0;

    friction_velocity::monin_obukhov_length(ur, thv, dthv, zldis, z0mg, um, obu);
  }
} // initialize_flux()

ACCELERATE
void stability_iteration(const LandType& Land, const int& frac_veg_nosno, const double& forc_hgt_t_patch,
                         const double& forc_hgt_u_patch, const double& forc_hgt_q_patch, const double& z0mg,
                         const double& zldis, const double& displa, const double& dth, const double& dqh,
                         const double& ur, const double& forc_q, const double& forc_th, const double& thv,
                         double& z0hg, double& z0qg, double& obu, double& um, double& temp1, double& temp2,
                         double& temp12m, double& temp22m, double& ustar)
{
  static constexpr int niters{3};      // number of iterations
  static constexpr double beta{1.0};   // coefficient of convective velocity [-]
  static constexpr double zii{1000.0}; // convective boundary height [m]

  double tstar;                     // temperature scaling parameter
  double qstar;                     // moisture scaling parameter
  double thvstar;                   // virtual potential temperature scaling parameter
  double wc;                        // convective velocity [m/s]
  double zeta;                      // dimensionless height used in Monin-Obukhov theory

  for (int i = 0; i < niters; i++) {

    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {

      // friction velocity calls
      friction_velocity::friction_velocity_wind(forc_hgt_u_patch, displa, um, obu, z0mg, ustar);
      friction_velocity::friction_velocity_temp(forc_hgt_t_patch, displa, obu, z0hg, temp1);
      friction_velocity::friction_velocity_humidity(forc_hgt_q_patch, forc_hgt_t_patch, displa, obu, z0hg, z0qg, temp1,
                                                    temp2);
      friction_velocity::friction_velocity_temp2m(obu, z0hg, temp12m);
      friction_velocity::friction_velocity_humidity2m(obu, z0hg, z0qg, temp12m, temp22m);

      tstar = temp1 * dth;
      qstar = temp2 * dqh;
      thvstar = tstar * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * qstar;
      z0hg = z0mg / exp(0.13 * pow((ustar * z0mg / 1.5e-5), 0.45));
      z0qg = z0hg;
      zeta =
          zldis * ELMconst::VKC() * ELMconst::GRAV() * thvstar / (pow(ustar, 2.0) * thv); // dimensionless height used in Monin-Obukhov theory

      if (zeta >= 0.0) { // stable
        zeta = std::min(2.0, std::max(zeta, 0.01));
        um = std::max(ur, 0.1);
      } else { // unstable
        zeta = std::max(-100.0, std::min(zeta, -0.01));
        wc = beta * pow((-ELMconst::GRAV() * ustar * thvstar * zii / thv), 0.333);
        um = std::sqrt(ur * ur + wc * wc);
      }
      obu = zldis / zeta;
    }
  }
} // stability_iteration()

template <typename ArrayD1>
ACCELERATE
void compute_flux(const LandType& Land, const int& frac_veg_nosno, const int& snl, const double& forc_rho,
                  const double& soilbeta, const double& dqgdT, const double& htvp, const double& t_h2osfc,
                  const double& qg_snow, const double& qg_soil, const double& qg_h2osfc, const ArrayD1 t_soisno,
                  const double& forc_pbot, const double& dth, const double& dqh, const double& temp1,
                  const double& temp2, const double& temp12m, const double& temp22m, const double& ustar,
                  const double& forc_q, const double& thm, double& cgrnds, double& cgrndl, double& cgrnd,
                  double& eflx_sh_grnd, double& eflx_sh_tot, double& eflx_sh_snow, double& eflx_sh_soil,
                  double& eflx_sh_h2osfc, double& qflx_evap_soi, double& qflx_evap_tot, double& qflx_ev_snow,
                  double& qflx_ev_soil, double& qflx_ev_h2osfc, double& t_ref2m, double& q_ref2m,
                  double& rh_ref2m)
{
  using ELMdims::nlevsno;
  
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
    raih = forc_rho * ELMconst::CPAIR() / rah;

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
    eflx_sh_snow = -raih * (thm - t_soisno(nlevsno() - snl));
    eflx_sh_soil = -raih * (thm - t_soisno(nlevsno()));
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
    qsat(t_ref2m, forc_pbot, e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT);

    rh_ref2m = std::min(100.0, (q_ref2m / qsat_ref2m * 100.0));

    //if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
    //  // t_ref2m_r  [double] Rural 2 m height surface air temperature (Kelvin)
    //  // rh_ref2m_r [double] Rural 2 m height surface relative humidity (%)

    //  rh_ref2m_r = rh_ref2m;
    //  t_ref2m_r = t_ref2m;
    //}
  }
} // compute_flux()

} // namespace ELM::bareground_fluxes
