#pragma once

#include "QSat.h"

namespace ELM {
/* ComputeFlux_BG()
DESCRIPTION:
calculated bare ground water and energy fluxes

INPUTS:
Land                       [LandType] struct containing information about landtype
frac_veg_nosno             [int] fraction of vegetation not covered by snow (0 OR 1) [-]
snl                        [int] number of snow layers
forc_rho                   [double] density (kg/m**3)
soilbeta                   [double] soil wetness relative to field capacity
dqgdT                      [double] temperature derivative of "qg"
htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]
t_h2osfc                   [double] surface water temperature
qg_snow                    [double] specific humidity at snow surface [kg/kg]
qg_soil                    [double] specific humidity at soil surface [kg/kg]
qg_h2osfc                  [double] specific humidity at h2osfc surface [kg/kg]
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
forc_pbot                  [double] atmospheric pressure (Pa)
dth                        [double] diff of virtual temp. between ref. height and surface
dqh                        [double] diff of humidity between ref. height and surface
temp1                      [double] relation for potential temperature profile
temp2                      [double] relation for specific humidity profile
temp12m                    [double] relation for potential temperature profile applied at 2-m
temp22m                    [double] relation for specific humidity profile applied at 2-m
ustar                      [double] friction velocity [m/s]
forc_q                     [double] atmospheric specific humidity (kg/kg)
thm                        [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)

OUTPUTS:
cgrnds                     [double] deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
cgrndl                     [double] deriv of soil latent heat flux wrt soil temp [w/m**2/k]
cgrnd                      [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
eflx_sh_grnd               [double] sensible heat flux from ground (W/m**2) [+ to atm]
eflx_sh_tot                [double] total sensible heat flux (W/m**2) [+ to atm]
eflx_sh_snow               [double] sensible heat flux from snow (W/m**2) [+ to atm]
eflx_sh_soil               [double] sensible heat flux from soil (W/m**2) [+ to atm]
eflx_sh_h2osfc             [double] sensible heat flux from soil (W/m**2) [+ to atm]
qflx_evap_soi              [double] soil evaporation (mm H2O/s) (+ = to atm)
qflx_evap_tot              [double] qflx_evap_soi + qflx_evap_can + qflx_tran_veg
qflx_ev_snow               [double] evaporation flux from snow (W/m**2) [+ to atm]
qflx_ev_soil               [double] evaporation flux from soil (W/m**2) [+ to atm]
qflx_ev_h2osfc             [double] evaporation flux from h2osfc (W/m**2) [+ to atm]
t_ref2m                    [double]  2 m height surface air temperature (Kelvin)
t_ref2m_r                  [double]  Rural 2 m height surface air temperature (Kelvin)
q_ref2m                    [double]  2 m height surface specific humidity (kg/kg)
rh_ref2m_r                 [double]  Rural 2 m height surface relative humidity (%)
rh_ref2m                   [double]  2 m height surface relative humidity (%)
*/
template <class dArray_type>
void ComputeFlux_BG(const LandType &Land, const int &frac_veg_nosno, const int &snl, const double &forc_rho,
                    const double &soilbeta, const double &dqgdT, const double &htvp, const double &t_h2osfc,
                    const double &qg_snow, const double &qg_soil, const double &qg_h2osfc, const dArray_type t_soisno,
                    const double &forc_pbot, const double &dth, const double &dqh, const double &temp1,
                    const double &temp2, const double &temp12m, const double &temp22m, const double &ustar,
                    const double &forc_q, const double &thm, double &cgrnds, double &cgrndl, double &cgrnd,
                    double &eflx_sh_grnd, double &eflx_sh_tot, double &eflx_sh_snow, double &eflx_sh_soil,
                    double &eflx_sh_h2osfc, double &qflx_evap_soi, double &qflx_evap_tot, double &qflx_ev_snow,
                    double &qflx_ev_soil, double &qflx_ev_h2osfc, double &t_ref2m, double &t_ref2m_r, double &q_ref2m,
                    double &rh_ref2m, double &rh_ref2m_r) {

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
