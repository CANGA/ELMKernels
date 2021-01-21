/* functions derived from BareGroundFluxesMod.F90
Compute sensible and latent fluxes and their derivatives with respect
to ground temperature using ground temperatures from previous time step.

Note:
the "displa" displacement height variable used here is local to these functions.
It is not the same as the displa calculated in CanopyTemperature 
*/

#include <algorithm>
#include <cmath>
#include "clm_constants.h"
#include "landtype.h"
#include "frictionvelocity.h"
#include "qsat.h"

class BareGroundFluxes {

private:
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

  // local vars
  double frac_veg_nosno;
  double forc_q;
  double forc_th;
  double forc_hgt_u_patch;
  double thm;
  double thv;
  double z0mg;


public:

/* BareGroundFluxes::InitializeFlux
DESCRIPTION:
initialize fluxes and call MoninObukIni()

INPUTS:
Land                [LandType] struct containing information about landtype 
frac_veg_nosno_in   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
forc_u              [double] atmospheric wind speed in east direction (m/s)
forc_v              [double] atmospheric wind speed in north direction (m/s)
forc_q_in           [double] atmospheric specific humidity (kg/kg)
forc_th_in          [double] atmospheric potential temperature (Kelvin)
forc_hgt_u_patch_in [double] observational height of wind at pft level [m]
thm_in              [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
thv_in              [double] virtual potential temperature (kelvin)
t_grnd              [double] ground temperature (Kelvin)
qg                  [double] ground specific humidity [kg/kg]         
z0mg_in             [double] roughness length over ground, momentum [m]

OUTPUTS:
dlrad            [double] downward longwave radiation below the canopy [W/m2]
ulrad            [double] upward longwave radiation above the canopy [W/m2]  
*/
  void InitializeFlux(
    const LandType& Land,
    const int& frac_veg_nosno_in,
    const double& forc_u,
    const double& forc_v,
    const double& forc_q_in,
    const double& forc_th_in,
    const double& forc_hgt_u_patch_in,
    const double& thm_in,
    const double& thv_in,
    const double& t_grnd,
    const double& qg,
    const double& z0mg_in,
    double& dlrad,
    double& ulrad)
  {
    frac_veg_nosno = frac_veg_nosno_in;
    forc_q = forc_q_in;
    forc_th = forc_th_in;
    forc_hgt_u_patch = forc_hgt_u_patch_in;
    thm = thm_in;
    thv = thv_in;
    z0mg = z0mg_in;

    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {
      ur = std::max(1.0, std::sqrt(forc_u * forc_u + forc_v * forc_v));
      dth = thm - t_grnd;
      dqh = forc_q - qg;
      zldis = forc_hgt_u_patch;
  
      double dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
      displa = 0.0;
      dlrad  = 0.0;
      ulrad  = 0.0;

      MoninObukIni(ur, thv, dthv, zldis, z0mg, um, obu);
    }
  } // InitializeFlux



/* BareGroundFluxes::StabilityIteration()
DESCRIPTION:
calculate Monin-Obukhov length and wind speed

INPUTS:
Land             [LandType] struct containing information about landtype
forc_hgt_t_patch [double] observational height of temperature at pft level [m]
forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
z0mg             [double] roughness length over ground, momentum [m]
zii              [double] convective boundary height [m]
beta             [double] coefficient of convective velocity [-]

OUTPUTS:
z0hg             [double] roughness length over ground, sensible heat [m]
z0qg             [double] roughness length over ground, latent heat [m]
*/
  void StabilityIteration(
    const LandType& Land,
    const double& forc_hgt_t_patch,
    const double& forc_hgt_q_patch,
    const double& z0mg,
    const double& zii,
    const double& beta,
    double& z0hg,
    double& z0qg) {  

    double tstar;   // temperature scaling parameter
    double qstar;   // moisture scaling parameter
    double thvstar; // virtual potential temperature scaling parameter
    double wc ;     // convective velocity [m/s]
    double zeta;    // dimensionless height used in Monin-Obukhov theory
    static const int niters = 3; // number of iterations
  
    for (int i = 0; i < niters; i++) {
  
      FrictionVelocityWind(forc_hgt_u_patch, displa, um, obu, z0mg, ustar);
      FrictionVelocityTemperature(forc_hgt_t_patch, displa, obu, z0hg, temp1);
      FrictionVelocityHumidity(forc_hgt_q_patch, forc_hgt_t_patch, displa, obu, z0hg, z0qg, temp1, temp2);
      FrictionVelocityTemperature2m(obu, z0hg, temp12m);
      FrictionVelocityHumidity2m(obu, z0hg, z0qg, temp12m, temp22m);
  
      if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {
        tstar = temp1 * dth;
        qstar = temp2 * dqh;
        thvstar = tstar * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * qstar;
        z0hg = z0mg / exp(0.13 * pow((ustar * z0mg / 1.5e-5), 0.45));
        z0qg = z0hg;
        zeta = zldis * vkc * grav * thvstar / (pow(ustar, 2.0) * thv); // dimensionless height used in Monin-Obukhov theory
  
        if (zeta >= 0.0) { // stable
          zeta = std::min(2.0, std::max(zeta, 0.01));
          um = std::max(ur, 0.1);
        } else { // unstable
          zeta = std::max(-100.0, std::min(zeta, -0.01));
          wc = beta * pow((-grav * ustar * thvstar * zii / thv), 0.333);
          um = std::sqrt(ur * ur + wc * wc);
        }
        obu = zldis / zeta;
      }
    }
  } // StabilityIteration



/* BareGroundFluxes::ComputeFlux
DESCRIPTION:
calculated bare ground water and energy fluxes

INPUTS:
Land                       [LandType] struct containing information about landtype
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
forc_pbot                  [double]  atmospheric pressure (Pa)

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
  void ComputeFlux(
    const LandType& Land,
    const int& snl,
    const double& forc_rho,
    const double& soilbeta,
    const double& dqgdT,
    const double& htvp,
    const double& t_h2osfc,
    const double& qg_snow,
    const double& qg_soil,
    const double& qg_h2osfc,
    const double *t_soisno,
    const double& forc_pbot,
    double& cgrnds,  
    double& cgrndl,
    double& cgrnd,
    double& eflx_sh_grnd,
    double& eflx_sh_tot,
    double& eflx_sh_snow,
    double& eflx_sh_soil,
    double& eflx_sh_h2osfc,
    double& qflx_evap_soi,
    double& qflx_evap_tot,
    double& qflx_ev_snow,
    double& qflx_ev_soil,
    double& qflx_ev_h2osfc,
    double& t_ref2m,
    double& t_ref2m_r,
    double& q_ref2m,
    double& rh_ref2m,
    double& rh_ref2m_r) {
    
    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {
      
      double rah, raw, raih, raiw;
      double e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT; 
      // Determine aerodynamic resistances
      rah  = 1.0 / (temp1 * ustar);
      raw  = 1.0 / (temp2 * ustar);
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
      cgrnd  = cgrnds + htvp * cgrndl;
  
      // Surface fluxes of momentum, sensible and latent heat
      // using ground temperatures from previous time step
      eflx_sh_grnd  = -raih * dth;
      eflx_sh_tot   = eflx_sh_grnd;
  
      // compute sensible heat fluxes individually
      eflx_sh_snow   = -raih * (thm - t_soisno[nlevsno-snl]);
      eflx_sh_soil   = -raih * (thm - t_soisno[nlevsno]);
      eflx_sh_h2osfc = -raih * (thm - t_h2osfc);
  
      // water fluxes from soil
      qflx_evap_soi  = -raiw * dqh;
      qflx_evap_tot  = qflx_evap_soi;
  
      // compute latent heat fluxes individually
      qflx_ev_snow   = -raiw * (forc_q - qg_snow);
      qflx_ev_soil   = -raiw * (forc_q - qg_soil);
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
  } // ComputeFlux

};
