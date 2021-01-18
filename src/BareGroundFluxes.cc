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
  double ustar;   // friction velocity [m/s]

public:

/* BareGroundFluxes::InitializeFlux
DESCRIPTION:
initialize fluxes and call MoninObukIni()

INPUTS:
Land             [LandType] struct containing information about landtype 
frac_veg_nosno   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
forc_u           [double] atmospheric wind speed in east direction (m/s)
forc_v           [double] atmospheric wind speed in north direction (m/s)
forc_q           [double] atmospheric specific humidity (kg/kg)
forc_th          [double] atmospheric potential temperature (Kelvin)
forc_hgt_u_patch [double] observational height of wind at pft level [m]
thm              [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
thv              [double] virtual potential temperature (kelvin)
t_grnd           [double] ground temperature (Kelvin)
qg               [double] ground specific humidity [kg/kg]         
z0mg             [double] roughness length over ground, momentum [m]

OUTPUTS:
dlrad            [double] downward longwave radiation below the canopy [W/m2]
ulrad            [double] upward longwave radiation above the canopy [W/m2]  
*/
  void InitializeFlux(
    const LandType& Land,
    const int& frac_veg_nosno,
    const double& forc_u,
    const double& forc_v,
    const double& forc_q,
    const double& forc_th,
    const double& forc_hgt_u_patch,
    const double& thm,
    const double& thv,
    const double& t_grnd,
    const double& qg,
    const double& z0mg,
    double& dlrad,
    double& ulrad)
  {
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
frac_veg_nosno   [int] fraction of vegetation not covered by snow (0 OR 1) [-] 
forc_hgt_u_patch [double] observational height of wind at pft level [m]
forc_hgt_t_patch [double] observational height of temperature at pft level [m]
forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
z0mg             [double] roughness length over ground, momentum [m]
zii              [double] convective boundary height [m]
forc_q           [double] atmospheric specific humidity (kg/kg)
forc_th          [double] atmospheric potential temperature (Kelvin)
beta             [double] coefficient of convective velocity [-]

OUTPUTS:
z0hg             [double] roughness length over ground, sensible heat [m]
z0qg             [double] roughness length over ground, latent heat [m]
*/
  void StabilityIteration(
    const LandType& Land,
    const int& frac_veg_nosno,
    const double& forc_hgt_u_patch,
    const double& forc_hgt_t_patch,
    const double& forc_hgt_q_patch,
    const double& z0mg,
    const double& zii,
    const double& forc_q,
    const double& forc_th,
    const double& thv,
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
frac_veg_nosno             [int] fraction of vegetation not covered by snow (0 OR 1) [-]
snl                        [int] number of snow layers
forc_rho                   [double] density (kg/m**3)
soilbeta                   [double] soil wetness relative to field capacity
dqgdT                      [double] temperature derivative of "qg"
htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]
forc_q                     [double] atmospheric specific humidity (kg/kg)
thm                        [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
t_h2osfc                   [double] surface water temperature
qg_snow                    [double] specific humidity at snow surface [kg/kg]
qg_soil                    [double] specific humidity at soil surface [kg/kg]
qg_h2osfc                  [double] specific humidity at h2osfc surface [kg/kg]
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)

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
*/
  void ComputeFlux(
    const LandType& Land,
    const int& frac_veg_nosno,
    const int& snl,
    const double& forc_rho,
    const double& soilbeta,
    const double& dqgdT,
    const double& htvp,
    const double& forc_q,
    const double& thm,
    const double& t_h2osfc,
    const double& qg_snow,
    const double& qg_soil,
    const double& qg_h2osfc,
    const double *t_soisno,
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
    double& qflx_ev_h2osfc) {
    
    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {
      
      double rah, raw, raih, raiw;
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
    }
  } // ComputeFlux

};
