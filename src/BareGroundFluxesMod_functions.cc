/* functions from BareGroundFluxesMod.F90
Compute sensible and latent fluxes and their derivatives with respect
to ground temperature using ground temperatures from previous time step.

Note:
the "displa" displacement height variable used here is local to these functions.
It is not the same as the displa calculated in CanopyTemperature 


*/
#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

void MoninObukIni(
  const double& ur,
  const double& thv,
  const double& dthv,
  const double& zldis,
  const double& z0m,
  double& um,
  double& obu);

/*

INPUTS:
lakpoi           [bool] true => landunit is a lake point
urbpoi           [bool] true => landunit is an urban point
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
ur               [double] wind speed at reference height [m/s]
dth              [double] diff of virtual temp. between ref. height and surface
dqh              [double] diff of humidity between ref. height and surface
*/
void BGFluxInit(
  const bool& lakpoi,
  const bool& urbpoi,
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

  double& ur,
  double& dth,
  double& dqh,
  double& zldis,
  double& um,
  double& obu)
{
  if (!lakpoi && !urbpoi && frac_veg_nosno == 0) {
    ur = std::max(1.0, std::sqrt(forc_u * forc_u + forc_v * forc_v));
    dth = thm - t_grnd;
    dqh = forc_q - qg;
    zldis = forc_hgt_u_patch;

    double dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
    MoninObukIni(ur, thv, dthv, zldis, z0mg, um, obu);
  }
}



/*
forc_hgt_u_patch [double] observational height of wind at pft level [m]
um               [double] wind speed including the stability effect [m/s]
obu              [double] monin-obukhov length (m)
z0m              [double] roughness length over vegetation, momentum [m]



*/

void BGStabilityIteration(
  const int& niters,
  const bool& lakpoi,
  const bool& urbpoi,

const double& forc_hgt_u_patch,
  const double& um,
  )
{  
  double tstar;   // temperature scaling parameter
  double qstar;   // moisture scaling parameter
  double thvstar; // virtual potential temperature scaling parameter
  double wc ;     // convective velocity [m/s]
  double displa = 0.0; // displacement height, 0 for bare ground

  for (int i = 0; i < niters; i++) {

    FrictionVelocityWind(forc_hgt_u_patch, displa, um, obu, z0mg, ustar);
    FrictionVelocityTemperature(forc_hgt_t_patch, displa, obu, z0hg, temp1);
    FrictionVelocityHumidity(forc_hgt_q_patch, forc_hgt_t_patch, displa, obu, z0hg, z0qg, temp1, temp2);

    if (!lakpoi && !urbpoi && frac_veg_nosno == 0) {
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
}


void BGFluxes(
  )
{
  if (!lakpoi && !urbpoi && frac_veg_nosno == 0) {
    
    // Determine aerodynamic resistances
    ram  = 1.0 / (ustar * ustar / um);
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
    ram1 = ram;  // pass value to global variable (driver)

    // Output to pft-level data structures
    // Derivative of fluxes with respect to ground temperature
    cgrnds = raih;
    cgrndl = raiw * dqgdT;
    cgrnd  = cgrnds + htvp * cgrndl;

    // Surface fluxes of momentum, sensible and latent heat
    // using ground temperatures from previous time step
    taux          = -forc_rho * forc_u / ram;
    tauy          = -forc_rho * forc_v / ram;
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
}
