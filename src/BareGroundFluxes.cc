/* functions derived from BareGroundFluxesMod.F90
   Compute sensible and latent fluxes and their derivatives with respect
   to ground temperature using ground temperatures from previous time step.

   Note:
   the "displa" displacement height variable used here is local to these functions.
   It is not the same as the displa calculated in CanopyTemperature
*/

#include "BareGroundFluxes.h"
#include "FrictionVelocity.h"
#include "QSat.h"

namespace ELM {

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
void BareGroundFluxes::InitializeFlux(const LandType &Land, const int &frac_veg_nosno_in, const double &forc_u,
                                      const double &forc_v, const double &forc_q_in, const double &forc_th_in,
                                      const double &forc_hgt_u_patch_in, const double &thm_in, const double &thv_in,
                                      const double &t_grnd, const double &qg, const double &z0mg_in, double &dlrad,
                                      double &ulrad) {
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
    dlrad = 0.0;
    ulrad = 0.0;

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
void BareGroundFluxes::StabilityIteration(const LandType &Land, const double &forc_hgt_t_patch,
                                          const double &forc_hgt_q_patch, const double &z0mg, const double &zii,
                                          const double &beta, double &z0hg, double &z0qg) {

  double tstar;                // temperature scaling parameter
  double qstar;                // moisture scaling parameter
  double thvstar;              // virtual potential temperature scaling parameter
  double wc;                   // convective velocity [m/s]
  double zeta;                 // dimensionless height used in Monin-Obukhov theory
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
      zeta =
          zldis * vkc * grav * thvstar / (pow(ustar, 2.0) * thv); // dimensionless height used in Monin-Obukhov theory

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

} // namespace ELM
