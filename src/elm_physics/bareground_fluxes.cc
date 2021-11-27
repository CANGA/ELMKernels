/* functions derived from Bareground_fluxesMod.F90
   Compute sensible and latent fluxes and their derivatives with respect
   to ground temperature using ground temperatures from previous time step.

   Note:
   the "displa" displacement height variable used here is local to these functions.
   It is not the same as the displa calculated in CanopyTemperature
*/

#include "bareground_fluxes.h"
#include "friction_velocity.h"
#include "qsat.h"

namespace ELM {
namespace bareground_fluxes {

void initialize_flux(const LandType &Land, const int &frac_veg_nosno, const double &forc_u, const double &forc_v,
                       const double &forc_q, const double &forc_th, const double &forc_hgt_u_patch, const double &thm,
                       const double &thv, const double &t_grnd, const double &qg, const double &z0mg, double &dlrad,
                       double &ulrad, double &zldis, double &displa, double &dth, double &dqh, double &obu, double &ur,
                       double &um) {

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

void stability_iteration(const LandType &Land, const int &frac_veg_nosno, const double &forc_hgt_t_patch,
                           const double &forc_hgt_u_patch, const double &forc_hgt_q_patch, const double &z0mg,
                           const double &zldis, const double &displa,
                           const double &dth, const double &dqh, const double &ur, const double &forc_q,
                           const double &forc_th, const double &thv, double &z0hg, double &z0qg, double &obu,
                           double &um, double &temp1, double &temp2, double &temp12m, double &temp22m, double &ustar) {

  double tstar;                     // temperature scaling parameter
  double qstar;                     // moisture scaling parameter
  double thvstar;                   // virtual potential temperature scaling parameter
  double wc;                        // convective velocity [m/s]
  double zeta;                      // dimensionless height used in Monin-Obukhov theory
  static const int niters = 3;      // number of iterations
  static const double beta = 1.0;   // coefficient of convective velocity [-]
  static const double zii = 1000.0; //convective boundary height [m]

  for (int i = 0; i < niters; i++) {

    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno == 0) {

      // friction velocity calls
      friction_velocity::friction_velocity_wind(forc_hgt_u_patch, displa, um, obu, z0mg, ustar);
      friction_velocity::friction_velocity_temp(forc_hgt_t_patch, displa, obu, z0hg, temp1);
      friction_velocity::friction_velocity_humidity(forc_hgt_q_patch, forc_hgt_t_patch, displa, obu, z0hg, z0qg, temp1, temp2);
      friction_velocity::friction_velocity_temp2m(obu, z0hg, temp12m);
      friction_velocity::friction_velocity_humidity2m(obu, z0hg, z0qg, temp12m, temp22m);

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
} // stability_iteration()

} // namespace bareground_fluxes
} // namespace ELM
