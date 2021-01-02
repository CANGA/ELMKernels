/* functions from FrictionVelocityMod.F90
FrictionVelocity:
Split into 3 functions - wind, temperature, humidity
!!DON'T NEED 2-meter variables, they are only used as feedbacks to ATM model and don't affect LSM and subsurface calculations!
Also don't need u10 variables, they are used for dust model

The friction velocity scheme is based on the work of Zeng et al. (1998):
Intercomparison of bulk aerodynamic algorithms for the computation
of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
Vol. 11, 2628-2644.
*/

#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

/* MoninObukIni()
DESCRIPTION: Initialization of the Monin-Obukhov length. The scheme is based on the work of 
Zeng et al. (1998): Intercomparison of bulk aerodynamic algorithms for the computation of sea 
surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11, 2628-2644.

INPUTS:
ur        [double]  wind speed at reference height [m/s]
thv       [double]  virtual potential temperature (kelvin)
dthv      [double]  diff of vir. poten. temp. between ref. height and surface
zldis     [double]  reference height "minus" zero displacement height [m]
z0m       [double]  roughness length, momentum [m]

OUTPUTS:
um       [double]  wind speed including the stability effect [m/s]
obu      [double]  monin-obukhov length (m)

*/
void MoninObukIni(
  const double& ur,
  const double& thv,
  const double& dthv,
  const double& zldis,
  const double& z0m,
  double& um,
  double& obu)
{
  double rib;          // bulk Richardson number
  double zeta;         // dimensionless height used in Monin-Obukhov theory
  // double ustar = 0.06; // friction velocity [m/s] -- in CLM, but not used
  double wc = 0.5;     // convective velocity [m/s]

  // wind speed
  if (dthv >= 0.0) {
    um = std::max(ur, 0.1);
  } else {
    um = std::sqrt(ur * ur + wc * wc);
  }

  rib = grav * zldis * dthv / (thv * um * um);

  if (rib >= 0.0) {  // neutral or stable
    zeta = rib * std::log(zldis / z0m) / (1.0 - 5.0 * std::min(rib, 0.19));
    zeta = std::min(2.0, std::max(zeta, 0.01));
  } else { // unstable
    zeta = rib * std::log(zldis / z0m);
    zeta = std::max(-100.0, std::min(zeta, -0.01));
  }

  // monin-obukhov length
  obu = zldis / zeta;
}


/* StabilityFunc1() - used in FrictionVelocityWind() */
double StabilityFunc1(double zeta)
{
  double chik, chik2, retval;

  chik2 = std::sqrt(1.0 - 16.0 * zeta);
  chik = std::sqrt(chik2);
  retval = 2.0 * std::log((1.0 + chik) * 0.5) + std::log((1.0 + chik2) * 0.5) - 2.0 * atan(chik) + SHR_CONST_PI * 0.5;
  return retval;
}

/* StabilityFunc2() - used in FrictionVelocityTemperature() and FrictionVelocityHumidity() */
double StabilityFunc2(double zeta)
{
  double chik2, retval;

  chik2 = std::sqrt(1.0 - 16.0 * zeta);
  retval = 2.0 * std::log((1.0 + chik2) * 0.5);
  return retval;
}


/* FrictionVelocityWind()
DESCRIPTION:
Calculation of the friction velocity of surface boundary layer.
INPUTS:
forc_hgt_u_patch [double] observational height of wind at pft level [m]
displa           [double] displacement height (m)
um               [double] wind speed including the stability effect [m/s]
obu              [double] monin-obukhov length (m)
z0m              [double] roughness length over vegetation, momentum [m]

OUTPUTS:
ustar            [double] friction velocity [m/s]
*/
void FrictionVelocityWind(
  const double& forc_hgt_u_patch,
  const double& displa,
  const double& um,
  const double& obu,
  const double& z0m,

  double& ustar)
{
  double zetam = 1.574; // transition point of flux-gradient relation (wind profile)
  double zldis = forc_hgt_u_patch - displa; // reference height "minus" zero displacement heght [m]
  double zeta = zldis / obu; // dimensionless height used in Monin-Obukhov theory

  if (zeta < (-zetam)) {
    ustar = vkc * um / (std::log(-zetam * obu / z0m) - StabilityFunc1(-zetam) + 
      StabilityFunc1(z0m / obu) + 1.14 * (pow((-zeta), 0.333) - pow(zetam, 0.333)));
  } else if (zeta < 0.0) {
    ustar = vkc * um / (std::log(zldis / z0m) - StabilityFunc1(zeta) + StabilityFunc1(z0m / obu));
  } else if (zeta <=  1.0) {
    ustar = vkc * um / (std::log(zldis / z0m) + 5.0 * zeta - 5.0 * z0m / obu);
  } else {
    ustar = vkc * um / (std::log(obu / z0m) + 5.0 - 5.0 * z0m / obu + (5.0 * std::log(zeta) + zeta - 1.0));
  }
}

/* FrictionVelocityTemperature()
DESCRIPTION:
Calculation of the relation for potential temperature of surface boundary layer.
INPUTS:
forc_hgt_t_patch [double] observational height of temperature at pft level [m]
displa           [double] displacement height (m)
obu              [double] monin-obukhov length (m)
z0h              [double] roughness length over vegetation, sensible heat [m]


OUTPUTS:
temp1            [double] relation for potential temperature profile

*/
void FrictionVelocityTemperature(
  const double& forc_hgt_t_patch,
  const double& displa,
  const double& obu,
  const double& z0h,

  double& temp1)
{
  double zetat = 0.465; // transition point of flux-gradient relation (temp. profile)
  double zldis = forc_hgt_t_patch - displa; // reference height "minus" zero displacement heght [m]
  double zeta = zldis / obu; // dimensionless height used in Monin-Obukhov theory

  if (zeta < (-zetat)) {
    temp1 = vkc / (std::log(-zetat * obu / z0h) - StabilityFunc2(-zetat) + 
      StabilityFunc2(z0h / obu) + 0.8 * (pow(zetat, -0.333) - pow((-zeta), -0.333)));
  } else if (zeta < 0.0) {
    temp1 = vkc / (std::log(zldis / z0h) - StabilityFunc2(zeta) + StabilityFunc2(z0h / obu));
  } else if (zeta <=  1.0) {
    temp1 = vkc / (std::log(zldis / z0h) + 5.0 * zeta - 5.0 * z0h / obu);
  } else {
    temp1 = vkc / (std::log(obu / z0h) + 5.0 - 5.0 * z0h / obu + (5.0 * std::log(zeta) + zeta - 1.0));
  }
}


/* FrictionVelocityHumidity()
DESCRIPTION:
Calculation of the relation for potential humidity of surface boundary layer.
INPUTS:
forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
forc_hgt_t_patch [double] observational height of temperature at pft level 
displa           [double] displacement height (m)
obu              [double] monin-obukhov length (m)
z0h              [double] roughness length over vegetation, sensible heat [m]
z0q              [double] roughness length over vegetation, latent heat [m]
temp1            [double] relation for potential temperature profile

OUTPUTS:
temp2            [double] relation for specific humidity profile
*/
void FrictionVelocityHumidity(
  const double& forc_hgt_q_patch,
  const double& forc_hgt_t_patch,
  const double& displa,
  const double& obu,
  const double& z0h,
  const double& z0q,
  const double& temp1,

  double& temp2)
{
  double zetat = 0.465; // transition point of flux-gradient relation (temp. profile)

  if (forc_hgt_q_patch == forc_hgt_t_patch && z0q == z0h) {
    temp2 = temp1;
  } else {
    double zldis = forc_hgt_q_patch - displa; // reference height "minus" zero displacement heght [m]
    double zeta = zldis / obu; // dimensionless height used in Monin-Obukhov theory
    if (zeta < (-zetat)) {
      temp2 = vkc / (std::log(-zetat * obu / z0q) - StabilityFunc2(-zetat) + 
        StabilityFunc2(z0q / obu) + 0.8 * (pow(zetat, -0.333) - pow((-zeta), -0.333)));
    } else if (zeta < 0.0) {
      temp2 = vkc / (std::log(zldis / z0q) - StabilityFunc2(zeta) + StabilityFunc2(z0q / obu));
    } else if (zeta <=  1.0) {
      temp2 = vkc / (std::log(zldis / z0q) + 5.0 * zeta - 5. * z0q / obu);
    } else {
      temp2 = vkc / (std::log(obu / z0q) + 5.0 - 5.0 * z0q / obu + (5.0 * std::log(zeta) + zeta - 1.0));
    }
  }
}

