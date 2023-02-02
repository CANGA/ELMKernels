/* functions from FrictionVelocityMod.F90
friction_velocity:
Split into 5 functions - wind, temperature, humidity, 2-m temp, 2-m humidity
Don't need u10 variables, they are used for dust model

The friction velocity scheme is based on the work of Zeng et al. (1998):
Intercomparison of bulk aerodynamic algorithms for the computation
of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
Vol. 11, 2628-2644.
*/
#pragma once

namespace ELM::friction_velocity {

/* StabilityFunc1() - used in friction_velocity_wind() */
ACCELERATE
double StabilityFunc1(const double& zeta)
{
  const double chik2{std::sqrt(1.0 - 16.0 * zeta)};
  double chik = std::sqrt(chik2);
  double retval = 2.0 * std::log((1.0 + chik) * 0.5) + std::log((1.0 + chik2) * 0.5) - 2.0 * atan(chik) + ELMconst::ELM_PI() * 0.5;
  return retval;
}

/* StabilityFunc2() - used in friction_velocity_temp() and friction_velocity_humidity() */
ACCELERATE
double StabilityFunc2(const double& zeta)
{
  const double chik2{std::sqrt(1.0 - 16.0 * zeta)};
  double retval = 2.0 * std::log((1.0 + chik2) * 0.5);
  return retval;
}

ACCELERATE
void monin_obukhov_length(const double& ur, const double& thv, const double& dthv, const double& zldis,
                          const double& z0m, double& um, double& obu)
{
  static constexpr double wc{0.5}; // convective velocity [m/s]

  // wind speed
  if (dthv >= 0.0) {
    um = std::max(ur, 0.1);
  } else {
    um = std::sqrt(ur * ur + wc * wc);
  }

  const double rib{ELMconst::GRAV() * zldis * dthv / (thv * um * um)}; // bulk Richardson number
  double zeta; // dimensionless height used in Monin-Obukhov theory
  if (rib >= 0.0) { // neutral or stable
    zeta = rib * std::log(zldis / z0m) / (1.0 - 5.0 * std::min(rib, 0.19));
    zeta = std::min(2.0, std::max(zeta, 0.01));
  } else { // unstable
    zeta = rib * std::log(zldis / z0m);
    zeta = std::max(-100.0, std::min(zeta, -0.01));
  }

  // monin-obukhov length
  obu = zldis / zeta;
}

ACCELERATE
void friction_velocity_wind(const double& forc_hgt_u_patch, const double& displa, const double& um,
                                const double& obu, const double& z0m, double& ustar)
{
  using ELMconst::VKC;
  static constexpr double zetam{1.574};          // transition point of flux-gradient relation (wind profile)
  const double zldis{forc_hgt_u_patch - displa}; // reference height "minus" zero displacement heght [m]
  const double zeta{zldis / obu};                // dimensionless height used in Monin-Obukhov theory

  if (zeta < (-zetam)) {
    ustar = VKC() * um /
            (std::log(-zetam * obu / z0m) - StabilityFunc1(-zetam) + StabilityFunc1(z0m / obu) +
             1.14 * (pow((-zeta), 0.333) - pow(zetam, 0.333)));
  } else if (zeta < 0.0) {
    ustar = VKC() * um / (std::log(zldis / z0m) - StabilityFunc1(zeta) + StabilityFunc1(z0m / obu));
  } else if (zeta <= 1.0) {
    ustar = VKC() * um / (std::log(zldis / z0m) + 5.0 * zeta - 5.0 * z0m / obu);
  } else {
    ustar = VKC() * um / (std::log(obu / z0m) + 5.0 - 5.0 * z0m / obu + (5.0 * std::log(zeta) + zeta - 1.0));
  }
}

ACCELERATE
void friction_velocity_temp(const double& forc_hgt_t_patch, const double& displa, const double& obu,
                            const double& z0h, double& temp1)
{
  using ELMconst::VKC;
  static constexpr double zetat{0.465};                     // transition point of flux-gradient relation (temp. profile)
  const double zldis{forc_hgt_t_patch - displa}; // reference height "minus" zero displacement heght [m]
  const double zeta{zldis / obu};                // dimensionless height used in Monin-Obukhov theory

  if (zeta < (-zetat)) {
    temp1 = VKC() / (std::log(-zetat * obu / z0h) - StabilityFunc2(-zetat) + StabilityFunc2(z0h / obu) +
                   0.8 * (pow(zetat, -0.333) - pow((-zeta), -0.333)));
  } else if (zeta < 0.0) {
    temp1 = VKC() / (std::log(zldis / z0h) - StabilityFunc2(zeta) + StabilityFunc2(z0h / obu));
  } else if (zeta <= 1.0) {
    temp1 = VKC() / (std::log(zldis / z0h) + 5.0 * zeta - 5.0 * z0h / obu);
  } else {
    temp1 = VKC() / (std::log(obu / z0h) + 5.0 - 5.0 * z0h / obu + (5.0 * std::log(zeta) + zeta - 1.0));
  }
}

ACCELERATE
void friction_velocity_humidity(const double& forc_hgt_q_patch, const double& forc_hgt_t_patch,
                                    const double& displa, const double& obu, const double& z0h, const double& z0q,
                                    const double& temp1, double& temp2)
{
  using ELMconst::VKC;
  static constexpr double zetat{0.465}; // transition point of flux-gradient relation (temp. profile)

  if (forc_hgt_q_patch == forc_hgt_t_patch && z0q == z0h) {
    temp2 = temp1;
  } else {
    double zldis = forc_hgt_q_patch - displa; // reference height "minus" zero displacement heght [m]
    double zeta = zldis / obu;                // dimensionless height used in Monin-Obukhov theory
    if (zeta < (-zetat)) {
      temp2 = VKC() / (std::log(-zetat * obu / z0q) - StabilityFunc2(-zetat) + StabilityFunc2(z0q / obu) +
                     0.8 * (pow(zetat, -0.333) - pow((-zeta), -0.333)));
    } else if (zeta < 0.0) {
      temp2 = VKC() / (std::log(zldis / z0q) - StabilityFunc2(zeta) + StabilityFunc2(z0q / obu));
    } else if (zeta <= 1.0) {
      temp2 = VKC() / (std::log(zldis / z0q) + 5.0 * zeta - 5. * z0q / obu);
    } else {
      temp2 = VKC() / (std::log(obu / z0q) + 5.0 - 5.0 * z0q / obu + (5.0 * std::log(zeta) + zeta - 1.0));
    }
  }
}

ACCELERATE
void friction_velocity_temp2m(const double& obu, const double& z0h, double& temp12m)
{
  using ELMconst::VKC;
  const double zldis{2.0 + z0h};
  const double zeta{zldis / obu};
  const double zetat{0.465}; // transition point of flux-gradient relation (temp. profile)

  if (zeta < -zetat) {
    temp12m = VKC() / (std::log(-zetat * obu / z0h) - StabilityFunc2(-zetat) + StabilityFunc2(z0h / obu) +
                     0.8 * (pow(zetat, -0.333) - pow(-zeta, -0.333)));
  } else if (zeta < 0.0) {
    temp12m = VKC() / (std::log(zldis / z0h) - StabilityFunc2(zeta) + StabilityFunc2(z0h / obu));
  } else if (zeta <= 1.0) {
    temp12m = VKC() / (std::log(zldis / z0h) + 5.0 * zeta - 5.0 * z0h / obu);
  } else {
    temp12m = VKC() / (std::log(obu / z0h) + 5.0 - 5.0 * (z0h / obu) + (5.0 * std::log(zeta) + zeta - 1.0));
  }
}

ACCELERATE
void friction_velocity_humidity2m(const double& obu, const double& z0h, const double& z0q, const double& temp12m,
                                      double& temp22m)
{
  using ELMconst::VKC;
  
  if (z0q == z0h) {
    temp22m = temp12m;
  } else {
    const double zldis{2.0 + z0q};
    const double zeta{zldis / obu};
    const double zetat{0.465}; // transition point of flux-gradient relation (temp. profile)
    if (zeta < -zetat) {
      temp22m = VKC() / (std::log(-zetat * obu / z0q) - StabilityFunc2(-zetat) + StabilityFunc2(z0q / obu) +
                       0.8 * (pow(zetat, -0.333) - pow(-zeta, -0.333)));
    } else if (zeta < 0.0) {
      temp22m = VKC() / (std::log(zldis / z0q) - StabilityFunc2(zeta) + StabilityFunc2(z0q / obu));
    } else if (zeta <= 1.0) {
      temp22m = VKC() / (std::log(zldis / z0q) + 5.0 * zeta - 5.0 * z0q / obu);
    } else {
      temp22m = VKC() / (std::log(obu / z0q) + 5.0 - 5.0 * z0q / obu + (5.0 * std::log(zeta) + zeta - 1.0));
    }
  }
}

} // namespace ELM::friction_velocity
