// functions from SurfaceResistanceMod.F90

#include "elm_constants.h"
#include "LandType.h"
#include <algorithm>
#include <cmath>

namespace ELM {
namespace surface_resistance {

double getlblcef(const double &rho, const double &temp) {
  const double C = 120.0;      // K
  const double T0 = 291.25;    // K
  const double mu0 = 18.27e-6; // Pa s
  const double prandtl = 0.72;
  // compute the kinetic viscosity
  double mu = mu0 * (T0 + C) / (temp + C) * pow(temp / T0, 1.5) / rho; // m^2 s^-1
  double diffh2o = 0.229e-4 * pow(temp / 273.15, 1.75);                // m^2 s^-1
  double sc = mu / diffh2o;                                            // schmidt number
  double result = 2.0 / vkc * pow(sc / prandtl, 2.0 / 3.0);
  return result;
}

} // namespace surface_resistance
} // namespace ELM
