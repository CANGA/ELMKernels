
#pragma once

namespace ELM {

ACCELERATE
double init_topo_slope(double raw_topo_slope)
{
  // check for near zero slopes, set minimum value
  return std::max(raw_topo_slope, 0.2);
}

ACCELERATE
double init_melt_factor(int ltype, double topo_std)
{
  double n_melt;
  if (ltype == LND::istice_mec) {
    /* ice_mec columns already account for subgrid topographic variability through
    their use of multiple elevation classes; thus, to avoid double-accounting for
    topographic variability in these columns, we ignore topo_std and use a value
    of n_melt that assumes little topographic variability within the column */
    n_melt = 10.0;
  } else {
    n_melt = 200.0 / std::max(10.0, topo_std);
  }
  return n_melt;
}

ACCELERATE
double init_micro_sigma(double topo_slope)
{
  // microtopographic parameter, units are meters (try smooth function of slope)
  const double slopebeta = 3.0;
  const double slopemax = 0.4;
  const double slope0 = pow(slopemax, (-1.0 / slopebeta));
  return pow((topo_slope + slope0), -slopebeta);
}

} // namespace ELM
