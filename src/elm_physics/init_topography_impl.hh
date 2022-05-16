
#pragma once

namespace ELM {

ACCELERATE
void init_topo_slope(double& topo_slope)
{
  // check for near zero slopes, set minimum value
  topo_slope = std::max(topo_slope, 0.2);
}

ACCELERATE
void init_micro_topo(const int& ltype, const double& topo_slope, const double& topo_std, double& n_melt,
                    double& micro_sigma)
{
  if (ltype == LND::istice_mec) {
    /* ice_mec columns already account for subgrid topographic variability through
    their use of multiple elevation classes; thus, to avoid double-accounting for
    topographic variability in these columns, we ignore topo_std and use a value
    of n_melt that assumes little topographic variability within the column */
    n_melt = 10.0;
  } else {
    n_melt = 200.0 / std::max(10.0, topo_std);
  }
  // microtopographic parameter, units are meters (try smooth function of slope)
  double slopebeta = 3.0;
  double slopemax = 0.4;
  double slope0 = pow(slopemax, (-1.0 / slopebeta));
  micro_sigma = pow((topo_slope + slope0), -slopebeta);
}

} // namespace ELM
