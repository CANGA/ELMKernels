/*
//from initVerticalMod
Initialize n_melt and micro_sigma for SCA calculations
INPUTS:
ltype       [integer] landunit type
topo_slope  [double] gridcell topographic slope -- this must be % slope
topo_std    [double] gridcell elevation standard deviation -- TOPOGRAPHIC STDdev (STD_ELEV) -- will need to determine how to calculate these - D4, D8, D2(1d)? Likely scale-dependent. For now I will initialze as 0.0.

OUTPUTS:
n_melt      [double] SCA shape parameter
micro_sigma [double] microtopography pdf sigma (m)
*/

#include <algorithm>
#include <cmath>
#include "clm_constants.h"

void TopoSlopes(const double& topo_slope_in,
    double& topo_slope_out)
{
    //check for near zero slopes, set minimum value
    topo_slope_out = std::max(topo_slope_in, 0.2);
}

void ColdStartMicroTopo(const int& ltype,
    const double& topo_slope,
    const double& topo_std,

    double& n_melt,
    double& micro_sigma)
{
    if (ltype == istice_mec) {
        /* ice_mec columns already account for subgrid topographic variability through
        their use of multiple elevation classes; thus, to avoid double-accounting for
        topographic variability in these columns, we ignore topo_std and use a value
        of n_melt that assumes little topographic variability within the column */
        n_melt = 10.0;
    } else {
        n_melt = 200.0 / std::max(10.0, topo_std);
    }
    //microtopographic parameter, units are meters (try smooth function of slope)
    double slopebeta = 3.0;
    double slopemax = 0.4;
    double slope0 = pow(slopemax, (-1.0 / slopebeta));
    micro_sigma = pow((topo_slope + slope0), -slopebeta);
}
