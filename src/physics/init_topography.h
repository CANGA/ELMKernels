/*! \file init_topography.h
\brief Functions derived from from initVerticalMod.F90
*/
#pragma once

#include "land_data.h"
#include "init_topography.h"
#include "elm_constants.h"
#include <algorithm>
#include <cmath>

#include "compile_options.hh"

namespace ELM {

/*! Set minimum slope value
/param[in]  raw_topo_slope  [double] raw_topo_slope (-)
\param[out] topo_slope      [double] processed slope (-)
*/
ACCELERATE
double init_topo_slope(double raw_topo_slope);

/*! Initialize n_melt for SCA calculations.
\param[in]  topo_std    [double] gridcell elevation standard deviation
\param[out] n_melt      [double] SCA shape parameter
*/
ACCELERATE
double init_melt_factor(int ltype, double topo_std);

/*! Initialize micro_sigma for SCA calculations.
\param[in]  topo_slope  [double] gridcell topographic slope -- this must be % slope
\param[out] micro_sigma [double] microtopography pdf sigma (m)
*/
ACCELERATE
double init_micro_sigma(double topo_slope);

} // namespace ELM

#include "init_topography_impl.hh"
