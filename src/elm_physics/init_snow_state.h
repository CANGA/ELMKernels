
#pragma once

#include "array.hh"
#include "elm_constants.h"
#include <algorithm>
#include <cmath>

#include "kokkos_includes.hh"

namespace ELM {

// from ColumnDataType.F90 and WaterStateType.F90
//-----------------------------------------------------------------------
// set cold-start initial values for select members of col_ws
//-----------------------------------------------------------------------
template <typename ArrayD1>
ACCELERATE
void init_snow_state(const bool& urbpoi, const int& snl, double& h2osno, double& int_snow, double& snow_depth,
                     double& h2osfc, double& h2ocan, double& frac_h2osfc, double& fwet, double& fdry, double& frac_sno,
                     ArrayD1 snw_rds);

// derived from initVerticalMod.F90
template <typename ArrayD1>
ACCELERATE
void init_snow_layers(const double& snow_depth, const bool& lakpoi, int& snl, ArrayD1 dz, ArrayD1 z, ArrayD1 zi);

} // namespace ELM

#include "init_snow_state_impl.hh"
