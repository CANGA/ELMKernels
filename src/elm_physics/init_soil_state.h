
#pragma once

#include "array.hh"
#include "elm_constants.h"
#include "land_data.h"

#include <algorithm>
#include <assert.h>
#include <cmath>

#include "kokkos_includes.hh"

namespace ELM::init_soil_state {

// from ColumnDataType.F90
//-----------------------------------------------------------------------
// set cold-start initial values for select members of col_es
//-----------------------------------------------------------------------
template <typename ArrayD1>
ACCELERATED
void init_soil_temp(const LandType& Land, const int& snl, ArrayD1 t_soisno, double& t_grnd);

// from ColumnDataType.F90 and WaterStateType.F90
//--------------------------------------------
// Set soil water
//--------------------------------------------
// volumetric water is set first and liquid content and ice lens are obtained
// NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
// and urban pervious road (other urban columns have zero soil water)
template <typename ArrayD1>
ACCELERATED
void init_soilh2o_state(const LandType& Land, const int& snl, const ArrayD1& watsat, const ArrayD1& t_soisno,
                        const ArrayD1& dz, ArrayD1 h2osoi_vol, ArrayD1 h2osoi_liq, ArrayD1 h2osoi_ice);

/*
DESCRIPTION:
compute root profile for soil water uptake using equation from Zeng 2001, J. Hydrometeorology

INPUTS:
vtype [int] vegetation type
roota_par [double] CLM rooting distribution parameter [1/m]
rootb_par [double] CLM rooting distribution parameter [1/m]

OUTPUT:
rootfr_road_perv[nlevgrnd]   [double] fraction of roots in each soil layer
*/
template <typename ArrayD1>
ACCELERATED
void init_vegrootfr(const int& vtype, const double& roota_par, const double& rootb_par, const ArrayD1& zi,
                    ArrayD1 rootfr);

} // namespace ELM::init_soil_state

#include "init_soil_state_impl.hh"
