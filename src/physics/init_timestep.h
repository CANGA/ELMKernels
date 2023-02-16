/*! \file init_timestep.hh
\brief Function derived from clm_driver.F90
*/

#pragma once

#include "array.hh"
#include "elm_constants.h"

#include "compile_options.hh"

namespace ELM {

/*! Initialize variables at the beginning of each time cycle.

\param[in]  lakpoi                       [bool] true if lake point
\param[in]  veg_active                   [bool] true if active vegetation present
\param[in]  frac_veg_nosno_alb           [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  snl                          [int] number of snow layers
\param[in]  h2osno                       [double] snow water (mm H2O)
\param[in]  h2osoi_ice[nlevgrnd()+nlevsno()] [double] ice lens (kg/m2)
\param[in]  h2osoi_liq[nlevgrnd()+nlevsno()] [double] liquid water (kg/m2)
\param[out] do_capsnow                   [int] true => do snow capping
\param[out] frac_veg_nosno               [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[out] frac_iceold                  [double] fraction of ice relative to the tot water
*/
template <typename ArrayD1>
ACCELERATE
void init_timestep(const bool& lakpoi, const bool& veg_active,
                   const int& frac_veg_nosno_alb, const int& snl,
                   const double& h2osno, const ArrayD1 h2osoi_ice,
                   const ArrayD1 h2osoi_liq, int& do_capsnow,
                   int& frac_veg_nosno, ArrayD1 frac_iceold);

} // namespace ELM

#include "init_timestep_impl.hh"
