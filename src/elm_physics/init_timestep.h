/*! \file init_timestep.hh
\brief Function derived from clm_driver.F90
*/

#pragma once

#include "elm_constants.h"
#include "array.hh"

namespace ELM {

/*! Initialize variables at the beginning of each time cycle.

\param[in]  lakpoi                       [bool] true if lake point
\param[in]  h2osno                       [double] snow water (mm H2O)
\param[in]  veg_active                   [bool] true if active vegetation present
\param[in]  snl                          [double] number of snow layers
\param[in]  h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
\param[in]  h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
\param[in]  frac_veg_nosno_alb           [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[out] h2osno_old                   [double] snow water (mm H2O) at previous time step
\param[out] do_capsnow                   [bool] true => do snow capping
\param[out] eflx_bot                     [double] heat flux from beneath soil/ice column (W/m**2)
\param[out] qflx_glcice                  [double] flux of new glacier ice (mm H2O/s) [+ = ice grows]
\param[out] frac_veg_nosno               [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[out] frac_iceold                  [double] fraction of ice relative to the tot water
*/
template <class ArrayD1>
void init_timestep(const bool &lakpoi, const double &h2osno, const bool &veg_active, const int &snl,
                  const ArrayD1 h2osoi_ice, const ArrayD1 h2osoi_liq, const int &frac_veg_nosno_alb, double &h2osno_old,
                  bool &do_capsnow, double &eflx_bot, double &qflx_glcice, int &frac_veg_nosno, ArrayD1 frac_iceold);

} // namespace ELM

#include "init_timestep_impl.hh"

