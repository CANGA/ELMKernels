/* clm_drv_init() from clm_driver.F90
INPUTS:
lakpoi [bool] true if lake point
h2osno [double] snow water (mm H2O)
veg_active [bool] true if active vegetation present
snl [double] number of snow layers
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
frac_veg_nosno_alb [int] fraction of vegetation not covered by snow (0 OR 1) [-]


OUTPUTS:
h2osno_old [double] snow water (mm H2O) at previous time step
do_capsnow [bool] true => do snow capping
eflx_bot [double] heat flux from beneath soil/ice column (W/m**2)
qflx_glcice [double] flux of new glacier ice (mm H2O/s) [+ = ice grows]
frac_veg_nosno [int] fraction of vegetation not covered by snow (0 OR 1) [-]
frac_iceold [double] fraction of ice relative to the tot water
*/

#pragma once

#include "clm_constants.h"

namespace ELM {

template <class dArray_type>
void InitTimestep(const bool &lakpoi, const double &h2osno, const bool &veg_active, const int &snl,
                  const dArray_type h2osoi_ice, const dArray_type h2osoi_liq, const int &frac_veg_nosno_alb,
                  double &h2osno_old, bool &do_capsnow, double &eflx_bot, double &qflx_glcice, int &frac_veg_nosno,
                  dArray_type frac_iceold) {
  // Save snow mass at previous time step
  h2osno_old = h2osno;

  // Decide whether to cap snow
  if (h2osno > h2osno_max) {
    do_capsnow = true;
  } else {
    do_capsnow = false;
  }

  // Reset flux from beneath soil/ice column
  eflx_bot = 0.0;

  // Initialize qflx_glcice everywhere, to zero.
  qflx_glcice = 0.0;

  if (veg_active) {
    frac_veg_nosno = frac_veg_nosno_alb;
  } else {
    frac_veg_nosno = 0.0;
  }

  // Initialize set of previous time-step variables
  // Ice fraction of snow at previous time step
  if (!lakpoi) {
    for (int i = 0; i < nlevsno; i++) {
      if (i <= nlevsno - snl) {
        frac_iceold[i] = h2osoi_ice[i] / (h2osoi_liq[i] + h2osoi_ice[i]);
      }
    }
  }
}

} // namespace ELM
