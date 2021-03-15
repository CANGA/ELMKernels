// clm_drv_init() from clm_driver.F90

#pragma once

#include "ELMConstants.h"

namespace ELM {

template <class ArrayD1>
void InitTimestep(const bool &lakpoi, const double &h2osno, const bool &veg_active, const int &snl,
                  const ArrayD1 h2osoi_ice, const ArrayD1 h2osoi_liq, const int &frac_veg_nosno_alb, double &h2osno_old,
                  bool &do_capsnow, double &eflx_bot, double &qflx_glcice, int &frac_veg_nosno, ArrayD1 frac_iceold) {
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
      if (i >= nlevsno - snl) {
        frac_iceold[i] = h2osoi_ice[i] / (h2osoi_liq[i] + h2osoi_ice[i]);
      }
    }
  }
}

} // namespace ELM
