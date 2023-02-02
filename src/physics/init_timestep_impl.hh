// clm_drv_init() from clm_driver.F90

#pragma once

namespace ELM {

template <class ArrayD1>
ACCELERATE
void init_timestep(const bool& lakpoi, const bool& veg_active,
                   const int& frac_veg_nosno_alb, const int& snl,
                   const double& h2osno, const ArrayD1 h2osoi_ice,
                   const ArrayD1 h2osoi_liq, bool& do_capsnow,
                   int& frac_veg_nosno, ArrayD1 frac_iceold)
{
  using ELMdims::nlevsno;
  
  // Decide whether to cap snow
  if (h2osno > ELMconst::H2OSNO_MAX()) {
    do_capsnow = true;
  } else {
    do_capsnow = false;
  }

  if (veg_active) {
    frac_veg_nosno = frac_veg_nosno_alb;
  } else {
    frac_veg_nosno = 0;
  }

  // Initialize set of previous time-step variables
  // Ice fraction of snow at previous time step
  if (!lakpoi) {
    for (int i = 0; i < nlevsno(); i++) {
      if (i >= nlevsno() - snl) {
        frac_iceold[i] = h2osoi_ice[i] / (h2osoi_liq[i] + h2osoi_ice[i]);
      }
    }
  }
}

} // namespace ELM
