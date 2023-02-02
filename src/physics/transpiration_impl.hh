


#pragma once

namespace ELM::trans {

/*
Generic routine to apply transpiration as a sink condition that
is vertically distributed over the soil column. Should be
applicable to any Richards solver that is not coupled to plant
hydraulics.

vegetation transpiration (mm H2O/s) (+ = to atm)
*/
template <typename ArrayD1>
ACCELERATE
void transpiration(const bool& veg_active, 
                   const double& qflx_tran_veg,
                   const ArrayD1 rootr,
                   ArrayD1 qflx_rootsoi)
{
  using ELMdims::nlevsoi;
  if (veg_active) {
    for (int i = 0; i < nlevsoi(); ++i) {
      qflx_rootsoi(i) = rootr(i) * qflx_tran_veg;
    }
  }
}

} // namespace ELM::trans
