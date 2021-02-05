#pragma once

#include "clm_constants.h"
#include "vegproperties.h"

namespace ELM {

void Photosynthesis(const VegProperties &veg, const int &vtype, const int &nrad, const double &forc_pbot,
                    const double &t_veg, const double &t10, const double &esat_tv, const double &eair,
                    const double &oair, const double &cair, const double &rb, const double &btran,
                    const double &dayl_factor, const double thm, const double tlai_z[nlevcan], const double &vcmaxcint,
                    const double par_z[nlevcan], const double lai_z[nlevcan], double &rs);

} // namespace ELM
