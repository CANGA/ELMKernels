
#pragma once

#include "array.hh"
#include "elm_constants.h"
#include <cmath>
#include <algorithm>

#include "kokkos_includes.hh"

namespace ELM {

ACCELERATE
void pedotransfer(const double& pct_sand, const double& pct_clay, double& watsat, double& bsw, double& sucsat,
                  double& xksat);

ACCELERATE
void soil_hydraulic_params(const double& pct_sand, const double& pct_clay, const double& zsoi, const double& om_frac,
                           double& watsat, double& bsw, double& sucsat, double& watdry, double& watopt, double& watfc,
                           double& tkmg, double& tkdry, double& csol);

template <typename ArrayD1>
ACCELERATE
void init_soil_hydraulics(const ArrayD1 pct_sand, const ArrayD1 pct_clay, const ArrayD1 organic, const ArrayD1 zsoi,
                          ArrayD1 watsat, ArrayD1 bsw, ArrayD1 sucsat, ArrayD1 watdry, ArrayD1 watopt, ArrayD1 watfc,
                          ArrayD1 tkmg, ArrayD1 tkdry, ArrayD1 csol);

} // namespace ELM

#include "soil_texture_hydraulic_model_impl.hh"
