
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {
  void kokkos_canopy_hydrology(ELMStateType& S,
                              const double& model_dt_secs,
                              const Utils::Date& time_plus_half_dt_secs);

  void kokkos_frac_wet(ELMStateType& S);
}
