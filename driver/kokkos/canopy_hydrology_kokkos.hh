
#pragma once

#include "elm_state.h"
#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {
  void kokkos_canopy_hydrology(const std::shared_ptr<ELMStateType>& S,
                               const double& dtime);
}
