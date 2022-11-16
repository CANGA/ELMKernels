
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {
  
  void kokkos_soil_temperature(ELMStateType& S,
                               const double& dtime);

}
