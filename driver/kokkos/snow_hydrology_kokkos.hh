
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

  void kokkos_snow_hydrology(ELMStateType& S,
                             const double& dtime,
                             const ELM::Utils::Date& time_plus_half_dt);

}
