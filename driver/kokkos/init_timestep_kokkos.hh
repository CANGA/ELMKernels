
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

void kokkos_init_timestep(ELMStateType& S,
                          const double dtime,
                          const Utils::Date& current,
                          const std::string& fname_surfdata);

}
