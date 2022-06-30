
#pragma once

#include "elm_state.h"
#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

void kokkos_init_timestep(const std::shared_ptr<ELMStateType>& S);

}
