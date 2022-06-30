
#pragma once

#include "elm_state.h"
#include "pft_data.h"
#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

  void kokkos_canopy_temperature(const std::shared_ptr<ELMStateType>& S,
                                 const std::shared_ptr<ELM::PFTData<ViewD1, ViewD2>>& pft_data);

}
