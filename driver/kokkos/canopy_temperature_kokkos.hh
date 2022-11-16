
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

  void kokkos_canopy_temperature(ELMStateType& S,
                                 ELM::PFTData<ViewD1>& pft_data);

}
