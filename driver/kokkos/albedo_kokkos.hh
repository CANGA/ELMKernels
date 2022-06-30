
#pragma once

#include "elm_state.h"
#include "aerosol_physics.h"
#include "snicar_data.h"
#include "pft_data.h"

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

void kokkos_albedo_snicar(const std::shared_ptr<ELMStateType>& S,
                          const std::shared_ptr<ELM::AerosolConcentrations<ViewD2>>& aerosol_concentrations,
                          const std::shared_ptr<ELM::SnicarData<ViewD1, ViewD2, ViewD3>>& snicar_data,
                          const std::shared_ptr<ELM::PFTData<ViewD1, ViewD2>>& pft_data);
}
