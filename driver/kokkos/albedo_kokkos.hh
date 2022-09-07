
#pragma once

#include "elm_state.h"
#include "aerosol_physics.h"
#include "snicar_data.h"
#include "pft_data.h"

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

void kokkos_albedo_snicar(ELMStateType& S,
                               ELM::AerosolConcentrations<ViewD2>& aerosol_concentrations,
                               ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data,
                               ELM::PFTData<ViewD1>& pft_data);
}
