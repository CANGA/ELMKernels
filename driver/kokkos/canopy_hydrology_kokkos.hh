
#pragma once

#include "atm_data.h"
#include "elm_state.h"
#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {
  void kokkos_canopy_hydrology(ELMStateType& S, AtmDataManager<ViewD1, ViewD2, AtmForcType::PREC>& forc_PREC,
                                   const double& model_dt_secs, const Utils::Date& time_plus_half_dt_secs);

  void kokkos_frac_wet(ELMStateType& S);
}
