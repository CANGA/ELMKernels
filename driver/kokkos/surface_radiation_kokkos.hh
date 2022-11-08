
#pragma once

#include "atm_data.h"
#include "elm_state.h"
#include "elm_constants.h"
#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

  void kokkos_surface_radiation(ELMStateType& S, AtmDataManager<ViewD1, ViewD2, AtmForcType::FSDS>& forc_FLDS,
                                const double& model_dt_days, const Utils::Date& time_plus_half_dt);

}
