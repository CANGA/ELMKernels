
#pragma once

#include "compile_options.hh"
#include "data_types.hh"
#include "atm_data.h"

namespace ELM {

  void read_forcing(ELMStateType& S,
                    const ELM::Utils::Date& current);

  void get_forcing(ELMStateType& S,
                   const double& model_dt,
                   const ELM::Utils::Date& time_plus_half_dt);

} // namespace ELM
