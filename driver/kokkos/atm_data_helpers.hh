
#pragma once

#include "compile_options.hh"
#include "data_types.hh"
#include "atm_data.h"

namespace ELM {

  void read_forcing(std::shared_ptr<ELM::AtmForcObjects<ViewD1, ViewD2>> atm_forcing,
                    const Utils::DomainDecomposition<2>& dd,
                    const ELM::Utils::Date& current,
                    const int atm_nsteps);

  void get_forcing(std::shared_ptr<ELM::AtmForcObjects<ViewD1, ViewD2>> atm_forcing,
                   ELMStateType& S,
                   const double& model_dt, const ELM::Utils::Date& time_plus_half_dt);

} // namespace ELM
