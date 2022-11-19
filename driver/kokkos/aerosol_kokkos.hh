
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

// convenience function to invoke aerosol deposition source functor
void invoke_aerosol_source(ELMStateType& S,
                           const double& dtime,
                           const Utils::Date& model_time);
// convenience function to invoke aerosol mass and concen functor
void invoke_aerosol_concen_and_mass(ELMStateType& S,
                                    const double& dtime);

} // namespace ELM
