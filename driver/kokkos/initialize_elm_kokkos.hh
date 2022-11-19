

#pragma once


#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

void initialize_kokkos_elm(
  ELMStateType& S,
  const std::string& fname_surfdata,
  const std::string& fname_param,
  const std::string& fname_snicar,
  const std::string& fname_snowage,
  const std::string& fname_aerosol);

}
