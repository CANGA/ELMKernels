
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

#include "phenology_data.h"

namespace ELM {

std::unordered_map<std::string, h_ViewD2> get_phen_host_views(
  const ELM::PhenologyDataManager<ViewD2>& phen_data);

void update_phenology(
  ELMStateType& S,
  const ELM::Utils::Date& current,
  const std::string& fname_surfdata);

} // namespace ELM
