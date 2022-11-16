
#include "compile_options.hh"
#include "data_types.hh"

#include "phenology_data.h"

namespace ELM {

std::unordered_map<std::string, h_ViewD2> get_phen_host_views(
  const ELM::PhenologyDataManager<ViewD2> *phen_data);

void update_phenology(
  ELM::PhenologyDataManager<ViewD2> *phen_data,
  std::unordered_map<std::string, h_ViewD2>& host_phen_views,
  ELMStateType* S,
  const h_ViewI1 vtype,
  const ELM::Utils::Date& current,
  const std::string& fname_surfdata);

} // namespace ELM
