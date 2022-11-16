
#include <string>
#include <unordered_map>

#include "invoke_kernel.hh"
#include "phenology_kokkos.hh"

// these phenology data functions will stay here until I decide what to do about host-device transfer
std::unordered_map<std::string, h_ViewD2>
ELM::get_phen_host_views(const ELM::PhenologyDataManager<ViewD2> *phen_data)
{
  std::unordered_map<std::string, h_ViewD2> phen_host_views;
  phen_host_views["MONTHLY_LAI"] = Kokkos::create_mirror_view(phen_data->mlai);
  phen_host_views["MONTHLY_SAI"] = Kokkos::create_mirror_view(phen_data->msai);
  phen_host_views["MONTHLY_HEIGHT_TOP"] = Kokkos::create_mirror_view(phen_data->mhtop);
  phen_host_views["MONTHLY_HEIGHT_BOT"] = Kokkos::create_mirror_view(phen_data->mhbot);
  return phen_host_views;
}


void ELM::update_phenology(ELM::PhenologyDataManager<ViewD2> *phen_data,
  std::unordered_map<std::string, h_ViewD2>& host_phen_views,
  ELMStateType *S,
  const h_ViewI1 vtype,
  const ELM::Utils::Date& current,
  const std::string& fname_surfdata)
{
  // copy device data to host
  // copying entire views is likely inefficient, but it's currently necessary
  // could be eliminated by shifting older month indices in parallel kernel
  // and reading new data into a mirror of a subview (or a subview of a mirror?)
  // then we would only need one copy from host view into the device view
  // instead of the two we currently have
  // will fix later - too infrequently run (once per month) to cause concern
  if (phen_data->need_data()) {
    NS::deep_copy(host_phen_views["MONTHLY_LAI"], phen_data->mlai);
    NS::deep_copy(host_phen_views["MONTHLY_SAI"], phen_data->msai);
    NS::deep_copy(host_phen_views["MONTHLY_HEIGHT_TOP"], phen_data->mhtop);
    NS::deep_copy(host_phen_views["MONTHLY_HEIGHT_BOT"], phen_data->mhbot);
  }
  // reads three months of data on first call
  // after first call, read new data if phen_data.need_new_data_ == true
  auto phen_updated = phen_data->read_data(host_phen_views, fname_surfdata, current, vtype); // if needed
  // copy host views to device
  // could be made more efficient, see above
  if (phen_updated) {
    NS::deep_copy(phen_data->mlai, host_phen_views["MONTHLY_LAI"]);
    NS::deep_copy(phen_data->msai, host_phen_views["MONTHLY_SAI"]);
    NS::deep_copy(phen_data->mhtop, host_phen_views["MONTHLY_HEIGHT_TOP"]);
    NS::deep_copy(phen_data->mhbot, host_phen_views["MONTHLY_HEIGHT_BOT"]);
  }
  // run parallel kernel to process phenology data
  phen_data->get_data(current, S->snow_depth,
                     S->frac_sno, S->vtype, S->elai, S->esai,
                     S->htop, S->hbot, S->tlai, S->tsai,
                     S->frac_veg_nosno_alb);
}
