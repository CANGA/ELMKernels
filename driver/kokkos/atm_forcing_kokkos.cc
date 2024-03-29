
#include "invoke_kernel.hh"
#include "atm_forcing_kokkos.hh"

namespace {
  template<ELM::AtmForcType ftype>
  using atm_forc_util = ELM::AtmDataManager<ViewD1, ViewD2, ftype>;

  template <ELM::AtmForcType ftype>
  atm_forc_util<ftype> create_forc_util(const std::string& filename,
                                      const ELM::Utils::Date &file_start_time,
                                      const int ntimes, const int ncells)
  { return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }

  template <ELM::AtmForcType ftype>
  void read_atm_data(std::shared_ptr<ELM::AtmDataManager<ViewD1, ViewD2, ftype>> atm_data,
                   const ELM::Utils::DomainDecomposition<2>& dd,
                   const ELM::Utils::Date& model_time,
                   const size_t& ntimes)
  {
    auto fptr = atm_data.get();
    auto h_data = NS::create_mirror_view(atm_data->data);
    atm_data->read_atm_forcing(h_data, dd, model_time, ntimes);
    if (atm_data->data.extent(0) != h_data.extent(0) || atm_data->data.extent(1) != h_data.extent(1))
      NS::resize(atm_data->data, h_data.extent(0), h_data.extent(1));
    NS::deep_copy(atm_data->data, h_data);
  }

} // anonymous namespace


void ELM::read_forcing(ELMStateType& S, const ELM::Utils::Date& current)
{
  auto fptr = S.atm_forcing.get();
  read_atm_data(fptr->forc_TBOT, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_PBOT, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_QBOT, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_FLDS, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_FSDS, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_PREC, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_WIND, S.dd, current, S.atm_nsteps);
  read_atm_data(fptr->forc_ZBOT, S.dd, current, S.atm_nsteps);
}


void ELM::get_forcing(ELMStateType& S,
                      const double& model_dt,
                      const ELM::Utils::Date& time_plus_half_dt)
{
  auto fptr = S.atm_forcing.get();
  fptr->forc_TBOT->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_tbot, S.forc_thbot);
  fptr->forc_PBOT->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_pbot);
  fptr->forc_QBOT->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_tbot, S.forc_pbot, S.forc_qbot);
  fptr->forc_FLDS->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_pbot, S.forc_qbot, S.forc_tbot, S.forc_lwrad);
  fptr->forc_FSDS->get_atm_forcing(model_dt, time_plus_half_dt, S.coszen, S.forc_solai, S.forc_solad);
  fptr->forc_PREC->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_tbot, S.forc_rain, S.forc_snow);
  fptr->forc_WIND->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_u, S.forc_v);
  fptr->forc_ZBOT->get_atm_forcing(model_dt, time_plus_half_dt, S.forc_hgt, S.forc_hgt_u_patch,
                                       S.forc_hgt_t_patch,S.forc_hgt_q_patch);

  // calculate constitutive air properties
//  ELM::atm_forcing_physics::ConstitutiveAirProperties
//    compute_air_props(
//      S.forc_qbot, S.forc_pbot,
//      S.forc_tbot, S.forc_vp,
//      S.forc_rho, S.forc_po2,
//      S.forc_pco2);
//  apply_parallel_for(compute_air_props, "ConstitutiveAirProperties", S.forc_pbot.extent(0));

}

