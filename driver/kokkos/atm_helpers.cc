
#include "kokkos_includes.hh"
#include "kokkos_types.hh"
#include "elm_constants.h"
#include "atm_helpers.hh"

namespace {
  template<ELM::AtmForcType ftype>
  using atm_forc_util = ELM::AtmDataManager<ViewD1, ViewD2, ftype>;

  template <ELM::AtmForcType ftype>
  atm_forc_util<ftype> create_forc_util(const std::string& filename,
                                      const ELM::Utils::Date &file_start_time,
                                      const int ntimes, const int ncells)
  { return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }

  template <ELM::AtmForcType ftype>
  void read_atm_data(ELM::AtmDataManager<ViewD1, ViewD2, ftype>& atm_data,
                   const ELM::Utils::DomainDecomposition<2>& dd,
                   const ELM::Utils::Date& model_time,
                   const size_t& ntimes)
 {
    auto h_data = Kokkos::create_mirror_view(atm_data.data);
    atm_data.read_atm_forcing(h_data, dd, model_time, ntimes);
    if (atm_data.data.extent(0) != h_data.extent(0) || atm_data.data.extent(1) != h_data.extent(1))
      NS::resize(atm_data.data, h_data.extent(0), h_data.extent(1));
    Kokkos::deep_copy(atm_data.data, h_data);
  }

} // anonymous namespace


ELM::AtmForcObjects::AtmForcObjects(const std::string& filename,
                                    const ELM::Utils::Date &file_start_time,
                                    const int atm_nsteps, const int ncells)
    :
      forc_TBOT(create_forc_util<AtmForcType::TBOT>(filename, file_start_time, atm_nsteps, ncells)),
      forc_PBOT(create_forc_util<AtmForcType::PBOT>(filename, file_start_time, atm_nsteps, ncells)),
      forc_QBOT(create_forc_util<AtmForcType::QBOT>(filename, file_start_time, atm_nsteps, ncells)),
      forc_FLDS(create_forc_util<AtmForcType::FLDS>(filename, file_start_time, atm_nsteps, ncells)),
      forc_FSDS(create_forc_util<AtmForcType::FSDS>(filename, file_start_time, atm_nsteps, ncells)),
      forc_PREC(create_forc_util<AtmForcType::PREC>(filename, file_start_time, atm_nsteps, ncells)),
      forc_WIND(create_forc_util<AtmForcType::WIND>(filename, file_start_time, atm_nsteps, ncells)),
      forc_ZBOT(create_forc_util<AtmForcType::ZBOT>(filename, file_start_time, atm_nsteps, ncells))
    {}


void ELM::read_forcing(const std::shared_ptr<ELM::AtmForcObjects>& atm_forcing,
                  const Utils::DomainDecomposition<2>& dd,
                  const ELM::Utils::Date& current,
                  const int atm_nsteps)
{
  read_atm_data(atm_forcing->forc_TBOT, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_PBOT, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_QBOT, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_FLDS, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_FSDS, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_PREC, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_WIND, dd, current, atm_nsteps);
  read_atm_data(atm_forcing->forc_ZBOT, dd, current, atm_nsteps);
}


void ELM::get_forcing(const std::shared_ptr<ELM::AtmForcObjects>& atm_forcing,
                      const std::shared_ptr<ELM::ELMState<ViewI1, ViewI2, ViewD1, ViewD2, ViewD3, ViewPSN1>>& S,
                      const double& model_dt, const ELM::Utils::Date& time_plus_half_dt)
{
  atm_forcing->forc_TBOT.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_tbot, S->forc_thbot);
  atm_forcing->forc_PBOT.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_pbot);
  atm_forcing->forc_QBOT.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_tbot, S->forc_pbot, S->forc_qbot, S->forc_rh);
  atm_forcing->forc_FLDS.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_pbot, S->forc_qbot, S->forc_tbot, S->forc_lwrad);
  atm_forcing->forc_FSDS.get_atm_forcing(model_dt, time_plus_half_dt, S->coszen, S->forc_solai, S->forc_solad);
  atm_forcing->forc_PREC.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_tbot, S->forc_rain, S->forc_snow);
  atm_forcing->forc_WIND.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_u, S->forc_v);
  atm_forcing->forc_ZBOT.get_atm_forcing(model_dt, time_plus_half_dt, S->forc_hgt, S->forc_hgt_u, S->forc_hgt_t,  S->forc_hgt_q);
}

