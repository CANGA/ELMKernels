
#pragma once

#include "atm_data.h"
#include "elm_state.h"

namespace ELM {

  // need to put a check in to automatically choose 
  // RH or QBOT
  // will currently error if a file with QBOT instead of RH is used 
  struct AtmForcObjects {
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::TBOT> forc_TBOT;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::PBOT> forc_PBOT;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::RH> forc_QBOT;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::FLDS> forc_FLDS;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::FSDS> forc_FSDS;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::PREC> forc_PREC;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::WIND> forc_WIND;
    ELM::AtmDataManager<ViewD1, ViewD2, AtmForcType::ZBOT> forc_ZBOT;

    AtmForcObjects(const std::string& filename,
                   const ELM::Utils::Date &file_start_time,
                   const int ntimes, const int ncells);
    ~AtmForcObjects() = default;
  };

  void read_forcing(const std::shared_ptr<ELM::AtmForcObjects>& atm_forcing,
                    const Utils::DomainDecomposition<2>& dd,
                    const ELM::Utils::Date& current,
                    const int atm_nsteps);

  void get_forcing(const std::shared_ptr<ELM::AtmForcObjects>& atm_forcing,
                   const std::shared_ptr<ELM::ELMState<ViewI1, ViewI2, ViewD1, ViewD2, ViewD3, ViewPSN1>>& S,
                   const double& model_dt, const ELM::Utils::Date& time_plus_half_dt);

} // namespace ELM
