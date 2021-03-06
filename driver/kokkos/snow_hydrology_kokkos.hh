
#pragma once

#include "elm_state.h"
#include "pft_data.h"
#include "compile_options.hh"
#include "data_types.hh"
#include "aerosol_physics.h"
#include "aerosol_data.h"

namespace ELM {

  void kokkos_snow_hydrology(const std::shared_ptr<ELMStateType>& S,
                             const std::shared_ptr<ELM::AerosolMasses<ViewD2>>& aerosol_masses,
                             const std::shared_ptr<ELM::AerosolDataManager<ViewD1>>& aerosol_data,
                             const std::shared_ptr<ELM::SnwRdsTable<ViewD3>>& snw_rds_table,
                             const ELM::Utils::Date& time_plus_half_dt,
                             const double& dtime);

}
