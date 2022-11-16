
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

  void kokkos_snow_hydrology(ELMStateType& S,
                             ELM::AerosolMasses<ViewD2>& aerosol_masses,
                             ELM::AerosolDataManager<ViewD1>& aerosol_data,
                             const ELM::SnwRdsTable<ViewD3>& snw_rds_table,
                             const ELM::Utils::Date& time_plus_half_dt,
                             const double& dtime);

}
