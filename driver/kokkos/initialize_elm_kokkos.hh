

#pragma once


#include "data_types.hh"

#include "pft_data.h"
#include "atm_data.h"
#include "snicar_data.h"
#include "aerosol_data.h"

#include "elm_state.h"


namespace ELM {

void initialize_kokkos_elm(
  ELMStateType& S,
  SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data,
  SnwRdsTable<ViewD3>& snw_rds_table,
  PFTData<ViewD1>& pft_data,
  AerosolDataManager<ViewD1>& aerosol_data,
  const Utils::DomainDecomposition<2>& dd,
  const std::string& fname_surfdata,
  const std::string& fname_param,
  const std::string& fname_snicar,
  const std::string& fname_snowage,
  const std::string& fname_aerosol);

}
