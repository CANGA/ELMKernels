

#pragma once


#include "kokkos_types.hh"

#include "pft_data.h"
#include "atm_data.h"
#include "snicar_data.h"
#include "aerosol_data.h"
//#include "phenology_data.h"

#include "elm_state.hh"


namespace ELM {

void initialize_with_kokkos (
const std::shared_ptr<ELMState<ViewI1, ViewD1, ViewD2, ViewPSN1>>& S,
const std::shared_ptr<SnicarData<ViewD1, ViewD2, ViewD3>>& snicar_data,
const std::shared_ptr<SnwRdsTable<ViewD3>>& snw_rds_table,
const std::shared_ptr<PFTData<ViewD1, ViewD2>>& pft_data,
const std::shared_ptr<AerosolDataManager<ViewD1>>& aerosol_data,
const Utils::DomainDecomposition<2>& dd,
const std::string& fname_surfdata,
const std::string& fname_param,
const std::string& fname_snicar,
const std::string& fname_snowage,
const std::string& fname_aerosol);

}
