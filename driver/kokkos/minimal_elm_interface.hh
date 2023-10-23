

#pragma once

namespace ELM {

using Utils::assign;
using ELMStateType = ELMState_ATS<ArrayB1, ArrayI1, ArrayI2, ArrayD1, ArrayD2, ArrayD3, ArrayPSN1>;

class MinimalInterface {

public:

  MinimalInterface(int ncols);
  ~ELMState_ATS() = default;

  void setup(int ncols);
  void initialize();
  bool advance(const Utils::Date& dt_start_date, const double& dt_seconds);
  void copyPrimaryVars(std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> > primary_vars);
  std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> > getPrimaryVars();

protected:
  std::shared_ptr<ELMStateType> S_{nullptr};
  std::shared_ptr<ELM::SnicarData<ViewD1, ViewD2, ViewD3>> snicar_data{nullptr};
  std::shared_ptr<ELM::SnwRdsTable<ViewD3>> snw_rds_table{nullptr};
  std::shared_ptr<ELM::PFTData<ViewD1>> pft_data{nullptr};
  std::shared_ptr<ELM::aero_data::AerosolFileInput<ArrayD1>> aero_input;
  std::shared_ptr<ELM::aero_data::AerosolMasses<ArrayD2>> aero_mass;
  std::shared_ptr<ELM::aero_data::AerosolConcentrations<ArrayD2>> aero_concen;
  std::shared_ptr<ELM::phen_data::PhenologyFileInput<ViewD1>> phen_data{nullptr};
  std::shared_ptr<ELM::atm_data::AtmosphereFileInput<ViewD1>> atm_forcing{nullptr};
  int ncols_;
  int cell_per_col_{ELM::ELMdims::nlevgrnd()};
  
  
  


};
