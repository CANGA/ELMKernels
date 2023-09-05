
#pragma once

#include <string>
#include <memory>

#include "data_types.hh"
#include "date_time.hh"


namespace ELM {
class ELMInterface {

public:

  ELMInterface(size_t ncols);

  void setup();
  bool advance(const ELM::Utils::Date& dt_start_date, const double& dt_seconds);
  void copyPrimaryVars(PrimaryVars<ViewI1, ViewD1, ViewD2>& primary_vars);
  std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> > getPrimaryVars();

protected:
  std::shared_ptr<ELMStateType> S_{nullptr};
  size_t ncols_ = 0, cell_per_col_ = 0;
  std::string fname_surfdata_ = "get rid of me";
  std::string fname_forc_ = "get rid of me, too";
};
} // namespace ELM
