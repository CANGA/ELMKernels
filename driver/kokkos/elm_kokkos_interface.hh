
#pragma once

#include <string>
#include <memory>

// conditional compilation options
#include "compile_options.hh"
#include "data_types.hh"

// utilities
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"
#include "helper_functions.hh"

// constants
#include "elm_constants.h"

// initialization routines
#include "initialize_elm_kokkos.hh"
#include "init_timestep_kokkos.hh"

// parallel physics
#include "aerosol_kokkos.hh"
#include "albedo_kokkos.hh"
#include "canopy_hydrology_kokkos.hh"
#include "surface_radiation_kokkos.hh"
#include "canopy_temperature_kokkos.hh"
#include "bareground_fluxes_kokkos.hh"
#include "canopy_fluxes_kokkos.hh"
#include "soil_temperature_kokkos.hh"
#include "snow_hydrology_kokkos.hh"
#include "surface_fluxes_kokkos.hh"
#include "conserved_quantity_kokkos.hh"



namespace ELM {

using Utils::assign;
using AtmForcType = AtmForcType;


class ELMInterface {

public:

  ELMInterface(size_t ncols);

  void setup();
  bool advance(const Utils::Date& dt_start_date, const double& dt_seconds);
  void copyPrimaryVars(std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> > primary_vars);
  std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> > getPrimaryVars();

protected:
  std::shared_ptr<ELMStateType> S_{nullptr};
  size_t ncols_ = 0, cell_per_col_ = 0;
  std::string fname_surfdata_ = "get rid of me";
  std::string fname_forc_ = "get rid of me, too";
  
  
  


};


} // namespace ELM
