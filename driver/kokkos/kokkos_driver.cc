
#include <iostream>
#include <string>
#include <unordered_map>

// utilities
#include "helper_functions.hh"

// constants
#include "elm_constants.h"

// initialization routines
#include "initialize_elm_kokkos.hh"
#include "init_timestep_kokkos.hh"

#include "conserved_quantity_kokkos.hh"

// conditional compilation options
#include "compile_options.hh"
#include "data_types.hh"

#include "elm_kokkos_interface.hh"

using ELM::Utils::assign;
using AtmForcType = ELM::AtmForcType;

int main(int argc, char **argv) {

{ // enclosing scope

  using namespace ELM::ELMdims;

  Kokkos::initialize(argc, argv);

  { // inner scope

    const int ncols = 1;
    int idx = 0; // hardwire for ncols = 1
    const int ntimes = 100;
    const double dtime = 1800.0;
    const double dtime_d = 1800.0 / 86400.0;
    const auto start = ELM::Utils::Date(2014, 1, 1);

    // instantiate the ELM interface
    ELM::ELMInterface elm_interface(ncols);
    // setup model - need to pass input data file paths as function args
    elm_interface.setup();

    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--*/
    /*             TIME LOOP                      */
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--*/
    ELM::Utils::Date current(start);

    for (int t = 0; t < ntimes; ++t) {

      // advance dtime seconds
      elm_interface.advance(current, dtime);

      std::cout << "Output at step: " << t << "\n";
      
      // get views of primary variables
      auto state_ptr = elm_interface.getPrimaryVars();
      
      // print depth-resolved variables
      for (int i = 0; i < nlevsno() + nlevgrnd(); ++i)
      std::cout << "column vars:  " << i <<
       "  t_soisno:  " << state_ptr->t_soisno(0, i) <<
       "  h2osoi_ice:  " << state_ptr->h2osoi_ice(0, i) <<
       "  h2osoi_liq:  " << state_ptr->h2osoi_liq(0, i) << "\n";

      // print column 
      std::cout << "snl " << state_ptr->snl(0) << "\n";
      std::cout << "snow_depth " << state_ptr->snow_depth(0) << "\n";
      std::cout << "frac_sno " << state_ptr->frac_sno(0) << "\n";
      std::cout << "int_snow " << state_ptr->int_snow(0) << "\n";
      std::cout << "h2ocan " << state_ptr->h2ocan(0) << "\n";
      std::cout << "h2osno " << state_ptr->h2osno(0) << "\n";
      std::cout << "h2osfc " << state_ptr->h2osfc(0) << "\n";
      std::cout << "t_grnd " << state_ptr->t_grnd(0) << "\n";
      std::cout << "t_h2osfc " << state_ptr->t_h2osfc(0) << "\n";
      std::cout << "t_h2osfc_bef " << state_ptr->t_h2osfc_bef(0) << "\n";

      current.increment_seconds(dtime);

    } // time loop
  } // inner scope

  Kokkos::finalize();
} // enclosing scope
return 0;
}
