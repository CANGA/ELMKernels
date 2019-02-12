#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>

#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"


int main(int argc, char ** argv)
{
  // dimensions
  const int n_months = 12;
  const int n_pfts = 17;
  const int n_grid_cells = 24;
  const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                       // day * half hour timestep

  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  const double dtime = 1800.0;

  // phenology state
  ELM::Utils::Matrix<> elai(n_grid_cells, n_pfts);
  ELM::Utils::Matrix<> esai(n_grid_cells, n_pfts);
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, elai, esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, elai, esai);

  // forcing state
  ELM::Utils::Matrix<> forc_rain(n_max_times, n_grid_cells);
  ELM::Utils::Matrix<> forc_snow(n_max_times, n_grid_cells);
  ELM::Utils::Matrix<> forc_air_temp(n_max_times, n_grid_cells);
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, forc_rain, forc_snow, forc_air_temp);

  ELM::Utils::Matrix<> forc_irrig(n_max_times, n_grid_cells); forc_irrig = 0.;
  
  // output state by the grid cell
  auto qflx_prec_intr = std::vector<double>(n_grid_cells, 0.);
  auto qflx_irrig = std::vector<double>(n_grid_cells, 0.);
  auto qflx_prec_grnd = std::vector<double>(n_grid_cells, 0.);
  auto qflx_snwcp_liq = std::vector<double>(n_grid_cells, 0.);
  auto qflx_snwcp_ice = std::vector<double>(n_grid_cells, 0.);
  auto qflx_snow_grnd_patch = std::vector<double>(n_grid_cells, 0.);
  auto qflx_rain_grnd = std::vector<double>(n_grid_cells, 0.);

  // output state by the pft
  auto h2o_can = ELM::Utils::Matrix<>(n_grid_cells, n_pfts); h2o_can = 0.;

  std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  auto min_max = std::minmax_element(h2o_can.begin(), h2o_can.end());
  std::cout << std::setprecision(16)
            << 0 << "\t" << std::accumulate(h2o_can.begin(), h2o_can.end(), 0.)
            << "\t" << *min_max.first
            << "\t" << *min_max.second << std::endl;
  
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    // grid cell and/or pft loop can be parallelized
    for (size_t g = 0; g != n_grid_cells; ++g) {
      for (size_t p = 0; p != n_pfts; ++p) {
        // NOTE: this currently punts on what to do with the qflx variables!
        // Surely they should be either accumulated or stored on PFTs as well.
        // --etc
        ELM::CanopyHydrologyKern1(dtime,
                forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                ltype, ctype, urbpoi, do_capsnow,
                elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                h2o_can(g,p), n_irrig_steps_left,
                qflx_prec_intr[g], qflx_irrig[g], qflx_prec_grnd[g],
                qflx_snwcp_liq[g], qflx_snwcp_ice[g],
                qflx_snow_grnd_patch[g], qflx_rain_grnd[g]);
        //printf("%i %i %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", g, p, forc_rain(t,g), forc_snow(t,g), elai(g,p), esai(g,p), h2o_can(g,p), qflx_prec_intr[g]);
      }
    }

    auto min_max = std::minmax_element(h2o_can.begin(), h2o_can.end());
    std::cout << std::setprecision(16)
              << t+1 << "\t" << std::accumulate(h2o_can.begin(), h2o_can.end(), 0.)
              << "\t" << *min_max.first
              << "\t" << *min_max.second << std::endl;

  }
  return 0;
}
