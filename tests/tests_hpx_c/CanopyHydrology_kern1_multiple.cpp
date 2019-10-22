#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>
//#include <mpi.h>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/parallel_for_loop.hpp>
#include <hpx/components/containers/partitioned_vector/partitioned_vector_view.hpp>
#include <hpx/include/partitioned_vector.hpp>
#include <chrono>
#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"
using namespace std::chrono; 

namespace ELM {
namespace Utils {

static const int n_months = 12;
static const int n_pfts = 17;
static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
// day * half hour timestep
static const int n_grid_cells = 24;

using MatrixState = MatrixStatic<n_grid_cells, n_pfts>;
using MatrixForc = MatrixStatic<n_max_times,n_grid_cells>;

} // namespace
} // namespace

HPX_REGISTER_PARTITIONED_VECTOR(double);

int main(int argc, char ** argv)
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;

  using Vec = hpx::partitioned_vector<float>;
  using View_2D = hpx::partitioned_vector_view<float,2>;
 
  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  double dtime = 1800.0;

  // phenology state
  ELM::Utils::MatrixState elai;
  ELM::Utils::MatrixState esai;

  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, elai, esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, elai, esai);

  // forcing state
  ELM::Utils::MatrixForc forc_rain;
  ELM::Utils::MatrixForc forc_snow;
  ELM::Utils::MatrixForc forc_air_temp;

  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, forc_rain, forc_snow, forc_air_temp);

 ELM::Utils::MatrixForc forc_irrig; forc_irrig = 0.;
  
  // output state by the grid cell
  
  auto qflx_prec_intr = ELM::Utils::MatrixState();
  auto qflx_irrig = ELM::Utils::MatrixState();
  auto qflx_prec_grnd = ELM::Utils::MatrixState();
  auto qflx_snwcp_liq = ELM::Utils::MatrixState();
  auto qflx_snwcp_ice = ELM::Utils::MatrixState();
  auto qflx_snow_grnd_patch = ELM::Utils::MatrixState();
  auto qflx_rain_grnd = ELM::Utils::MatrixState();


  // output state by the pft
  auto h2o_can = ELM::Utils::MatrixState(); h2o_can = 0.;

  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_multiple.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  auto min_max = std::minmax_element(h2o_can.begin(), h2o_can.end());
  soln_file << std::setprecision(16)
          << 0 << "\t" << std::accumulate(h2o_can.begin(), h2o_can.end(), 0.) 
          << "\t" << *min_max.first
          << "\t" << *min_max.second << std::endl;

  std::cout << std::setprecision(16)
          << 0 << "\t" << std::accumulate(h2o_can.begin(), h2o_can.end(), 0.) 
          << "\t" << *min_max.first
          << "\t" << *min_max.second << std::endl;
  
  auto start = high_resolution_clock::now();
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    // grid cell and/or pft loop can be parallelized
    hpx::parallel::for_loop(hpx::parallel::execution::par.with(
                                hpx::parallel::execution::static_chunk_size()),
                            0, n_grid_cells, [&](const size_t g) {
      hpx::parallel::for_loop(hpx::parallel::execution::par.with(
                                hpx::parallel::execution::static_chunk_size()),
                            0, n_pfts, [&](const size_t p) {
    //for (size_t g = 0; g != n_grid_cells; ++g) {
    //for (size_t p = 0; p != n_pfts; ++p) {
        // NOTE: this currently punts on what to do with the qflx variables!
        // Surely they should be either accumulated or stored on PFTs as well.
        // --etc
        ELM::CanopyHydrology_Interception(dtime,
                forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                ltype, ctype, urbpoi, do_capsnow,
                elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                h2o_can(g,p), n_irrig_steps_left,
                qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p));
          });
    });

    
    auto min_max = std::minmax_element(h2o_can.begin(), h2o_can.end());
    std::cout << std::setprecision(16)
              << t+1 << "\t" << std::accumulate(h2o_can.begin(), h2o_can.end(), 0.)
              << "\t" << *min_max.first
              << "\t" << *min_max.second << std::endl;
    soln_file << std::setprecision(16)
              << t+1 << "\t" << std::accumulate(h2o_can.begin(), h2o_can.end(), 0.)
              << "\t" << *min_max.first
              << "\t" << *min_max.second << std::endl;

  } soln_file.close();

  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl;
  return hpx::finalize();
}
