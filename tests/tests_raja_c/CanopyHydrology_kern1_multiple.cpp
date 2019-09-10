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
#include <fstream>
#include <algorithm>
#include "memoryManager.hpp"
#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"

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


int main(int RAJA_UNUSED_ARG(argc), char **RAJA_UNUSED_ARG(argv[]))
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;
  
  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  double dtime = 1800.0;
  RAJA::RangeSegment row_range(0, n_grid_cells);
  RAJA::RangeSegment col_range(0, n_pfts);
 
  // phenology state
 
  double* elai = memoryManager::allocate<double>(n_grid_cells * n_pfts);
  double* esai = memoryManager::allocate<double>(n_grid_cells * n_pfts);

  const int DIM = 2;
  RAJA::View<double, RAJA::Layout<DIM> > h_elai( elai, n_months, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_esai( esai, n_months, n_pfts );
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, h_elai, h_esai);

  // forcing state
 
  double* forc_rain = memoryManager::allocate<double>( n_max_times * n_grid_cells );
  double* forc_snow = memoryManager::allocate<double>( n_max_times * n_grid_cells );
  double* forc_air_temp = memoryManager::allocate<double>( n_max_times * n_grid_cells );

  RAJA::View<double, RAJA::Layout<DIM> > h_forc_rain( forc_rain,n_max_times,n_grid_cells );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_snow(forc_snow,n_max_times,n_grid_cells);
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_air_temp( forc_air_temp,n_max_times,n_grid_cells );

  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  
  double* forc_irrig = memoryManager::allocate<double>( n_max_times*n_grid_cells );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_irrig(forc_irrig,n_max_times,n_grid_cells );
  //ELM::Utils::MatrixForc forc_irrig; forc_irrig = 0.;
  
  // output state by the grid cell
  
  double* qflx_prec_intr= memoryManager::allocate<double>(n_grid_cells*n_pfts);
  double* qflx_irrig= memoryManager::allocate<double>(n_grid_cells*n_pfts);
  double* qflx_prec_grnd= memoryManager::allocate<double>(n_grid_cells*n_pfts);
  double* qflx_snwcp_liq= memoryManager::allocate<double>(n_grid_cells*n_pfts);
  double* qflx_snwcp_ice = memoryManager::allocate<double>(n_grid_cells*n_pfts);
  double* qflx_snow_grnd_patch= memoryManager::allocate<double>(n_grid_cells*n_pfts);
  double* qflx_rain_grnd= memoryManager::allocate<double>(n_grid_cells*n_pfts  );

  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_prec_intr(qflx_prec_intr, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_irrig(qflx_irrig, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_prec_grnd(qflx_prec_grnd, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snwcp_liq(qflx_snwcp_liq, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snwcp_ice( qflx_snwcp_ice, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snow_grnd_patch( qflx_snow_grnd_patch, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_rain_grnd(  qflx_rain_grnd, n_grid_cells, n_pfts  );

  // output state by the pft
  double* h2o_can= memoryManager::allocate<double>(n_grid_cells*n_pfts);
 
   //auto h2o_can = ELM::Utils::MatrixState(); 
  RAJA::View<double, RAJA::Layout<DIM> > h_h2o_can( h2o_can,n_grid_cells,n_pfts );
  std::memset(h2o_can, 0, n_grid_cells*n_pfts * sizeof(double));
  double* end = NULL;
  end = &h_h2o_can(n_grid_cells-1, n_pfts-1) ;

  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_multiple.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  std::cout << " check " << n_grid_cells*n_pfts << std::endl ;
  std::cout << "Size " << end+1 - &h_h2o_can(0,0) << std::endl ;
  auto min_max = std::minmax_element(&h_h2o_can(0,0), end+1);//h2o_can1.begin(), h2o_can1.end());
  soln_file << std::setprecision(16)
          << 0 << "\t" << std::accumulate(&h_h2o_can(0,0), end+1, 0.) //h2o_can1.begin(), h2o_can1.end(), 0.)
          << "\t" << *min_max.first
          << "\t" << *min_max.second << std::endl;

  std::cout << std::setprecision(16)
          << 0 << "\t" << std::accumulate(&h_h2o_can(0,0), end+1, 0.) //h2o_can1.begin(), h2o_can1.end(), 0.)
          << "\t" << *min_max.first
          << "\t" << *min_max.second << std::endl;



  
  // main loop
  // -- the timestep loop cannot/should not be parallelized

  //
  // Define the kernel execution policy
  //


  //
  // Define the kernel
  //

    for (size_t t = 0; t != n_times; ++t) {


    // using EXEC_POL1 =
    // RAJA::KernelPolicy<
    //   RAJA::statement::For<0, RAJA::omp_parallel_for_exec,  // row
    //     RAJA::statement::For<0, RAJA::loop_exec,            // col
    //       RAJA::statement::Lambda<0>
    //     >
    //   >
    // >;

    //   RAJA::kernel<EXEC_POL1>(RAJA::make_tuple(n_grid_cells,n_pfts),
    //   [=] (const size_t& g, const size_t& p) {
      // for (size_t g = 0; g != n_grid_cells; ++g) {
      //   for (size_t p = 0; p != n_pfts; ++p) {
        RAJA::forall<RAJA::loop_exec>( row_range, [=](int g) {

          RAJA::forall<RAJA::loop_exec>( col_range, [=](int p) {
          ELM::CanopyHydrology_Interception(dtime,
                    h_forc_rain(t,g), h_forc_snow(t,g), h_forc_irrig(t,g),
                    ltype, ctype, urbpoi, do_capsnow,
                    h_elai(g,p), h_esai(g,p), dewmx, frac_veg_nosno,
                    h_h2o_can(g,p), n_irrig_steps_left,
                    h_qflx_prec_intr(g,p), h_qflx_irrig(g,p), h_qflx_prec_grnd(g,p),
                    h_qflx_snwcp_liq(g,p), h_qflx_snwcp_ice(g,p),
                    h_qflx_snow_grnd_patch(g,p), h_qflx_rain_grnd(g,p));
        });
      //}
      });
    auto min_max = std::minmax_element(&h_h2o_can(0,0), end+1);//h2o_can1.begin(), h2o_can1.end());
    std::cout << std::setprecision(16)
              << t+1 << "\t" << std::accumulate(&h_h2o_can(0,0), end+1, 0.)//h2o_can1.begin(), h2o_can1.end(), 0.)
              << "\t" << *min_max.first
              << "\t" << *min_max.second << std::endl;
    soln_file << std::setprecision(16)
              << t+1 << "\t" << std::accumulate(&h_h2o_can(0,0), end+1, 0.)//h2o_can1.begin(), h2o_can1.end(), 0.)
              << "\t" << *min_max.first
              << "\t" << *min_max.second << std::endl;

  } soln_file.close();

  memoryManager::deallocate( elai);
  memoryManager::deallocate( esai);
  memoryManager::deallocate( forc_rain);
  memoryManager::deallocate( forc_snow);
  memoryManager::deallocate( forc_irrig);
  memoryManager::deallocate( forc_air_temp);
  memoryManager::deallocate( qflx_prec_intr);
  memoryManager::deallocate( qflx_irrig);
  memoryManager::deallocate( qflx_prec_grnd);
  memoryManager::deallocate( qflx_snwcp_liq);
  memoryManager::deallocate( qflx_snwcp_ice);
  memoryManager::deallocate( qflx_snow_grnd_patch);
  memoryManager::deallocate( qflx_rain_grnd);
  memoryManager::deallocate( h2o_can);
  
  return 0;
}
