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
#include "landunit_varcon.h"
#include "column_varcon.h" 
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

 
  // phenology state
  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  double* elai = new double [n_months * n_pfts];
  double* esai = new double [n_months * n_pfts];

  const int DIM = 2;
  RAJA::View<double, RAJA::Layout<DIM> > h_elai( elai, n_months, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_esai( esai, n_months, n_pfts );
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, h_elai, h_esai);

  // forcing state
  // ELM::Utils::MatrixForc forc_rain;
  // ELM::Utils::MatrixForc forc_snow;
  // ELM::Utils::MatrixForc forc_air_temp;

  double* forc_rain = new double [ n_max_times*n_grid_cells ];
  double* forc_snow = new double [ n_max_times*n_grid_cells ];
  double* forc_air_temp = new double [ n_max_times*n_grid_cells ];
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_rain( forc_rain,n_max_times,n_grid_cells );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_snow(forc_snow,n_max_times,n_grid_cells);
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_air_temp( forc_air_temp,n_max_times,n_grid_cells );
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  double* forc_irrig = new double [ n_max_times*n_grid_cells ];
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_irrig(forc_irrig,n_max_times,n_grid_cells );
  //ELM::Utils::MatrixForc forc_irrig; forc_irrig = 0.;
  
  // output state by the grid cell
  // auto qflx_prec_intr = std::array<double,n_grid_cells>();
  // auto qflx_irrig = std::array<double,n_grid_cells>();
  // auto qflx_prec_grnd = std::array<double,n_grid_cells>();
  // auto qflx_snwcp_liq = std::array<double,n_grid_cells>();
  // auto qflx_snwcp_ice = std::array<double,n_grid_cells>();
  // auto qflx_snow_grnd_patch = std::array<double,n_grid_cells>();
  // auto qflx_rain_grnd = std::array<double,n_grid_cells>();
  // auto qflx_prec_intr = ELM::Utils::MatrixState();
  // auto qflx_irrig = ELM::Utils::MatrixState();
  // auto qflx_prec_grnd = ELM::Utils::MatrixState();
  // auto qflx_snwcp_liq = ELM::Utils::MatrixState();
  // auto qflx_snwcp_ice = ELM::Utils::MatrixState();
  // auto qflx_snow_grnd_patch = ELM::Utils::MatrixState();
  // auto qflx_rain_grnd = ELM::Utils::MatrixState();
  double* qflx_prec_intr= new double [n_grid_cells*n_pfts];
  double* qflx_irrig= new double [n_grid_cells*n_pfts];
  double* qflx_prec_grnd= new double [n_grid_cells*n_pfts];
  double* qflx_snwcp_liq= new double [n_grid_cells*n_pfts];
  double* qflx_snwcp_ice = new double [n_grid_cells*n_pfts];
  double* qflx_snow_grnd_patch= new double [n_grid_cells*n_pfts];
  double* qflx_rain_grnd= new double [qflx_rain_grnd, n_grid_cells, n_pfts  ];
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_prec_intr(qflx_prec_intr, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_irrig(qflx_irrig, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_prec_grnd(qflx_prec_grnd, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snwcp_liq(qflx_snwcp_liq, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snwcp_ice( qflx_snwcp_ice, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snow_grnd_patch( qflx_snow_grnd_patch, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_rain_grnd(  qflx_rain_grnd, n_grid_cells, n_pfts  );

  // output state by the pft
  double* h2o_can= new double [n_grid_cells*n_pfts];
  //auto h2o_can = ELM::Utils::MatrixState(); 
  RAJA::View<double, RAJA::Layout<DIM> > h_h2o_can( h2o_can,n_grid_cells,n_pfts );
  double* end = &h_h2o_can(n_grid_cells-1, n_pfts-1) ;

  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_multiple.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
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
  for (size_t t = 0; t != n_times; ++t) {
    Kokkos::parallel_for("n_grid_cells", n_grid_cells, KOKKOS_LAMBDA (const size_t& g) {
      for (size_t p = 0; p != n_pfts; ++p) {
        ELM::CanopyHydrology_Interception(dtime,
                forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                ltype, ctype, urbpoi, do_capsnow,
                elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                h2o_can(g,p), n_irrig_steps_left,
                qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p));

                
      }
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
