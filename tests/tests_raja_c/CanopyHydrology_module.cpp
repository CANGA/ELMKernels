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
#include "CanopyHydrology_SnowWater_impl.hh"

namespace ELM {
namespace Utils {

static const int n_months = 12;
static const int n_pfts = 17;
static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                            // day * half hour timestep
static const int n_grid_cells = 24;
static const int n_levels_snow = 5;

} // namespace
} // namespace


int main(int RAJA_UNUSED_ARG(argc), char **RAJA_UNUSED_ARG(argv[]))
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;
  using ELM::Utils::n_levels_snow;
  
  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  const double dtime = 1800.0;

  // fixed magic parameters for SnowWater
  const double qflx_snow_melt = 0.;

  // fixed magic parameters for fracH2Osfc  
  const int oldfflag = 0;
  const double micro_sigma = 0.1;
  const double min_h2osfc = 1.0e-8;
  const double n_melt = 0.7;
  
  RAJA::RangeSegment row_range(0, n_grid_cells);
  RAJA::RangeSegment col_range(0, n_pfts);
  RAJA::RangeSegment snow_range(0, n_levels_snow);
                            
  // phenology input
  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  
  double* elai = memoryManager::allocate<double>(n_grid_cells * n_pfts);
  double* esai = memoryManager::allocate<double>(n_grid_cells * n_pfts);

  const int DIM = 2;
  //const int DIM1 = 1;

  RAJA::View<double, RAJA::Layout<DIM> > h_elai( elai, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_esai( esai, n_grid_cells, n_pfts );

  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, h_elai, h_esai);

  // forcing input
  double* forc_rain = memoryManager::allocate<double>( n_max_times * n_grid_cells );
  double* forc_snow = memoryManager::allocate<double>( n_max_times * n_grid_cells );
  double* forc_air_temp = memoryManager::allocate<double>( n_max_times * n_grid_cells );

  RAJA::View<double, RAJA::Layout<DIM> > h_forc_rain( forc_rain, n_max_times,n_grid_cells );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_snow( forc_snow, n_max_times,n_grid_cells );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_air_temp( forc_air_temp, n_max_times,n_grid_cells );

  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  
  double* forc_irrig = memoryManager::allocate<double>( n_max_times*n_grid_cells );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_irrig( forc_irrig, n_max_times,n_grid_cells );
  double qflx_floodg = 0.0;

  
  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // // need them to be for these kernels. --etc
  
  double* z = memoryManager::allocate<double>(n_grid_cells*n_levels_snow );
  double* zi= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );
  double* dz= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );

  RAJA::View<double, RAJA::Layout<DIM> > h_z( z, n_grid_cells, n_levels_snow );
  RAJA::View<double, RAJA::Layout<DIM> > h_zi( zi, n_grid_cells, n_levels_snow );
  RAJA::View<double, RAJA::Layout<DIM> > h_dz( dz, n_grid_cells, n_levels_snow );

  // state variables that require ICs and evolve (in/out)
  
  double* h2ocan= memoryManager::allocate<double>( n_grid_cells*n_pfts );
  std::memset(h2ocan, 0, n_grid_cells*n_pfts * sizeof(double));
  double* swe_old= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );
  double* h2osoi_liq= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );
  double* h2osoi_ice= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );
  double* t_soisno= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );
  double* frac_iceold= memoryManager::allocate<double>(n_grid_cells*n_levels_snow );

  RAJA::View<double, RAJA::Layout<DIM> > h_h2ocan( h2ocan, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_swe_old( swe_old, n_grid_cells, n_levels_snow );
  RAJA::View<double, RAJA::Layout<DIM> > h_h2osoi_liq( h2osoi_liq, n_grid_cells, n_levels_snow );
  RAJA::View<double, RAJA::Layout<DIM> > h_h2osoi_ice( h2osoi_ice, n_grid_cells, n_levels_snow  );
  RAJA::View<double, RAJA::Layout<DIM> > h_t_soisno( t_soisno, n_grid_cells, n_levels_snow );
  RAJA::View<double, RAJA::Layout<DIM> > h_frac_iceold( frac_iceold, n_grid_cells, n_levels_snow );

  
  double* t_grnd= memoryManager::allocate<double>(n_grid_cells );
  double* h2osno= memoryManager::allocate<double>(n_grid_cells );
  std::memset(h2osno, 0, n_grid_cells* sizeof(double));
  double* snow_depth= memoryManager::allocate<double>(n_grid_cells );
  int* snow_level= memoryManager::allocate<int>(n_grid_cells );

  //RAJA::Layout<1> layout = RAJA::make_offset_layout<1>({{0}}, {{n_grid_cells}});
  const RAJA::Layout<1> layout(n_grid_cells);

  RAJA::View<double, RAJA::Layout<1> > h_t_grnd(  t_grnd,layout);
  RAJA::View<double, RAJA::Layout<1> > h_h2osno( h2osno,layout);
  RAJA::View<double, RAJA::Layout<1> > h_snow_depth(  snow_depth,layout);
  RAJA::View<int, RAJA::Layout<1> > h_snow_level( snow_level,layout);

  
  double* h2osfc= memoryManager::allocate<double>(n_grid_cells );
  double* frac_h2osfc= memoryManager::allocate<double>(n_grid_cells );
  std::memset(frac_h2osfc, 0, n_grid_cells* sizeof(double));
  RAJA::View<double, RAJA::Layout<1> > h_h2osfc(  h2osfc,layout);
  RAJA::View<double, RAJA::Layout<1> > h_frac_h2osfc( frac_h2osfc,layout);

  
  // output fluxes by pft
  
  double* qflx_prec_intr= memoryManager::allocate<double>( n_grid_cells*n_pfts );
  double* qflx_irrig= memoryManager::allocate<double>( n_grid_cells*n_pfts );
  double* qflx_prec_grnd= memoryManager::allocate<double>( n_grid_cells*n_pfts );;
  double* qflx_snwcp_liq= memoryManager::allocate<double>( n_grid_cells*n_pfts );
  double* qflx_snwcp_ice = memoryManager::allocate<double>( n_grid_cells*n_pfts );
  double* qflx_snow_grnd_patch= memoryManager::allocate<double>( n_grid_cells*n_pfts );
  double* qflx_rain_grnd= memoryManager::allocate<double>( n_grid_cells*n_pfts );

  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_prec_intr( qflx_prec_intr, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_irrig( qflx_irrig, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_prec_grnd( qflx_prec_grnd, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snwcp_liq(  qflx_snwcp_liq, n_grid_cells, n_pfts);
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snwcp_ice( qflx_snwcp_ice, n_grid_cells, n_pfts  );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_snow_grnd_patch( qflx_snow_grnd_patch, n_grid_cells, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_qflx_rain_grnd(  qflx_rain_grnd, n_grid_cells, n_pfts  );

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  
  double* integrated_snow= memoryManager::allocate<double>(n_grid_cells );
  RAJA::View<double, RAJA::Layout<1> > h_integrated_snow(integrated_snow,layout);
  
  // output fluxes, state by the column
  
  double* qflx_snow_grnd_col= memoryManager::allocate<double>(n_grid_cells );
  double* qflx_snow_h2osfc= memoryManager::allocate<double>(n_grid_cells );
  double* qflx_h2osfc2topsoi= memoryManager::allocate<double>(n_grid_cells );
  double* qflx_floodc= memoryManager::allocate<double>(n_grid_cells );

  RAJA::View<double, RAJA::Layout<1> > h_qflx_snow_grnd_col(  qflx_snow_grnd_col,layout);
  RAJA::View<double, RAJA::Layout<1> > h_qflx_snow_h2osfc( qflx_snow_h2osfc,layout);
  RAJA::View<double, RAJA::Layout<1> > h_qflx_h2osfc2topsoi(  qflx_h2osfc2topsoi,layout);
  RAJA::View<double, RAJA::Layout<1> > h_qflx_floodc( qflx_floodc,layout);

  
  double* frac_sno_eff= memoryManager::allocate<double>(n_grid_cells );
  double* frac_sno= memoryManager::allocate<double>(n_grid_cells );

  RAJA::View<double, RAJA::Layout<1> > h_frac_sno_eff(  frac_sno_eff,layout);
  RAJA::View<double, RAJA::Layout<1> > h_frac_sno( frac_sno,layout);
  


  double* end1 = NULL;
  end1 = &h_h2ocan(n_grid_cells-1, n_pfts-1) ;
  double* end2 = NULL;
  end2 = &h_h2osno(n_grid_cells-1) ;
  double* end3 = NULL;
  end3 = &h_frac_h2osfc(n_grid_cells-1) ;
  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_module.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;
  std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;
  auto min_max_water = std::minmax_element(&h_h2ocan(0,0), end1+1);
  auto sum_water = std::accumulate(&h_h2ocan(0,0), end1+1, 0.);

  auto min_max_snow = std::minmax_element(&h_h2osno(0), end2+1);
  auto sum_snow = std::accumulate(&h_h2osno(0), end2+1, 0.);

  auto min_max_frac_sfc = std::minmax_element(&h_frac_h2osfc(0), end3+1);
  auto avg_frac_sfc = std::accumulate(&h_frac_h2osfc(0), end3+1, 0.) / (end3+1 - &h_frac_h2osfc(0));

  soln_file << std::setprecision(16)
            << 0 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
            << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
            << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;

  std::cout << std::setprecision(16)
            << 0 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
            << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
            << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;

  // using POL = RAJA::KernelPolicy<
  //           RAJA::statement::For<1, RAJA::loop_exec,
  //             RAJA::statement::InitLocalMem<RAJA::cpu_tile_mem, RAJA::ParamList<0, 1>,
  //               RAJA::statement::For<0, RAJA::loop_exec,
  //                 RAJA::statement::Lambda<0>
  //               >,
  //               RAJA::statement::For<0, RAJA::loop_exec,
  //                 RAJA::statement::Lambda<1>
  //               >
  //             >
  //           >
  //         >;
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    // grid cell and/or pft loop can be parallelized
    // for (size_t g = 0; g != n_grid_cells; ++g) {

    //   // PFT level operations
    using REDUCE_POL1 = RAJA::seq_reduce;
    RAJA::ReduceSum<REDUCE_POL1, int> seq_sum(0);

    RAJA::forall<RAJA::loop_exec>( row_range, [=](int g) {
    RAJA::forall<RAJA::loop_exec>( col_range, [=](int p) {
    //   for (size_t p = 0; p != n_pfts; ++p) {
      // RAJA::kernel_param<POL> (RAJA::RangeSegment(0,n_pfts),RAJA::make_tuple(RAJA::RangeSegment(0,n_grid_cells)),

      // [=] (size_t g, size_t p) {
        //
        // Calculate interception
        //
        // NOTE: this currently punts on what to do with the qflx variables!
        // Surely they should be either accumulated or stored on PFTs as well.
        // --etc
        ELM::CanopyHydrology_Interception(dtime,
                    h_forc_rain(t,g), h_forc_snow(t,g), h_forc_irrig(t,g),
                    ltype, ctype, urbpoi, do_capsnow,
                    h_elai(g,p), h_esai(g,p), dewmx, frac_veg_nosno,
                    h_h2ocan(g,p), n_irrig_steps_left,
                    h_qflx_prec_intr(g,p), h_qflx_irrig(g,p), h_qflx_prec_grnd(g,p),
                    h_qflx_snwcp_liq(g,p), h_qflx_snwcp_ice(g,p),
                    h_qflx_snow_grnd_patch(g,p), h_qflx_rain_grnd(g,p));
        //printf("%i %i %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n"(g, p, forc_rain(t,g), forc_snow(t,g), elai(g,p), esai(g,p), h2ocan(g,p), qflx_prec_intr(g));

        //
        // Calculate fraction of LAI that is wet vs dry.
        //
        // FIXME: this currently punts on what to do with the fwet/fdry variables.
        // Surely they should be something, as such this is dead code.
        // By the PFT?
        // --etc
        double fwet = 0., fdry = 0.;
        ELM::CanopyHydrology_FracWet(frac_veg_nosno, h_h2ocan(g,p), h_elai(g,p), h_esai(g,p), dewmx, fwet, fdry);

        seq_sum += h_qflx_snow_grnd_patch(g,p); 
      }); // end PFT loop
      // Column level operations

      // NOTE: this is effectively an accumulation kernel/task! --etc
      h_qflx_snow_grnd_col(g) = seq_sum ; 

      // Calculate ?water balance? on the snow column, adding throughfall,
      // removing melt, etc.
      //
      // local outputs
      int newnode;
      ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
              ltype, ctype, urbpoi, do_capsnow, oldfflag,
              h_forc_air_temp(t,g), h_t_grnd(g),
              h_qflx_snow_grnd_col(g), qflx_snow_melt, n_melt, h_frac_h2osfc(g),
              h_snow_depth(g), h_h2osno(g), h_integrated_snow(g), &h_swe_old(g*n_levels_snow,0),
              &h_h2osoi_liq(g*n_levels_snow,0), &h_h2osoi_ice(g*n_levels_snow,0), &h_t_soisno(g*n_levels_snow,0), &h_frac_iceold(g*n_levels_snow,0),
              h_snow_level(g), &h_dz(g*n_levels_snow,0), &h_z(g*n_levels_snow,0), &h_zi(g*n_levels_snow,0), newnode,
              h_qflx_floodc(g), h_qflx_snow_h2osfc(g), h_frac_sno_eff(g), h_frac_sno(g));

      // Calculate Fraction of Water to the Surface?
      //
      // FIXME: Fortran black magic... h2osoi_liq is a vector, but the
      // interface specifies a single double.  For now passing the 0th
      // entry. --etc
      ELM::CanopyHydrology_FracH2OSfc(dtime, min_h2osfc, ltype, micro_sigma,
              h_h2osno(g), h_h2osfc(g), h_h2osoi_liq(g,0), h_frac_sno(g), h_frac_sno_eff(g),
              h_qflx_h2osfc2topsoi(g), h_frac_h2osfc(g));
      
    }); // end grid cell loop

    
    auto min_max_water = std::minmax_element(&h_h2ocan(0,0), end1+1);
    auto sum_water = std::accumulate(&h_h2ocan(0,0), end1+1, 0.);

    auto min_max_snow = std::minmax_element(&h_h2osno(0), end2+1);
    auto sum_snow = std::accumulate(&h_h2osno(0), end2+1, 0.);

    auto min_max_frac_sfc = std::minmax_element(&h_frac_h2osfc(0), end3+1);
    auto avg_frac_sfc = std::accumulate(&h_frac_h2osfc(0), end3+1, 0.) / (end3+1 - &h_frac_h2osfc(0));
     
    std::cout << std::setprecision(16)
              << t+1 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
              << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
              << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;
                          
    soln_file << std::setprecision(16)
              << t+1 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
              << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
              << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;

  } // end timestep loop
  soln_file.close();
  
  memoryManager::deallocate( elai);
  memoryManager::deallocate( esai);
  memoryManager::deallocate( forc_rain);
  memoryManager::deallocate( forc_snow);
  memoryManager::deallocate( forc_air_temp);
  memoryManager::deallocate( forc_irrig);
  memoryManager::deallocate( z);
  memoryManager::deallocate( zi);
  memoryManager::deallocate( dz);
  memoryManager::deallocate( h2ocan);
  memoryManager::deallocate( swe_old);
  memoryManager::deallocate( h2osoi_liq);
  memoryManager::deallocate( h2osoi_ice);
  memoryManager::deallocate( t_soisno);
  memoryManager::deallocate( frac_iceold);
  memoryManager::deallocate( t_grnd);
  memoryManager::deallocate( h2osno);
  memoryManager::deallocate( snow_depth);
  memoryManager::deallocate( snow_level);
  memoryManager::deallocate( h2osfc);
  memoryManager::deallocate( frac_h2osfc);
  memoryManager::deallocate( qflx_prec_intr);
  memoryManager::deallocate( qflx_irrig);
  memoryManager::deallocate( qflx_prec_grnd);
  memoryManager::deallocate( qflx_snwcp_liq);
  memoryManager::deallocate( qflx_snwcp_ice);
  memoryManager::deallocate( qflx_snow_grnd_patch);
  memoryManager::deallocate( qflx_rain_grnd);
  memoryManager::deallocate( integrated_snow);
  memoryManager::deallocate( qflx_snow_grnd_col);
  memoryManager::deallocate( qflx_snow_h2osfc);
  memoryManager::deallocate( qflx_h2osfc2topsoi);
  memoryManager::deallocate( qflx_floodc);
  memoryManager::deallocate( frac_sno_eff);
  memoryManager::deallocate( frac_sno);

  return 0;
}
