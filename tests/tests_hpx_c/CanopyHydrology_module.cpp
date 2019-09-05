#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/parallel_for_loop.hpp>
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
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>
//#include <mpi.h>
#include <chrono>
#include "utils.hh"
#include "readers.hh"


#include "CanopyHydrology.hh"
#include "CanopyHydrology_SnowWater_impl.hh"
using namespace std::chrono; 


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


int main(int argc, char ** argv)
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;
  using ELM::Utils::n_levels_snow;
  // int myrank, numprocs;
  // double mytime, maxtime, mintime, avgtime;

  // MPI_Init(&argc,&argv);
  // MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  // MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  // MPI_Barrier(MPI_COMM_WORLD);
  
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
                               
  // phenology input
  double* elai = new double[n_grid_cells * n_pfts];
  double* esai = new double[n_grid_cells * n_pfts];
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, elai, esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, elai, esai);

  // forcing input
  double* forc_rain = new double[ n_max_times * n_grid_cells ];
  double* forc_snow = new double[ n_max_times * n_grid_cells ];
  double* forc_air_temp = new double[ n_max_times * n_grid_cells ];
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, forc_rain, forc_snow, forc_air_temp);
  double* forc_irrig = new double[ n_max_times * n_grid_cells ];
  double qflx_floodg = 0.0;

  
  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // need them to be for these kernels. --etc
  auto z = new double[n_grid_cells*n_levels_snow];
  auto zi = new double[n_grid_cells*n_levels_snow];
  auto dz = new double[n_grid_cells*n_levels_snow];

  // state variables that require ICs and evolve (in/out)
  double* h2ocan = new double[n_grid_cells * n_pfts]; 
  double* swe_old = new double[n_grid_cells*n_levels_snow];
  double* h2osoi_liq = new double[n_grid_cells*n_levels_snow];
  double* h2osoi_ice = new double[n_grid_cells*n_levels_snow];
  double* t_soisno = new double[n_grid_cells*n_levels_snow];
  double* frac_iceold = new double[n_grid_cells*n_levels_snow];
  double* t_grnd = new double[n_grid_cells];
  double* h2osno = new double[n_grid_cells]; 
  double* snow_depth = new double[n_grid_cells];
  int* snow_level = new int[n_grid_cells]; // note this tracks the snow_depth

  double* h2osfc = new double[n_grid_cells];
  double* frac_h2osfc = new double[n_grid_cells]; 

  
  // output fluxes by pft
  double* qflx_prec_intr = new double[n_grid_cells * n_pfts];
  double* qflx_irrig = new double[n_grid_cells * n_pfts];
  double* qflx_prec_grnd = new double[n_grid_cells * n_pfts];
  double* qflx_snwcp_liq = new double[n_grid_cells * n_pfts];
  double* qflx_snwcp_ice = new double[n_grid_cells * n_pfts];
  double* qflx_snow_grnd_patch = new double[n_grid_cells * n_pfts];
  double* qflx_rain_grnd = new double[n_grid_cells * n_pfts];

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  double* integrated_snow = new double[n_grid_cells];
  
  // output fluxes, state by the column
  double* qflx_snow_grnd_col = new double[n_grid_cells];
  double* qflx_snow_h2osfc = new double[n_grid_cells];
  double* qflx_h2osfc2topsoi = new double[n_grid_cells];
  double* qflx_floodc = new double[n_grid_cells];

  double* frac_sno_eff = new double[n_grid_cells];
  double* frac_sno = new double[n_grid_cells];

  double* end = &h2ocan[n_grid_cells, n_pfts] ;
  double* end2 = &h2osno[n_grid_cells-1] ;
  double* end3 = &frac_h2osfc[n_grid_cells-1] ;

  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_module.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;
  auto min_max_water = std::minmax_element(&h2ocan[0,0], end+1);
  auto sum_water = std::accumulate(&h2ocan[0,0], end+1, 0.);

  auto min_max_snow = std::minmax_element(&h2osno[0], end2+1);
  auto sum_snow = std::accumulate(&h2osno[0], end2+1, 0.);

  auto min_max_frac_sfc = std::minmax_element(&frac_h2osfc[0], end3+1);
  auto avg_frac_sfc = std::accumulate(&frac_h2osfc[0], end3+1, 0.) / (end3+1 - &frac_h2osfc[0]);
  
  soln_file << std::setprecision(16)
            << 0 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
            << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
            << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;
  auto start = high_resolution_clock::now();
  // mytime = MPI_Wtime();
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    // grid cell and/or pft loop can be parallelized
    for (size_t g = 0; g != n_grid_cells; ++g) {
      double sum;

      // PFT level operations
      for (size_t p = 0; p != n_pfts; ++p) {
        //
        // Calculate interception
        //
        // NOTE: this currently punts on what to do with the qflx variables!
        // Surely they should be either accumulated or stored on PFTs as well.
        // --etc
        ELM::CanopyHydrology_Interception(dtime,
                forc_rain[t,g], forc_snow[t,g], forc_irrig[t,g],
                ltype, ctype, urbpoi, do_capsnow,
                elai[g,p], esai[g,p], dewmx, frac_veg_nosno,
                h2ocan[g,p], n_irrig_steps_left,
                qflx_prec_intr[g,p], qflx_irrig[g,p], qflx_prec_grnd[g,p],
                qflx_snwcp_liq[g,p], qflx_snwcp_ice[g,p],
                qflx_snow_grnd_patch[g,p], qflx_rain_grnd[g,p]);
        //printf("%i %i %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", g, p, forc_rain[t,g], forc_snow[t,g], elai[g,p], esai[g,p], h2ocan[g,p], qflx_prec_intr[g]);

        //
        // Calculate fraction of LAI that is wet vs dry.
        //
        // FIXME: this currently punts on what to do with the fwet/fdry variables.
        // Surely they should be something, as such this is dead code.
        // By the PFT?
        // --etc
        double fwet = 0., fdry = 0.;
        ELM::CanopyHydrology_FracWet(frac_veg_nosno, h2ocan[g,p], elai[g,p], esai[g,p], dewmx, fwet, fdry);

        sum += qflx_snow_grnd_patch[g,p];
      } // end PFT loop

      qflx_snow_grnd_col[g] = sum ;

      // Column level operations

      // NOTE: this is effectively an accumulation kernel/task! --etc
      

      // Calculate ?water balance? on the snow column, adding throughfall,
      // removing melt, etc.
      //
      // local outputs
      int newnode;
      ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
              ltype, ctype, urbpoi, do_capsnow, oldfflag,
              forc_air_temp[t,g], t_grnd[g],
              qflx_snow_grnd_col[g], qflx_snow_melt, n_melt, frac_h2osfc[g],
              snow_depth[g], h2osno[g], integrated_snow[g], swe_old[g],
              h2osoi_liq[g], h2osoi_ice[g], t_soisno[g], frac_iceold[g],
              snow_level[g], dz[g], z[g], zi[g], newnode,
              qflx_floodc[g], qflx_snow_h2osfc[g], frac_sno_eff[g], frac_sno[g]);

      // Calculate Fraction of Water to the Surface?
      //
      // FIXME: Fortran black magic... h2osoi_liq is a vector, but the
      // interface specifies a single double.  For now passing the 0th
      // entry. --etc
      ELM::CanopyHydrology_FracH2OSfc(dtime, min_h2osfc, ltype, micro_sigma,
              h2osno[g], h2osfc[g], h2osoi_liq[g,0], frac_sno[g], frac_sno_eff[g],
              qflx_h2osfc2topsoi[g], frac_h2osfc[g]);
      
    } // end grid cell loop

    auto min_max_water = std::minmax_element(&h2ocan[0,0], end+1);
    auto sum_water = std::accumulate(&h2ocan[0,0], end+1, 0.);

    auto min_max_snow = std::minmax_element(&h2osno[0], end2+1);
    auto sum_snow = std::accumulate(&h2osno[0], end2+1, 0.);

    auto min_max_frac_sfc = std::minmax_element(&frac_h2osfc[0], end3+1);
    auto avg_frac_sfc = std::accumulate(&frac_h2osfc[0], end3+1, 0.) / (end3+1 - &frac_h2osfc[0]);
                  
    soln_file << std::setprecision(16)
              << t+1 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
              << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
              << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;
  } // end timestep loop
  // mytime = MPI_Wtime() - mytime;
  auto stop = high_resolution_clock::now();
  // std::cout <<"Timing from node "<< myrank  << " is "<< mytime << "seconds." << std::endl;

  // MPI_Reduce(&mytime, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
  // MPI_Reduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
  // MPI_Reduce(&mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
  // if (myrank == 0) {
  // avgtime /= numprocs;
  // std::cout << "Min: "<< mintime <<  ", Max: " << maxtime << ", Avg: " <<avgtime << std::endl;
  // }

  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl;
  // MPI_Finalize();
  return 0;
}
