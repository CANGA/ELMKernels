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
#include <mpi.h>
#include <chrono>

#include "utils.hh"
#include "array.hh"
#include "readers_decl.hh"


#include "CanopyHydrology.hh"
#include "CanopyHydrology_SnowWater_impl.hh"
using namespace std::chrono; 

int main(int argc, char ** argv)
{
  // MPI_Init, etc
  int myrank, n_procs;
  double mytime, maxtime, mintime, avgtime;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Barrier(MPI_COMM_WORLD);

  // get ranks in x, y
  std::size_t nx_procs, ny_procs;
  std::tie(nx_procs, ny_procs) =
      ELM::Utils::get_domain_decomposition(n_procs, argc, argv);

  // NOTE: _global indicates values that are across all ranks.  The absence of
  // global means the variable is spatially local.
  const std::size_t start_year = 2014;
  const std::size_t start_month = 1;
  const std::size_t n_months = 12;
  const std::size_t n_pfts = 17;
  
  const std::string files = "location_of_data";
  
  const auto problem_dims = ELM::IO::get_dimensions(files, start_year, start_month, n_months);
  const std::size_t n_times = std::get<0>(problem_dims);
  const std::size_t nx_global = std::get<1>(problem_dims);
  const std::size_t ny_global = std::get<2>(problem_dims);

  // domain decomposition
  assert(nx_global % nx_procs == 0 && "Currently expect perfectly divisible decomposition.");
  assert(ny_global % ny_procs == 0 && "Currently expect perfectly divisible decomposition.");

  // -- number of local grid cells per process
  const std::size_t nx_local = nx_global / nx_procs;
  const std::size_t ny_local = ny_global / ny_procs;
  const std::size_t n_grid_cells = nx_local * ny_local;

  // -- where am i on the process grid?
  const std::size_t i_proc = myrank % nx_procs;
  const std::size_t j_proc = myrank / nx_procs;

  // -- where do my local unknowns start globally
  const std::size_t i_begin_global = i_proc * nx_local;
  const std::size_t j_begin_global = j_proc * ny_local;

  // allocate storage and initialize phenology input data
  // -- allocate
  ELM::Utils::Array<double,3> elai(n_months, n_grid_cells, n_pfts);
  ELM::Utils::Array<double,3> esai(n_months, n_grid_cells, n_pfts);

  {
    // -- reshape to fit the files, creating a view into elai/esai
    auto elai4D = ELM::Utils::reshape(elai, std::array<std::size_t,4>{n_months, nx_local, ny_local, n_pfts});
    auto esai4D = ELM::Utils::reshape(esai, std::array<std::size_t,4>{n_months, nx_local, ny_local, n_pfts});

    // -- read
    ELM::IO::read_phenology(files, "ELAI", start_year, start_month, i_begin_global, j_begin_global, elai4D);
    ELM::IO::read_phenology(files, "ESAI", start_year, start_month, i_begin_global, j_begin_global, esai4D);
  }
  
  // allocate storage and initialize forcing input data
  // -- allocate
  ELM::Utils::Array<double,2> forc_rain(n_times, n_grid_cells); // NOTE (etc): order uncertain?
  ELM::Utils::Array<double,2> forc_snow(n_times, n_grid_cells); // NOTE (etc): order uncertain?
  ELM::Utils::Array<double,2> forc_air_temp(n_times, n_grid_cells); // NOTE (etc): order uncertain?
  ELM::Utils::Array<double,2> forc_irrig(n_times, n_grid_cells, 0.);
  double qflx_floodg = 0.0;

  {
    // -- reshape to fit the files, creating a view into forcing arrays
    auto forc_rain3D = ELM::Utils::reshape(forc_rain, std::array<std::size_t,4>{n_times, nx_local, ny_local});
    auto forc_snow3D = ELM::Utils::reshape(forc_snow, std::array<std::size_t,4>{n_times, nx_local, ny_local});
    auto forc_air_temp3D = ELM::Utils::reshape(forc_air_temp, std::array<std::size_t,4>{n_times, nx_local, ny_local});

    // -- read
    ELM::IO::read_forcing(files, "RAIN", start_year, start_month, i_begin_global, j_begin_global, forc_rain3D);
    ELM::IO::read_forcing(files, "SNOW", start_year, start_month, i_begin_global, j_begin_global, forc_snow3D);
    ELM::IO::read_forcing(files, "AIR_TEMP", start_year, start_month, i_begin_global, j_begin_global, forc_air_temp3D);
  }    

  
  // fixed magic parameters for now
  const int n_levels_snow = 5;
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
                               
  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // need them to be for these kernels. --etc
  auto z = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto zi = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto dz = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);

  // state variables that require ICs and evolve (in/out)
  auto h2ocan = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto swe_old = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto h2osoi_liq = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto h2osoi_ice = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto t_soisno = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto frac_iceold = ELM::Utils::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto t_grnd = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto h2osno = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto snow_depth = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto snow_level = ELM::Utils::Array<int,1>(n_grid_cells, 0); // note this tracks the snow_depth

  auto h2osfc = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto frac_h2osfc = ELM::Utils::Array<double,1>(n_grid_cells, 0.);

  // output fluxes by pft
  auto qflx_prec_intr = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_irrig = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_prec_grnd = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_snwcp_liq = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_snwcp_ice = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_snow_grnd_patch = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_rain_grnd = ELM::Utils::Array<double,2>(n_grid_cells, n_pfts, 0.);

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  auto integrated_snow = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  
  // output fluxes, state by the column
  auto qflx_snow_grnd_col = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto qflx_snow_h2osfc = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto qflx_h2osfc2topsoi = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto qflx_floodc = ELM::Utils::Array<double,1>(n_grid_cells, 0.);

  auto frac_sno_eff = ELM::Utils::Array<double,1>(n_grid_cells, 0.);
  auto frac_sno = ELM::Utils::Array<double,1>(n_grid_cells, 0.);

#ifdef DEBUG
  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_module.soln");
  {
    soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;
    auto min_max_water = std::minmax_element(h2ocan.begin(), h2ocan.end());
    auto sum_water = std::accumulate(h2ocan.begin(), h2ocan.end(), 0.);

    auto min_max_snow = std::minmax_element(h2osno.begin(), h2osno.end());
    auto sum_snow = std::accumulate(h2osno.begin(), h2osno.end(), 0.);

    auto min_max_frac_sfc = std::minmax_element(frac_h2osfc.begin(), frac_h2osfc.end());
    auto avg_frac_sfc = std::accumulate(frac_h2osfc.begin(), frac_h2osfc.end(), 0.) / (frac_h2osfc.end() - frac_h2osfc.begin());
      
    soln_file << std::setprecision(16)
              << 0 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
              << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
              << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;
  }
#endif

  auto start = high_resolution_clock::now();
  mytime = MPI_Wtime();
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {
    // NOTE (etc): check me... is this correct/reasonable?
    int i_month = (int) std::floor((double) t / (365.0 * 8 / 12));

    // grid cell and/or pft loop can be parallelized
    for (size_t g = 0; g != n_grid_cells; ++g) {

      // PFT level operations
      for (size_t p = 0; p != n_pfts; ++p) {
        //
        // Calculate interception
        //
        // NOTE: this currently punts on what to do with the qflx variables!
        // Surely they should be either accumulated or stored on PFTs as well.
        // --etc
        ELM::CanopyHydrology_Interception(dtime,
                forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                ltype, ctype, urbpoi, do_capsnow,
                elai(i_month,g,p), esai(i_month,g,p), dewmx, frac_veg_nosno,
                h2ocan(g,p), n_irrig_steps_left,
                qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p));

        //
        // Calculate fraction of LAI that is wet vs dry.
        //
        // FIXME: this currently punts on what to do with the fwet/fdry variables.
        // Surely they should be something, as such this is dead code.
        // By the PFT?
        // --etc
        double fwet = 0., fdry = 0.;
        ELM::CanopyHydrology_FracWet(frac_veg_nosno, h2ocan(g,p), elai(i_month,g,p), esai(i_month,g,p), dewmx, fwet, fdry);
      } // end PFT loop

      // Column level operations

      // NOTE: this is effectively an accumulation kernel/task! --etc
      qflx_snow_grnd_col[g] = std::accumulate(qflx_snow_grnd_patch[g].begin(),
              qflx_snow_grnd_patch[g].end(), 0.);

      // Calculate ?water balance? on the snow column, adding throughfall,
      // removing melt, etc.
      //
      // local outputs
      int newnode;
      ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
              ltype, ctype, urbpoi, do_capsnow, oldfflag,
              forc_air_temp(t,g), t_grnd(g),
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
              h2osno[g], h2osfc[g], h2osoi_liq[g][0], frac_sno[g], frac_sno_eff[g],
              qflx_h2osfc2topsoi[g], frac_h2osfc[g]);
      
    } // end grid cell loop
#ifdef DEBUG
    auto min_max_water = std::minmax_element(h2ocan.begin(), h2ocan.end());
    auto sum_water = std::accumulate(h2ocan.begin(), h2ocan.end(), 0.);

    auto min_max_snow = std::minmax_element(h2osno.begin(), h2osno.end());
    auto sum_snow = std::accumulate(h2osno.begin(), h2osno.end(), 0.);

    auto min_max_frac_sfc = std::minmax_element(frac_h2osfc.begin(), frac_h2osfc.end());
    auto avg_frac_sfc = std::accumulate(frac_h2osfc.begin(), frac_h2osfc.end(), 0.) / (frac_h2osfc.end() - frac_h2osfc.begin());
                  
    soln_file << std::setprecision(16)
              << t+1 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
              << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
              << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;
#endif
  } // end timestep loop
  auto stop = high_resolution_clock::now();

  mytime = MPI_Wtime() - mytime;
  std::cout <<"Timing from node "<< myrank  << " is "<< mytime << "seconds." << std::endl;

  MPI_Reduce(&mytime, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
  MPI_Reduce(&mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
  if (myrank == 0) {
    avgtime /= n_procs;
    std::cout << "Min: "<< mintime <<  ", Max: " << maxtime << ", Avg: " <<avgtime << std::endl;
  }
  MPI_Finalize();

  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl;
  return 0;
}
