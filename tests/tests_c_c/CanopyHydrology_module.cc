#include <iostream>
#include <fstream>
#include <tuple>
#include <functional>
#include <array>
#include <string>

#include "mpi.h"

#include "../utils/utils.hh"
#include "../utils/array.hh"
#include "../utils/read_input.hh"

#include "CanopyHydrology.hh"
#include "CanopyHydrology_SnowWater_impl.hh"

int main(int argc, char ** argv)
{
  // MPI_Init, etc
  int myrank, n_procs;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    std::cout << "CanopyHydrology: MPI" << std::endl
              << "====================" << std::endl
              << "Problem Setup" << std::endl
              << "--------------------" << std::endl
              << " n_procs = " << n_procs << std::endl;
  }
  auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
  if (myrank == 0) {
    std::cout << " proc_decomp: " << proc_decomp[0] << "," << proc_decomp[1] << std::endl;
  }

  // NOTE: _global indicates values that are across all ranks.  The absence of
  // global means the variable is spatially local.
  const auto start = ELM::Utils::Date(2014, 1, 1);
  const int n_months = 1;
  const int ticks_per_day = 8; // 3-hourly data
  const int write_interval = 8 * 12;
  
  const std::string dir_atm = ATM_DATA_LOCATION;
  const std::string dir_elm = ELM_DATA_LOCATION;

  // dimension: time, lat (ny), lon (nx)
  const std::string basename_forc("Precip3Hrly/clmforc.GSWP3.c2011.0.5x0.5.Prec.");
  const auto forc_dims = 
      ELM::IO::get_forcing_dimensions(MPI_COMM_WORLD, dir_atm, basename_forc, "PRECTmms", start, n_months);

  const std::string basename_phen("surfdata_360x720cru_simyr1850_c180216.nc");
  const auto phen_dims = 
      ELM::IO::get_phenology_dimensions(MPI_COMM_WORLD, dir_elm, basename_phen, "MONTHLY_LAI", start, n_months);

  assert(forc_dims[1] == phen_dims[2]); // n_lat_global
  assert(forc_dims[2] == phen_dims[3]); // n_lon_global
  
  const int n_times = forc_dims[0];
  const int n_pfts = phen_dims[1];

  // Create the domain decomposition.

  auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { forc_dims[1], forc_dims[2] },
          { myrank / proc_decomp[2], myrank % proc_decomp[2] });
  
  
  if (myrank == 0) {
    std::cout << " dimensions:" << std::endl
              << "   n_time = " << n_times << std::endl
              << "   global_n_lat,lon = " << dd.n_global[0] << "," << dd.n_global[1] << std::endl
              << "   nprocs_lat,lon = " << dd.n_procs[0] << "," << dd.n_procs[1] << std::endl
              << "   local_start = " << dd.start[0] << "," << dd.start[1] << std::endl
              << "   local_n_lat,lon = " << dd.n_local[0] << "," << dd.n_local[1] << std::endl;
  }
  auto n_grid_cells = dd.n_local[0] * dd.n_local[1];

  // allocate storage and initialize phenology input data
  // -- allocate
  ELM::Array<double,3> elai(n_months, n_grid_cells, n_pfts);
  ELM::Array<double,3> esai(n_months, n_grid_cells, n_pfts);
  {
    // -- read
    ELM::IO::read_and_reshape_phenology(dir_elm, basename_phen, "MONTHLY_LAI",
            start, n_months, dd, elai);
    if (myrank == 0) {
      std::cout << "File I/O" << std::endl
                << "--------------------" << std::endl
                << "  Phenology LAI read" << std::endl;
    }

    ELM::IO::read_and_reshape_phenology(dir_elm, basename_phen, "MONTHLY_SAI",
            start, n_months, dd, esai);
    if (myrank == 0) {
      std::cout << "  Phenology SAI read" << std::endl;
    }
  }
  
  // allocate storage and initialize forcing input data
  // -- allocate
  ELM::Array<double,2> forc_rain(n_times, n_grid_cells); 
  ELM::Array<double,2> forc_snow(n_times, n_grid_cells); 
  ELM::Array<double,2> forc_air_temp(n_times, n_grid_cells); 
  ELM::Array<double,2> forc_irrig(n_times, n_grid_cells, 0.);
  double qflx_floodg = 0.0;

  {
    // -- reshape to fit the files, creating a view into forcing arrays, and read directly (no copy needed)
    auto forc_rain3D = ELM::reshape(forc_rain, std::array<int,3>{n_times, (int) dd.n_local[0], (int) dd.n_local[1]});
    auto forc_snow3D = ELM::reshape(forc_snow, std::array<int,3>{n_times, (int) dd.n_local[0], (int) dd.n_local[1]});
    auto forc_air_temp3D = ELM::reshape(forc_air_temp, std::array<int,3>{n_times, (int) dd.n_local[0], (int) dd.n_local[1]});

    std::string basename("Precip3Hrly/clmforc.GSWP3.c2011.0.5x0.5.Prec.");
    ELM::IO::read_forcing(dir_atm, basename, "PRECTmms",
                          start, n_months, dd, forc_rain3D);
    if (myrank == 0) std::cout << "  Forcing precip read" << std::endl;

    basename="TPHWL3Hrly/clmforc.GSWP3.c2011.0.5x0.5.TPQWL.";
    ELM::IO::read_forcing(dir_atm, basename, "TBOT",
                          start, n_months, dd, forc_air_temp3D);
    if (myrank == 0) std::cout << "  Forcing air temperature read" << std::endl;
  }    
  ELM::Utils::convert_precip_to_rain_snow(forc_rain, forc_snow, forc_air_temp);
  if (myrank == 0) std::cout << "  Converted precip to rain + snow" << std::endl;
  
  if (myrank == 0) {
    std::cout << "Test Execution" << std::endl
              << "--------------" << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
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
  auto z = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto zi = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto dz = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);

  // state variables that require ICs and evolve (in/out)
  auto h2ocan = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto swe_old = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto h2osoi_liq = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto h2osoi_ice = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto t_soisno = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto frac_iceold = ELM::Array<double,2>(n_grid_cells, n_levels_snow, 0.);
  auto t_grnd = ELM::Array<double,1>(n_grid_cells, 0.);
  auto h2osno = ELM::Array<double,1>(n_grid_cells, 0.);
  auto snow_depth = ELM::Array<double,1>(n_grid_cells, 0.);
  auto snow_level = ELM::Array<int,1>(n_grid_cells, 0); // note this tracks the snow_depth

  auto h2osfc = ELM::Array<double,1>(n_grid_cells, 0.);
  auto frac_h2osfc = ELM::Array<double,1>(n_grid_cells, 0.);

  // output fluxes by pft
  auto qflx_prec_intr = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_irrig = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_prec_grnd = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_snwcp_liq = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_snwcp_ice = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_snow_grnd_patch = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);
  auto qflx_rain_grnd = ELM::Array<double,2>(n_grid_cells, n_pfts, 0.);

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  auto integrated_snow = ELM::Array<double,1>(n_grid_cells, 0.);
  
  // output fluxes, state by the column
  auto qflx_snow_grnd_col = ELM::Array<double,1>(n_grid_cells, 0.);
  auto qflx_snow_h2osfc = ELM::Array<double,1>(n_grid_cells, 0.);
  auto qflx_h2osfc2topsoi = ELM::Array<double,1>(n_grid_cells, 0.);
  auto qflx_floodc = ELM::Array<double,1>(n_grid_cells, 0.);

  auto frac_sno_eff = ELM::Array<double,1>(n_grid_cells, 0.);
  auto frac_sno = ELM::Array<double,1>(n_grid_cells, 0.);

  // for unit testing
  std::ofstream soln_file;
  {
    auto min_max_sum_water = ELM::Utils::min_max_sum(MPI_COMM_WORLD, h2ocan);
    auto min_max_sum_snow = ELM::Utils::min_max_sum(MPI_COMM_WORLD, h2osno);
    auto min_max_sum_surfacewater = ELM::Utils::min_max_sum(MPI_COMM_WORLD, frac_h2osfc);
    if (myrank == 0) std::cout << "  writing ts 0" << std::endl;

    if (myrank == 0) {
      soln_file.open("test_CanopyHydrology_module.soln");
      soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;

      soln_file << std::setprecision(16) << 0
                << "\t" << min_max_sum_water[2] << "\t" << min_max_sum_water[0] << "\t" << min_max_sum_water[1]
                << "\t" << min_max_sum_snow[2] << "\t" << min_max_sum_snow[0] << "\t" << min_max_sum_snow[1]
                << "\t" << min_max_sum_surfacewater[2] << "\t" << min_max_sum_surfacewater[0] << "\t" << min_max_sum_surfacewater[1];
    }
  }

  auto start_wc = ELM::Utils::Clock::time();
  ELM::Utils::Ticker ticker(start, ticks_per_day);

  // main loop
  // -- the timestep loop cannot be parallelized
  for (int t = 0; t != n_times; ++t) {
    int i_month = ELM::Utils::months_since(ticker.now(), ticker.start);

    // grid cell and/or pft loop can be parallelized
    for (int g = 0; g != n_grid_cells; ++g) {

      // PFT level operations
      for (int p = 0; p != n_pfts; ++p) {
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
    if (t % write_interval == 0) {
      auto min_max_sum_water = ELM::Utils::min_max_sum(MPI_COMM_WORLD, h2ocan);
      auto min_max_sum_snow = ELM::Utils::min_max_sum(MPI_COMM_WORLD, h2osno);
      auto min_max_sum_surfacewater = ELM::Utils::min_max_sum(MPI_COMM_WORLD, frac_h2osfc);
      if (myrank == 0) std::cout << "  writing ts " << t << std::endl;

      if (myrank == 0) {
        soln_file << std::setprecision(16)
                  << t << "\t" << min_max_sum_water[2] << "\t" << min_max_sum_water[0] << "\t" << min_max_sum_water[1]
                  << "\t" << min_max_sum_snow[2] << "\t" << min_max_sum_snow[0] << "\t" << min_max_sum_snow[1]
                  << "\t" << min_max_sum_surfacewater[2] << "\t" << min_max_sum_surfacewater[0] << "\t" << min_max_sum_surfacewater[1];
      }
    }
  } // end timestep loop

  auto stop_wc = ELM::Utils::Clock::time();
  auto times = ELM::Utils::Clock::min_max_mean(MPI_COMM_WORLD, stop_wc-start_wc);
  if (myrank == 0) {
    std::cout << "Timing: min: "<< times[0] <<  ", max: " << times[1]
              << ", mean: " << times[2] << std::endl;
  }

  // write final canopy water for kicks
  ELM::IO::reshape_and_write_grid_cell("./final.nc", "snow_depth", dd, snow_depth);
  
  if (myrank == 0) soln_file.close();

  MPI_Finalize();
  return 0;
}


