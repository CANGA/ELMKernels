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
#include <chrono>

#include "Teuchos_Comm.hpp"
#include "Kokkos_Core.hpp"

#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"
#include "CanopyHydrology_SnowWater_impl.hh"

using namespace std::chrono; 

int main(int argc, char ** argv)
{
  // NOTE: _global indicates values that are across all ranks.  The absence of
  // global means the variable is spatially local.
  std::size_t n_times = 365 * 8; // 1 year times 3-hourly timestep
  std::size_t n_months = 12;
  std::size_t nx_global = 360;
  std::size_t ny_global = 180;

  std::size_t nx_ranks = 3;
  std::size_t ny_ranks = 2;

  std::size_t nx = nx_global / nx_ranks; // assumes exact
  std::size_t ny = ny_global / ny_ranks; // assumes exact
  std::size_t n_grid_cells = nx * ny;

  // MPI_Init, etc
  // TODO: convert this to Teuchos::Comm
  int myrank, numprocs;
  double mytime, maxtime, mintime, avgtime;
  MPI_Init(&argc,&argv);
  Kokkos::initialize( argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Barrier(MPI_COMM_WORLD);
  assert(nx_ranks * ny_ranks == n_ranks && "Compile-time sizes set so that code must be run with 6 mpi processes.");
  
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
  using Array1 = Kokkos::View<double*>;
  using Array2 = Kokkos::View<double**>;
  using Array3 = Kokkos::View<double***>;
  using IArray1 = Kokkos::View<int*>;
  using IArray2 = Kokkos::View<int**>;

  // dimensionality... not clear this is right yet --etc
  Array3 elai("elai", n_grid_cells, n_pfts, n_months);
  Array3 esai("esai", n_grid_cells, n_pfts, n_months);
  auto h_elai = Kokkos::create_mirror_view(elai);
  auto h_esai = Kokkos::create_mirror_view(esai);

  // FIX ME (etc)
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  Kokkos::deep_copy(elai, h_elai);
  Kokkos::deep_copy(esai, h_esai);

  // forcing input
  Array2 forc_rain("forc_rain", n_times, n_grid_cells);
  Array2 forc_snow("forc_snow", n_times, n_grid_cells);
  Array2 forc_air_temp("forc_air_temp", n_times, n_grid_cells);
  auto h_forc_rain = Kokkos::create_mirror_view(forc_rain);
  auto h_forc_snow = Kokkos::create_mirror_view(forc_snow);
  auto h_forc_air_temp = Kokkos::create_mirror_view(forc_air_temp);
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  Kokkos::deep_copy(forc_rain, h_forc_rain);
  Kokkos::deep_copy(forc_snow, h_forc_snow);
  Kokkos::deep_copy(forc_air_temp, h_forc_air_temp);
  
  Array2 forc_irrig("forc_irrig", n_times, n_grid_cells);
  double qflx_floodg = 0.0;

  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // // need them to be for these kernels. --etc
  Array2 z("z", n_grid_cells, n_levels_snow);
  Array2 zi("zi", n_grid_cells, n_levels_snow);
  Array2 dz("dz", n_grid_cells, n_levels_snow);

  // state variables that require ICs and evolve (in/out)
  Array2 h2ocan("h2ocan", n_grid_cells, n_pfts);
  Array2 swe_old("swe_old", n_grid_cells, n_levels_snow);
  Array2 h2osoi_liq("h2osoi_liq", n_grid_cells, n_levels_snow);
  Array2 h2osoi_ice("h2osoi_ice", n_grid_cells, n_levels_snow);
  Array2 t_soisno("t_soisno", n_grid_cells, n_levels_snow);
  Array2 frac_iceold("frac_iceold", n_grid_cells, n_levels_snow);
  
  Array1 t_grnd("t_grnd", n_grid_cells);
  Array1 h2osno("h2osno", n_grid_cells);
  Array1 snow_depth("snow_depth", n_grid_cells);
  Array1 snow_level("snow_level", n_grid_cells);

  Array1 h2osfc("h2osfc", n_grid_cells);
  Array1 frac_h2osfc("frac_h2osfc", n_grid_cells);
  
  // output fluxes by pft
  Array2 qflx_prec_intr("qflx_prec_intr", n_grid_cells, n_pfts);
  Array2 qflx_irrig("qflx_irrig", n_grid_cells, n_pfts );
  Array2 qflx_prec_grnd("qflx_prec_grnd", n_grid_cells, n_pfts );
  Array2 qflx_snwcp_liq("qflx_snwcp_liq", n_grid_cells, n_pfts);
  Array2 qflx_snwcp_ice ("qflx_snwcp_ice ", n_grid_cells, n_pfts );
  Array2 qflx_snow_grnd_patch("qflx_snow_grnd_patch", n_grid_cells, n_pfts );
  Array2 qflx_rain_grnd("qflx_rain_grnd", n_grid_cells, n_pfts );

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  Array1 integrated_snow("integrated_snow", n_grid_cells);
  
  // output fluxes, state by the column
  Array1 qflx_snow_grnd_col("qflx_snow_grnd_col", n_grid_cells);
  Array1 qflx_snow_h2osfc("qflx_snow_h2osfc", n_grid_cells);
  Array1 qflx_h2osfc2topsoi("qflx_h2osfc2topsoi", n_grid_cells);
  Array1 qflx_floodc("qflx_floodc", n_grid_cells);
  
  Array1 frac_sno_eff("frac_sno_eff", n_grid_cells);
  Array1 frac_sno("frac_sno", n_grid_cells);

  // for unit testing
  auto min_max_sum_water = min_max_sum(comm, h2ocan);
  auto min_max_sum_snow = min_max_sum(comm, h20sno);
  auto min_max_sum_surfacewater = min_max_sum(comm, frac_h2osfc);

  if (myrank == 0) {
    std::ofstream soln_file;
    soln_file.open("test_CanopyHydrology_module.soln");
    soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;

    soln_file << std::setprecision(16)
              << 0 << "\t" << min_max_sum_water[2] << "\t" << min_max_sum_water[0] << "\t" << min_max_sum_water[1]
              << "\t" << min_max_sum_snow[2] << "\t" << min_max_sum_snow[0] << "\t" << min_max_sum_snow[1]
              << "\t" << min_max_sum_surfacewater[2] << "\t" << min_max_sum_surfacewater[0] << "\t" << min_max_sum_surfacewater[1];
  }

  Kokkos::Timer timer;
  auto start = high_resolution_clock::now();
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    // Column level operations
    // NOTE: this is effectively an accumulation kernel/task! --etc
    typedef Kokkos::TeamPolicy<> team_policy ;
    typedef typename team_policy::member_type team_type ;
    Kokkos::parallel_for(
        "pft_reduction",
        Kokkos::TeamPolicy<> (n_grid_cells, Kokkos::AUTO()),
        KOKKOS_LAMBDA (const team_type& team) {
          auto g = team.league_rank() ;
          double sum = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, n_pfts),
              [=] (const size_t& p, double& lsum) {
                ELM::CanopyHydrology_Interception(dtime,
                        forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                        ltype, ctype, urbpoi, do_capsnow,
                        elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                        h2ocan(g,p), n_irrig_steps_left,
                        qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                        qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                        qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p)); 

                double fwet = 0., fdry = 0.;
                ELM::CanopyHydrology_FracWet(frac_veg_nosno, h2ocan(g,p), elai(g,p), esai(g,p), dewmx, fwet, fdry);

                lsum += qflx_snow_grnd_patch(team.league_rank(),p);
              }, sum);
          qflx_snow_grnd_col(team.league_rank()) = sum ;

          int newnode;
          ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
                  ltype, ctype, urbpoi, do_capsnow, oldfflag,
                  forc_air_temp(t,g), t_grnd(g),
                  qflx_snow_grnd_col(g), qflx_snow_melt, n_melt, frac_h2osfc(g),
                  snow_depth(g), h2osno(g), integrated_snow(g),
                  Kokkos::subview(swe_old, g, Kokkos::ALL),
                  Kokkos::subview(h2osoi_liq, g, Kokkos::ALL),
                  Kokkos::subview(h2osoi_ice, g, Kokkos::ALL),
                  Kokkos::subview(t_soisno, g, Kokkos::ALL),
                  Kokkos::subview(frac_iceold, g, Kokkos::ALL),
                  snow_level(g),
                  Kokkos::subview(dz, g, Kokkos::ALL),
                  Kokkos::subview(z, g, Kokkos::ALL),
                  Kokkos::subview(zi, g, Kokkos::ALL), newnode,
                  qflx_floodc(g), qflx_snow_h2osfc(g), frac_sno_eff(g), frac_sno(g));
          
          // Calculate Fraction of Water to the Surface?
          //
          // FIXME: Fortran black magic... h2osoi_liq is a vector, but the
          // interface specifies a single double.  For now passing the 0th
          // entry. --etc
          ELM::CanopyHydrology_FracH2OSfc(dtime, min_h2osfc, ltype, micro_sigma,
                  h2osno(g), h2osfc(g), h2osoi_liq(g,0), frac_sno(g), frac_sno_eff(g),
                  qflx_h2osfc2topsoi(g), frac_h2osfc(g));
        });
         
    auto min_max_sum_water = min_max_sum(comm, h2ocan);
    auto min_max_sum_snow = min_max_sum(comm, h20sno);
    auto min_max_sum_surfacewater = min_max_sum(comm, frac_h2osfc);

    if (myrank == 0) {
      soln_file << std::setprecision(16)
                << 0 << "\t" << min_max_sum_water[2] << "\t" << min_max_sum_water[0] << "\t" << min_max_sum_water[1]
                << "\t" << min_max_sum_snow[2] << "\t" << min_max_sum_snow[0] << "\t" << min_max_sum_snow[1]
                << "\t" << min_max_sum_surfacewater[2] << "\t" << min_max_sum_surfacewater[0] << "\t" << min_max_sum_surfacewater[1];
    }
  } // end timestep loop
  soln_file.close();

  double time = timer.seconds();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl; 
}

Kokkos::finalize();
MPI_Finalize();
return 0;
}
