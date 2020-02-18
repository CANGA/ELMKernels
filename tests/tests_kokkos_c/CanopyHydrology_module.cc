#include <iostream>
#include <fstream>
#include <array>
#include <string>

#include "mpi.h"

#include "Kokkos_Core.hpp"

#include "../utils/utils.hh"
#include "../utils/array.hh"
#include "../utils/readers.hh"
#include "../utils/kokkos_utils.hh"

#include "CanopyHydrology.hh"
#include "CanopyHydrology_SnowWater_impl.hh"

#define UNIT_TEST 1

int main(int argc, char ** argv)
{
  // MPI_Init, etc
  MPI_Init(&argc,&argv);
  Kokkos::initialize( argc, argv);
  int myrank, n_procs;
  MPI_Comm_size(MPI_COMM_WORLD,&n_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    std::cout << "CanopyHydrology: Kokkos" << std::endl
              << "=======================" << std::endl
              << "Problem Setup" << std::endl
              << "--------------------" << std::endl
              << " n_procs = " << n_procs << std::endl;
  }

  // get ranks in x, y
  const auto domain_decom =
      ELM::Utils::get_domain_decomposition(n_procs, argc, argv);
  const int nx_procs = std::get<0>(domain_decom);
  const int ny_procs = std::get<1>(domain_decom);

  // NOTE: _global indicates values that are across all ranks.  The absence of
  // global means the variable is spatially local.
  const int start_year = 2014;
  const int start_month = 1;
  const int n_months = 3;
  const int n_pfts = 17;
  const int write_interval = 8 * 12;
  
  const std::string dir_atm = "ATM_DATA_LOCATION";
  const std::string dir_elm = "ELM_DATA_LOCATION";

  // dimension: time, lat (ny), lon (nx)
  const std::string basename1("Precip3Hrly/clmforc.GSWP3.c2011.0.5x0.5.Prec.");
  const auto problem_dims = 
      ELM::IO::get_dimensions(MPI_COMM_WORLD, dir_atm, basename1, start_year, start_month, n_months);
  const int n_times = std::get<0>(problem_dims);
  const int ny_global = std::get<1>(problem_dims);
  const int nx_global = std::get<2>(problem_dims);
  if (myrank == 0) {
    std::cout << " dimensions:" << std::endl
              << "   n_time = " << n_times << std::endl
              << "   n_lat = " << ny_global << std::endl
              << "   n_lon = " << nx_global << std::endl;
  }

  // domain decomposition
  assert(nx_global % nx_procs == 0 && "Currently expect perfectly divisible decomposition.");
  assert(ny_global % ny_procs == 0 && "Currently expect perfectly divisible decomposition.");

  // -- number of local grid cells per process
  const int nx_local = nx_global / nx_procs;
  const int ny_local = ny_global / ny_procs;
  const int n_grid_cells = nx_local * ny_local;
  if (myrank == 0) {
    std::cout << " domain decomposition = " << nx_procs << "," << ny_procs << std::endl
              << " local problem size = " << nx_local << "," << ny_local << std::endl;
  }

  // -- where am i on the process grid?
  const int i_proc = myrank % nx_procs;
  const int j_proc = myrank / nx_procs;

  // -- where do my local unknowns start globally
  int i_begin_global = i_proc * nx_local;
  int j_begin_global = j_proc * ny_local;

  // allocate storage and initialize phenology input data
  // -- allocate
  MPI_Barrier(MPI_COMM_WORLD);

  // allocate storage and initialize phenology input data
  // -- allocate on device
  Kokkos::View<double***> elai("elai", n_months, n_grid_cells, n_pfts);
  Kokkos::View<double***> esai("esai", n_months, n_grid_cells, n_pfts);
  {
    // -- host mirror
    auto h_elai = Kokkos::create_mirror_view(elai);
    auto h_esai = Kokkos::create_mirror_view(esai);

    const std::string basename("surfdata_360x720cru_simyr1850_c180216.nc");
    // -- read
    ELM::IO::read_and_reshape_phenology(MPI_COMM_WORLD, dir_elm, basename, "ELAI",
            start_year, start_month, 
            i_begin_global, j_begin_global, ny_local, nx_local, h_elai);
    if (myrank == 0) {
      std::cout << "File I/O" << std::endl
                << "--------------------" << std::endl
                << "  Phenology LAI read" << std::endl;
    }

    ELM::IO::read_and_reshape_phenology(MPI_COMM_WORLD, dir_elm, basename, "ESAI",
            start_year, start_month,
            i_begin_global, j_begin_global, ny_local, nx_local, h_esai);
    if (myrank == 0) {
      std::cout << "  Phenology SAI read" << std::endl;
    }

    // -- copy to device
    Kokkos::deep_copy(elai, h_elai);
    Kokkos::deep_copy(esai, h_esai);
  } // destroys host views, 3D arrays
  
  // allocate storage and initialize forcing input data
  //
  // NOTE, this may be too big, and we'll have to stage forcing data?
  //
  // -- allocate
  Kokkos::View<double**> forc_rain("forc_rain", n_times, n_grid_cells);
  Kokkos::View<double**> forc_snow("forc_snow", n_times, n_grid_cells);
  Kokkos::View<double**> forc_air_temp("forc_air_temp", n_times, n_grid_cells);
  Kokkos::View<double**> forc_irrig("forc_irrig", n_times, n_grid_cells);
  double qflx_floodg = 0.0;
  {
    // -- host views
    auto h_forc_rain = Kokkos::create_mirror_view(forc_rain);
    auto h_forc_snow = Kokkos::create_mirror_view(forc_snow);
    auto h_forc_air_temp = Kokkos::create_mirror_view(forc_air_temp);

    {
      // -- read
      std::string basename("Precip3Hrly/clmforc.GSWP3.c2011.0.5x0.5.Prec.");
      ELM::IO::read_and_reshape_forcing(MPI_COMM_WORLD, dir_atm, basename, "PRECIP",
              start_year, start_month, n_months,
              i_begin_global, j_begin_global, ny_local, nx_local, h_forc_rain);
      if (myrank == 0) std::cout << "  Forcing precip read" << std::endl;
      Kokkos::deep_copy(h_forc_snow, h_forc_rain);

      basename="TPHWL3Hrly/clmforc.GSWP3.c2011.0.5x0.5.TPQWL.";
      ELM::IO::read_and_reshape_forcing(MPI_COMM_WORLD, dir_atm, basename, "AIR_TEMP",
              start_year, start_month, n_months,
              i_begin_global, j_begin_global, ny_local, nx_local, h_forc_air_temp);
      if (myrank == 0) std::cout << "  Forcing air temperature read" << std::endl;
    }
    ELM::Utils::convert_precip_to_rain_snow(h_forc_rain,h_forc_snow,h_forc_air_temp);
    if (myrank == 0) std::cout << "  Converted precip to rain + snow" << std::endl;

    // -- copy to device
    Kokkos::deep_copy(forc_rain, h_forc_rain);
    Kokkos::deep_copy(forc_snow, h_forc_snow);
    Kokkos::deep_copy(forc_air_temp, h_forc_air_temp);
  } // destroys host views, Array2D objects

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    std::cout << "Test Execution" << std::endl
              << "--------------" << std::endl;
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
  // // need them to be for these kernels. --etc
  Kokkos::View<double**> z("z", n_grid_cells, n_levels_snow);
  Kokkos::View<double**> zi("zi", n_grid_cells, n_levels_snow);
  Kokkos::View<double**> dz("dz", n_grid_cells, n_levels_snow);

  // state variables that require ICs and evolve (in/out)
  Kokkos::View<double**> h2ocan("h2ocan", n_grid_cells, n_pfts);
  Kokkos::View<double**> swe_old("swe_old", n_grid_cells, n_levels_snow);
  Kokkos::View<double**> h2osoi_liq("h2osoi_liq", n_grid_cells, n_levels_snow);
  Kokkos::View<double**> h2osoi_ice("h2osoi_ice", n_grid_cells, n_levels_snow);
  Kokkos::View<double**> t_soisno("t_soisno", n_grid_cells, n_levels_snow);
  Kokkos::View<double**> frac_iceold("frac_iceold", n_grid_cells, n_levels_snow);
  
  Kokkos::View<double*> t_grnd("t_grnd", n_grid_cells);
  Kokkos::View<double*> h2osno("h2osno", n_grid_cells);
  Kokkos::View<double*> snow_depth("snow_depth", n_grid_cells);
  Kokkos::View<int*> snow_level("snow_level", n_grid_cells);

  Kokkos::View<double*> h2osfc("h2osfc", n_grid_cells);
  Kokkos::View<double*> frac_h2osfc("frac_h2osfc", n_grid_cells);
  
  // output fluxes by pft
  Kokkos::View<double**> qflx_prec_intr("qflx_prec_intr", n_grid_cells, n_pfts);
  Kokkos::View<double**> qflx_irrig("qflx_irrig", n_grid_cells, n_pfts );
  Kokkos::View<double**> qflx_prec_grnd("qflx_prec_grnd", n_grid_cells, n_pfts );
  Kokkos::View<double**> qflx_snwcp_liq("qflx_snwcp_liq", n_grid_cells, n_pfts);
  Kokkos::View<double**> qflx_snwcp_ice ("qflx_snwcp_ice ", n_grid_cells, n_pfts );
  Kokkos::View<double**> qflx_snow_grnd_patch("qflx_snow_grnd_patch", n_grid_cells, n_pfts );
  Kokkos::View<double**> qflx_rain_grnd("qflx_rain_grnd", n_grid_cells, n_pfts );

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  Kokkos::View<double*> integrated_snow("integrated_snow", n_grid_cells);
  
  // output fluxes, state by the column
  Kokkos::View<double*> qflx_snow_grnd_col("qflx_snow_grnd_col", n_grid_cells);
  Kokkos::View<double*> qflx_snow_h2osfc("qflx_snow_h2osfc", n_grid_cells);
  Kokkos::View<double*> qflx_h2osfc2topsoi("qflx_h2osfc2topsoi", n_grid_cells);
  Kokkos::View<double*> qflx_floodc("qflx_floodc", n_grid_cells);
  
  Kokkos::View<double*> frac_sno_eff("frac_sno_eff", n_grid_cells);
  Kokkos::View<double*> frac_sno("frac_sno", n_grid_cells);

#ifdef UNIT_TEST  
  // for unit testing
  std::ofstream soln_file;
  {
    auto min_max_sum_water = ELM::ELMKokkos::min_max_sum2(MPI_COMM_WORLD, h2ocan);
    auto min_max_sum_snow = ELM::ELMKokkos::min_max_sum1(MPI_COMM_WORLD, h2osno);
    auto min_max_sum_surfacewater = ELM::ELMKokkos::min_max_sum1(MPI_COMM_WORLD, frac_h2osfc);
    if (myrank == 0) std::cout << "  writing ts 0" << std::endl;

    if (myrank == 0) {
      soln_file.open("test_CanopyHydrology_module.soln");
      soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;

      soln_file << std::setprecision(16)
                << 0 << "\t" << min_max_sum_water[2] << "\t" << min_max_sum_water[0] << "\t" << min_max_sum_water[1]
                << "\t" << min_max_sum_snow[2] << "\t" << min_max_sum_snow[0] << "\t" << min_max_sum_snow[1]
                << "\t" << min_max_sum_surfacewater[2] << "\t" << min_max_sum_surfacewater[0] << "\t" << min_max_sum_surfacewater[1];
    }
  }
#endif

  auto start = ELM::Utils::Clock::time();

  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (int t = 0; t != n_times; ++t) {
    // data is 3hourly, so we have 8 data per day
    int i_month = ELM::Utils::month_from_day((int)(t/8), start_month) ;

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
                        elai(i_month, g,p), esai(i_month, g,p), dewmx, frac_veg_nosno,
                        h2ocan(g,p), n_irrig_steps_left,
                        qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                        qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                        qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p)); 

                double fwet = 0., fdry = 0.;
                ELM::CanopyHydrology_FracWet(frac_veg_nosno, h2ocan(g,p),
                        elai(i_month, g,p), esai(i_month, g,p),
                        dewmx, fwet, fdry);

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

#ifdef UNIT_TEST
    if (t % write_interval == 0) {
      auto min_max_sum_water = ELM::ELMKokkos::min_max_sum2(MPI_COMM_WORLD, h2ocan);
      auto min_max_sum_snow = ELM::ELMKokkos::min_max_sum1(MPI_COMM_WORLD, h2osno);
      auto min_max_sum_surfacewater = ELM::ELMKokkos::min_max_sum1(MPI_COMM_WORLD, frac_h2osfc);
      if (myrank == 0) std::cout << "  writing ts " << t << std::endl;

      if (myrank == 0) {
        soln_file << std::setprecision(16) << 0
                  << "\t" << min_max_sum_water[2] << "\t" << min_max_sum_water[0] << "\t" << min_max_sum_water[1]
                  << "\t" << min_max_sum_snow[2] << "\t" << min_max_sum_snow[0] << "\t" << min_max_sum_snow[1]
                  << "\t" << min_max_sum_surfacewater[2] << "\t" << min_max_sum_surfacewater[0] << "\t" << min_max_sum_surfacewater[1];
      }
    }
#endif
    
  } // end timestep loop


  auto stop = ELM::Utils::Clock::time();
  auto times = ELM::Utils::Clock::min_max_mean(MPI_COMM_WORLD, stop-start);
  if (myrank == 0) {
    std::cout << "Timing: min: "<< times[0] <<  ", max: " << times[1]
              << ", mean: " << times[2] << std::endl;
  }

#ifdef UNIT_TEST
  if (myrank == 0) soln_file.close();
#endif

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
