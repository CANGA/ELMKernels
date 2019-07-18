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
// #include <mpi.h>
#include <chrono>
#include <Kokkos_Core.hpp>
#include "utils.hh"
#include "readers.hh"
#include "landunit_varcon.h"
#include "column_varcon.h" 
#include "CanopyHydrology.hh"
using namespace std::chrono; 

namespace ELM {
namespace Utils {

static const int n_months = 12;
static const int n_pfts = 17;
static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                            // day * half hour timestep
static const int n_grid_cells = 24;

// using MatrixState = MatrixStatic<n_grid_cells, n_pfts>;
// using MatrixForc = MatrixStatic<n_max_times,n_grid_cells>;


} // namespace
} // namespace


int main(int argc, char ** argv)
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;
  using Kokkos::parallel_for;
  using Kokkos::TeamPolicy;
  using Kokkos::TeamThreadRange;
  
  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  double dtime = 1800.0;

  // int myrank, numprocs;
  // double mytime;

  // MPI_Init(&argc,&argv);
  // MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  // MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  // MPI_Barrier(MPI_COMM_WORLD);

  Kokkos::initialize( argc, argv );
  {

  typedef Kokkos::View<double**>  ViewMatrixType;
  //typedef Kokkos::Cuda ExecSpace;
  //typedef Kokkos::Cuda MemSpace;
  //typedef Kokkos::RangePolicy<ExecSpace> range_policy;
 
  ViewMatrixType elai( "elai", n_grid_cells, n_pfts );
  ViewMatrixType esai( "esai", n_grid_cells, n_pfts );
  ViewMatrixType::HostMirror h_elai = Kokkos::create_mirror_view( elai );
  ViewMatrixType::HostMirror h_esai = Kokkos::create_mirror_view( esai );

  // phenology state
  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, h_elai, h_esai);

  // forcing state
  // ELM::Utils::MatrixForc forc_rain;
  // ELM::Utils::MatrixForc forc_snow;
  // ELM::Utils::MatrixForc forc_air_temp;

  ViewMatrixType forc_rain( "forc_rain", n_max_times,n_grid_cells );
  ViewMatrixType forc_snow( "forc_snow", n_max_times,n_grid_cells );
  ViewMatrixType forc_air_temp( "forc_air_temp", n_max_times,n_grid_cells );
  ViewMatrixType::HostMirror h_forc_rain = Kokkos::create_mirror_view( forc_rain );
  ViewMatrixType::HostMirror h_forc_snow = Kokkos::create_mirror_view( forc_snow );
  ViewMatrixType::HostMirror h_forc_air_temp = Kokkos::create_mirror_view( forc_air_temp );
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  ViewMatrixType forc_irrig( "forc_irrig", n_max_times,n_grid_cells );
  ViewMatrixType::HostMirror h_forc_irrig = Kokkos::create_mirror_view( forc_irrig );
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
  ViewMatrixType qflx_prec_intr( "qflx_prec_intr", n_grid_cells, n_pfts );
  ViewMatrixType qflx_irrig( "qflx_irrig", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_prec_grnd( "qflx_prec_grnd", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_snwcp_liq( "qflx_snwcp_liq", n_grid_cells, n_pfts );
  ViewMatrixType qflx_snwcp_ice ( "qflx_snwcp_ice ", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_snow_grnd_patch( "qflx_snow_grnd_patch", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_rain_grnd( "qflx_rain_grnd", n_grid_cells, n_pfts  );
  ViewMatrixType::HostMirror h_qflx_prec_intr = Kokkos::create_mirror_view( qflx_prec_intr );
  ViewMatrixType::HostMirror h_qflx_irrig = Kokkos::create_mirror_view( qflx_irrig);
  ViewMatrixType::HostMirror h_qflx_prec_grnd = Kokkos::create_mirror_view( qflx_prec_grnd);
  ViewMatrixType::HostMirror h_qflx_snwcp_liq = Kokkos::create_mirror_view(  qflx_snwcp_liq);
  ViewMatrixType::HostMirror h_qflx_snwcp_ice = Kokkos::create_mirror_view( qflx_snwcp_ice   );
  ViewMatrixType::HostMirror h_qflx_snow_grnd_patch = Kokkos::create_mirror_view( qflx_snow_grnd_patch );
  ViewMatrixType::HostMirror h_qflx_rain_grnd = Kokkos::create_mirror_view(  qflx_rain_grnd  );

  // output state by the pft
  ViewMatrixType h2o_can( "h2o_can", n_grid_cells, n_pfts );
  //auto h2o_can = ELM::Utils::MatrixState(); 
  ViewMatrixType::HostMirror h_h2o_can = Kokkos::create_mirror_view( h2o_can );
  //h_h2o_can = 0.;
  //auto h2o_can1 = ELM::Utils::MatrixState(); 
  // Array<int64_t, 2> a = h_h2o_can;
  //const size_t n0 = h2o_can.extent_0 ();

  //   const int64_t begin0 = h_h2o_can.begin();
  //   const int64_t end0= h_h2o_can.end();
  Kokkos::deep_copy( elai, h_elai);
  Kokkos::deep_copy( esai, h_esai);
  Kokkos::deep_copy( forc_rain, h_forc_rain);
  Kokkos::deep_copy( forc_snow, h_forc_snow);
  //Kokkos::deep_copy( forc_irrig, h_forc_irrig);
  Kokkos::deep_copy( forc_air_temp, h_forc_air_temp);
  //Kokkos::deep_copy( qflx_prec_intr, h_qflx_prec_intr);
  // Kokkos::deep_copy( qflx_irrig, h_qflx_irrig);
  // Kokkos::deep_copy( qflx_prec_grnd, h_qflx_prec_grnd);
  // Kokkos::deep_copy( qflx_snwcp_liq, h_qflx_snwcp_liq);
  // Kokkos::deep_copy( qflx_snwcp_ice, h_qflx_snwcp_ice);
  // Kokkos::deep_copy( qflx_snow_grnd_patch, h_qflx_snow_grnd_patch);
  // Kokkos::deep_copy( qflx_rain_grnd, h_qflx_rain_grnd);
  // Kokkos::deep_copy( h2o_can, h_h2o_can);
 
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



  Kokkos::Timer timer;
  auto start = high_resolution_clock::now();
  // mytime = MPI_Wtime();
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    typedef typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> > MDPolicyType_2D;

    // Construct 2D MDRangePolicy: lower and upper bounds provided, tile dims defaulted
    MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{n_grid_cells,n_pfts}} );

    // // Execute parallel_for with rank 2 MDRangePolicy
    // Kokkos::parallel_for( "md2d", mdpolicy_2d, ELM::CanopyHydrology_Interception(dtime,
    //               forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
    //               ltype, ctype, urbpoi, do_capsnow,
    //               elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
    //               h2o_can(g,p), n_irrig_steps_left,
    //               qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
    //               qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
    //              qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p)); );
    Kokkos::parallel_for("md2d",mdpolicy_2d,KOKKOS_LAMBDA (const size_t& g, const size_t& p) { 
                ELM::CanopyHydrology_Interception(dtime,
                  forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                  ltype, ctype, urbpoi, do_capsnow,
                  elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                  h2o_can(g,p), n_irrig_steps_left,
                  qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                  qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                 qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p)); });

        //   parallel_for (TeamPolicy<> (0,n_grid_cells), KOKKOS_LAMBDA (const size_t& g)
        //   {
        //   parallel_for (TeamThreadRange (g, n_pfts),
        //     [=] (size_t p) {
        //    ELM::CanopyHydrology_Interception(dtime,
        //           forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
        //           ltype, ctype, urbpoi, do_capsnow,
        //           elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
        //           h2o_can(g,p), n_irrig_steps_left,
        //           qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
        //           qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
        //          qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p));
        //   });
        // });

  Kokkos::deep_copy( h_qflx_irrig, qflx_irrig);
  Kokkos::deep_copy( h_qflx_prec_intr, qflx_prec_intr);
  Kokkos::deep_copy( h_qflx_prec_grnd, qflx_prec_grnd);
  Kokkos::deep_copy( h_qflx_snwcp_liq, qflx_snwcp_liq);
  Kokkos::deep_copy( h_qflx_snwcp_ice, qflx_snwcp_ice);
  Kokkos::deep_copy( h_qflx_snow_grnd_patch, qflx_snow_grnd_patch);
  Kokkos::deep_copy( h_qflx_rain_grnd, qflx_rain_grnd);
  Kokkos::deep_copy( h_h2o_can, h2o_can);

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

  

  double time = timer.seconds();
  double Gbytes = 1.0e-9 * double( sizeof(double) * ( n_grid_cells + n_grid_cells * n_pfts + n_pfts ) );

  printf( "  n_pfts( %d ) n_grid_cells( %d ) n_times ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
          n_pfts, n_grid_cells, n_times, Gbytes * 1000, time, Gbytes * n_times / time );
  // mytime = MPI_Wtime() - mytime;
  auto stop = high_resolution_clock::now();
  // std::cout <<"Timing from node "<< myrank  << " is "<< mytime << "seconds." << std::endl;
  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl; 
  }
  Kokkos::finalize();
  return 0;
  // MPI_Finalize();
}
