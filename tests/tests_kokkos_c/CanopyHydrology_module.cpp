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
#include <mpi.h>
#include <chrono>
#include <Kokkos_Core.hpp>
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

using MatrixStatePFT = MatrixStatic<n_grid_cells, n_pfts>;
using MatrixStateSoilColumn = MatrixStatic<n_grid_cells, n_levels_snow>;
using MatrixForc = MatrixStatic<n_max_times,n_grid_cells>;
using VectorColumn = VectorStatic<n_grid_cells>;
using VectorColumnInt = VectorStatic<n_grid_cells,int>;

} // namespace
} // namespace


int main(int argc, char ** argv)
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;
  using ELM::Utils::n_levels_snow;
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
  const double dtime = 1800.0;

  // fixed magic parameters for SnowWater
  const double qflx_snow_melt = 0.;

  // fixed magic parameters for fracH2Osfc  
  const int oldfflag = 0;
  const double micro_sigma = 0.1;
  const double min_h2osfc = 1.0e-8;
  const double n_melt = 0.7;

  int myrank, numprocs;
  double mytime;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Barrier(MPI_COMM_WORLD);
  
  Kokkos::initialize( argc, argv );
  {                             
  // phenology input
  typedef Kokkos::View<double*>   ViewVectorType;
  typedef Kokkos::View<double**>  ViewMatrixType;
  typedef Kokkos::View<int**>  ViewMatrixType1;
  typedef Kokkos::View<int*>   ViewVectorType1;
  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  ViewMatrixType elai( "elai", n_grid_cells, n_pfts );
  ViewMatrixType esai( "esai", n_grid_cells, n_pfts );
  ViewMatrixType::HostMirror h_elai = Kokkos::create_mirror_view( elai );
  ViewMatrixType::HostMirror h_esai = Kokkos::create_mirror_view( esai );
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, h_elai, h_esai);

  // forcing input
  ViewMatrixType forc_rain( "forc_rain", n_max_times,n_grid_cells );
  ViewMatrixType forc_snow( "forc_snow", n_max_times,n_grid_cells );
  ViewMatrixType forc_air_temp( "forc_air_temp", n_max_times,n_grid_cells );
  ViewMatrixType::HostMirror h_forc_rain = Kokkos::create_mirror_view( forc_rain );
  ViewMatrixType::HostMirror h_forc_snow = Kokkos::create_mirror_view( forc_snow );
  ViewMatrixType::HostMirror h_forc_air_temp = Kokkos::create_mirror_view( forc_air_temp );
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  //ELM::Utils::MatrixForc forc_irrig; forc_irrig = 0.;
  ViewMatrixType forc_irrig( "forc_irrig", n_max_times,n_grid_cells );
  ViewMatrixType::HostMirror h_forc_irrig = Kokkos::create_mirror_view( forc_irrig );
  double qflx_floodg = 0.0;

  
  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // // need them to be for these kernels. --etc
  // auto z = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto zi = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto dz = ELM::Utils::MatrixStateSoilColumn(0.);
  ViewMatrixType z( "z", n_grid_cells, n_levels_snow );
  ViewMatrixType zi( "zi", n_grid_cells, n_levels_snow );
  ViewMatrixType dz( "dz", n_grid_cells, n_levels_snow );
  ViewMatrixType::HostMirror h_z = Kokkos::create_mirror_view( z );
  ViewMatrixType::HostMirror h_zi = Kokkos::create_mirror_view( zi );
  ViewMatrixType::HostMirror h_dz = Kokkos::create_mirror_view( dz );

  // state variables that require ICs and evolve (in/out)
  // auto h2ocan = ELM::Utils::MatrixStatePFT(0.);
  // auto swe_old = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto h2osoi_liq = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto h2osoi_ice = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto t_soisno = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto frac_iceold = ELM::Utils::MatrixStateSoilColumn(0.);
  ViewMatrixType h2ocan( "h2ocan", n_grid_cells, n_pfts );
  ViewMatrixType swe_old( "swe_old", n_grid_cells, n_levels_snow );
  ViewMatrixType h2osoi_liq( "h2osoi_liq", n_grid_cells, n_levels_snow );
  ViewMatrixType h2osoi_ice( "h2osoi_ice", n_grid_cells, n_levels_snow );
  ViewMatrixType t_soisno( "t_soisno", n_grid_cells, n_levels_snow );
  ViewMatrixType frac_iceold( "frac_iceold", n_grid_cells, n_levels_snow );
  ViewMatrixType::HostMirror h_h2ocan = Kokkos::create_mirror_view( h2ocan );
  ViewMatrixType::HostMirror h_swe_old = Kokkos::create_mirror_view( swe_old );
  ViewMatrixType::HostMirror h_h2osoi_liq = Kokkos::create_mirror_view( h2osoi_liq );
  ViewMatrixType::HostMirror h_h2osoi_ice = Kokkos::create_mirror_view( h2osoi_ice );
  ViewMatrixType::HostMirror h_t_soisno = Kokkos::create_mirror_view( t_soisno );
  ViewMatrixType::HostMirror h_frac_iceold = Kokkos::create_mirror_view( frac_iceold );

  // auto t_grnd = ELM::Utils::VectorColumn(0.);
  // auto h2osno = ELM::Utils::VectorColumn(0.);
  // auto snow_depth = ELM::Utils::VectorColumn(0.);
  // auto snl = ELM::Utils::VectorColumnInt(0.); // note this tracks the snow_depth
  ViewVectorType t_grnd( "t_grnd", n_grid_cells );
  ViewVectorType h2osno( "h2osno", n_grid_cells );
  ViewVectorType snow_depth( "snow_depth", n_grid_cells );
  ViewVectorType1 snow_level( "snow_level", n_grid_cells );
  ViewVectorType::HostMirror h_t_grnd = Kokkos::create_mirror_view(  t_grnd);
  ViewVectorType::HostMirror h_h2osno = Kokkos::create_mirror_view( h2osno);
  ViewVectorType::HostMirror h_snow_depth = Kokkos::create_mirror_view(  snow_depth);
  ViewVectorType1::HostMirror h_snow_level = Kokkos::create_mirror_view( snow_level);

  // auto h2osfc = ELM::Utils::VectorColumn(0.);
  // auto frac_h2osfc = ELM::Utils::VectorColumn(0.);
  ViewVectorType h2osfc( "h2osfc", n_grid_cells );
  ViewVectorType frac_h2osfc( "frac_h2osfc", n_grid_cells );
  ViewVectorType::HostMirror h_h2osfc = Kokkos::create_mirror_view(  h2osfc);
  ViewVectorType::HostMirror h_frac_h2osfc = Kokkos::create_mirror_view( frac_h2osfc);

  
  // output fluxes by pft
  // auto qflx_prec_intr = ELM::Utils::MatrixStatePFT();
  // auto qflx_irrig = ELM::Utils::MatrixStatePFT();
  // auto qflx_prec_grnd = ELM::Utils::MatrixStatePFT();
  // auto qflx_snwcp_liq = ELM::Utils::MatrixStatePFT();
  // auto qflx_snwcp_ice = ELM::Utils::MatrixStatePFT();
  // auto qflx_snow_grnd_patch = ELM::Utils::MatrixStatePFT();
  // auto qflx_rain_grnd = ELM::Utils::MatrixStatePFT();
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

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  //auto integrated_snow = ELM::Utils::VectorColumn(0.);
  ViewVectorType integrated_snow( "integrated_snow", n_grid_cells );
  ViewVectorType::HostMirror h_integrated_snow = Kokkos::create_mirror_view(  integrated_snow);
  
  // output fluxes, state by the column
  // auto qflx_snow_grnd_col = ELM::Utils::VectorColumn();
  // auto qflx_snow_h2osfc = ELM::Utils::VectorColumn();
  // auto qflx_h2osfc2topsoi = ELM::Utils::VectorColumn();
  // auto qflx_floodc = ELM::Utils::VectorColumn();
  ViewVectorType qflx_snow_grnd_col( "qflx_snow_grnd_col", n_grid_cells );
  ViewVectorType qflx_snow_h2osfc( "qflx_snow_h2osfc", n_grid_cells );
  ViewVectorType qflx_h2osfc2topsoi( "qflx_h2osfc2topsoi", n_grid_cells );
  ViewVectorType qflx_floodc( "qflx_floodc", n_grid_cells );
  ViewVectorType::HostMirror h_qflx_snow_grnd_col = Kokkos::create_mirror_view(  qflx_snow_grnd_col);
  ViewVectorType::HostMirror h_qflx_snow_h2osfc = Kokkos::create_mirror_view( qflx_snow_h2osfc);
  ViewVectorType::HostMirror h_qflx_h2osfc2topsoi = Kokkos::create_mirror_view(  qflx_h2osfc2topsoi);
  ViewVectorType::HostMirror h_qflx_floodc = Kokkos::create_mirror_view( qflx_floodc);

  // auto frac_sno_eff = ELM::Utils::VectorColumn();
  // auto frac_sno = ELM::Utils::VectorColumn();
  ViewVectorType frac_sno_eff( "frac_sno_eff", n_grid_cells );
  ViewVectorType frac_sno( "frac_sno", n_grid_cells );
  ViewVectorType::HostMirror h_frac_sno_eff = Kokkos::create_mirror_view(  frac_sno_eff);
  ViewVectorType::HostMirror h_frac_sno = Kokkos::create_mirror_view( frac_sno);
  

  // std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  // auto min_max = std::minmax_element(h2ocan.begin(), h2ocan.end());
  // std::cout << std::setprecision(16)
  //           << 0 << "\t" << std::accumulate(h2ocan.begin(), h2ocan.end(), 0.)
  //           << "\t" << *min_max.first
  //           << "\t" << *min_max.second << std::endl;
  
  Kokkos::deep_copy( elai, h_elai);
  Kokkos::deep_copy( esai, h_esai);
  Kokkos::deep_copy( forc_rain, h_forc_rain);
  Kokkos::deep_copy( forc_snow, h_forc_snow);
  Kokkos::deep_copy( forc_air_temp, h_forc_air_temp);
  Kokkos::deep_copy( forc_irrig, h_forc_irrig);
  Kokkos::deep_copy( z, h_z);
  Kokkos::deep_copy( zi, h_zi);
  Kokkos::deep_copy( dz, h_dz);
  Kokkos::deep_copy( h2ocan, h_h2ocan);
  Kokkos::deep_copy( swe_old, h_swe_old);
  Kokkos::deep_copy( h2osoi_liq, h_h2osoi_liq);
  Kokkos::deep_copy( h2osoi_ice, h_h2osoi_ice);
  Kokkos::deep_copy( t_soisno, h_t_soisno);
  Kokkos::deep_copy( frac_iceold, h_frac_iceold);
  Kokkos::deep_copy( t_grnd, h_t_grnd);
  Kokkos::deep_copy( h2osno, h_h2osno);
  Kokkos::deep_copy( snow_depth, h_snow_depth);
  Kokkos::deep_copy( snow_level, h_snow_level);
  Kokkos::deep_copy( h2osfc, h_h2osfc);
  Kokkos::deep_copy( frac_h2osfc, h_frac_h2osfc);
  Kokkos::deep_copy( qflx_prec_intr,h_qflx_prec_intr);
  Kokkos::deep_copy( qflx_irrig,h_qflx_irrig);
  Kokkos::deep_copy( qflx_prec_grnd,h_qflx_prec_grnd);
  Kokkos::deep_copy( qflx_snwcp_liq,h_qflx_snwcp_liq);
  Kokkos::deep_copy( qflx_snwcp_ice,h_qflx_snwcp_ice);
  Kokkos::deep_copy( qflx_snow_grnd_patch,h_qflx_snow_grnd_patch);
  Kokkos::deep_copy( qflx_rain_grnd,h_qflx_rain_grnd);
  Kokkos::deep_copy( integrated_snow,h_integrated_snow);
  Kokkos::deep_copy( qflx_snow_grnd_col, h_qflx_snow_grnd_col);
  Kokkos::deep_copy( qflx_snow_h2osfc, h_qflx_snow_h2osfc);
  Kokkos::deep_copy( qflx_h2osfc2topsoi, h_qflx_h2osfc2topsoi);
  Kokkos::deep_copy( qflx_floodc, h_qflx_floodc);
  Kokkos::deep_copy( frac_sno_eff, h_frac_sno_eff);
  Kokkos::deep_copy( frac_sno, h_frac_sno);

  double* end1 = &h_h2ocan(n_grid_cells-1, n_pfts-1) ;
  double* end2 = &h_h2osno(n_grid_cells-1) ;
  double* end3 = &h_frac_h2osfc(n_grid_cells-1) ;
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

   Kokkos::Timer timer;
   auto start = high_resolution_clock::now();
   mytime = MPI_Wtime();
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    typedef typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> > MDPolicyType_2D;
        
    // Construct 2D MDRangePolicy: lower and upper bounds provided, tile dims defaulted
    MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{n_grid_cells,n_pfts}} );
          
    Kokkos::parallel_for("md2d",mdpolicy_2d,KOKKOS_LAMBDA (const size_t& g, const size_t& p) { 
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
    });
      
      // Column level operations
      // NOTE: this is effectively an accumulation kernel/task! --etc
    typedef Kokkos::TeamPolicy<>              team_policy ;
    typedef typename team_policy::member_type team_type ;
    Kokkos::parallel_for (Kokkos::TeamPolicy<> (n_grid_cells, n_pfts),
                     KOKKOS_LAMBDA (const team_type& team) {
      
      double sum = 0;
      Kokkos::parallel_reduce (Kokkos::TeamThreadRange (team, team.team_size()),
        [=] (const size_t& p, double& lsum) {
        lsum += qflx_snow_grnd_patch(team.league_rank(),p);
        }, sum);
      qflx_snow_grnd_col(team.league_rank()) = sum ;
    });
    // Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space>(0,n_grid_cells), KOKKOS_LAMBDA (const size_t& g, double& upd) {
    // upd += qflx_snow_grnd_patch(g,p);
    // }, sum);
    // qflx_snow_grnd_col(g) = sum ;

         
      // Calculate ?water balance? on the snow column, adding throughfall,
      // removing melt, etc.
      //
      // local outputs

    Kokkos::parallel_for (n_grid_cells,
                     KOKKOS_LAMBDA (const size_t& g) {
      int newnode;
      ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
              ltype, ctype, urbpoi, do_capsnow, oldfflag,
              forc_air_temp(t,g), t_grnd(g),
              qflx_snow_grnd_col(g), qflx_snow_melt, n_melt, frac_h2osfc(g),
              snow_depth(g), h2osno(g), integrated_snow(g), Kokkos::subview(swe_old, g , Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, g , Kokkos::ALL), Kokkos::subview(h2osoi_ice, g , Kokkos::ALL), Kokkos::subview(t_soisno, g , Kokkos::ALL), Kokkos::subview(frac_iceold, g , Kokkos::ALL),
              snow_level(g), Kokkos::subview(dz, g , Kokkos::ALL), Kokkos::subview(z, g , Kokkos::ALL), Kokkos::subview(zi, g , Kokkos::ALL), newnode,
              qflx_floodc(g), qflx_snow_h2osfc(g), frac_sno_eff(g), frac_sno(g));
      
      // Calculate Fraction of Water to the Surface?
      //
      // FIXME: Fortran black magic... h2osoi_liq is a vector, but the
      // interface specifies a single double.  For now passing the 0th
      // entry. --etc
       ELM::CanopyHydrology_FracH2OSfc(dtime, min_h2osfc, ltype, micro_sigma,
              h2osno(g), h2osfc(g), h2osoi_liq(g,0), frac_sno(g), frac_sno_eff(g),
              qflx_h2osfc2topsoi(g), frac_h2osfc(g));
      
    }); // end grid cell loop

    
    // auto min_max = std::minmax_element(h2ocan.begin(), h2ocan.end());
    // std::cout << std::setprecision(16)
    //           << t+1 << "\t" << std::accumulate(h2ocan.begin(), h2ocan.end(), 0.)
    //           << "\t" << *min_max.first
    //           << "\t" << *min_max.second << std::endl;
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

  double time = timer.seconds();
  double Gbytes = 1.0e-9 * double( sizeof(double) * ( n_grid_cells + n_grid_cells * n_pfts + n_pfts ) );

  printf( "  n_pfts( %d ) n_grid_cells( %d ) n_times ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
          n_pfts, n_grid_cells, n_times, Gbytes * 1000, time, Gbytes * n_times / time );

  mytime = MPI_Wtime() - mytime;
  auto stop = high_resolution_clock::now();
  std::cout <<"Timing from node "<< myrank  << " is "<< mytime << "seconds." << std::endl;
  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl; 
  }
  Kokkos::finalize();
  return 0;
  MPI_Finalize();
}
