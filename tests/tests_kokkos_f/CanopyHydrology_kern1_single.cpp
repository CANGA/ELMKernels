#include <netcdf.h>
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
#include <Kokkos_Core.hpp>

#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"  

namespace ELM {
namespace Utils {

static const int n_months = 12;
static const int n_pfts = 17;
using MatrixState = MatrixStatic<n_months, n_pfts>;

static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                            // day * half hour timestep
using MatrixForc = MatrixStatic<n_max_times,1>;

} // namespace
} // namespace

int main(int argc, char ** argv)
{
  // dimensions
  const int n_months = 12;
  const int n_pfts = 17;
  const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                       // day * half hour timestep	

  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  const double dtime = 1800.0;

  Kokkos::initialize( argc, argv );
  {

  // phenology state
  typedef Kokkos::View<double**>  ViewMatrixType;
  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  ViewMatrixType elai( "elai", n_months, n_pfts );
  ViewMatrixType esai( "esai", n_months, n_pfts );
  ViewMatrixType::HostMirror h_elai = Kokkos::create_mirror_view( elai );
  ViewMatrixType::HostMirror h_esai = Kokkos::create_mirror_view( esai );
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
    

  // forcing state
  // ELM::Utils::MatrixForc forc_rain;
  // ELM::Utils::MatrixForc forc_snow;
  // ELM::Utils::MatrixForc forc_air_temp;
  ViewMatrixType forc_rain( "forc_rain", n_max_times,1 );
  ViewMatrixType forc_snow( "forc_snow", n_max_times,1 );
  ViewMatrixType forc_air_temp( "forc_air_temp", n_max_times,1 );
  ViewMatrixType::HostMirror h_forc_rain = Kokkos::create_mirror_view( forc_rain );
  ViewMatrixType::HostMirror h_forc_snow = Kokkos::create_mirror_view( forc_snow );
  ViewMatrixType::HostMirror h_forc_air_temp = Kokkos::create_mirror_view( forc_air_temp );
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 6, 1, h_forc_rain, h_forc_snow, h_forc_air_temp);

  


  double h2ocan = 0.0;
  double qflx_prec_intr = 0.;
  double qflx_irrig = 0.;
  double qflx_prec_grnd = 0.;
  double qflx_snwcp_liq = 0.;
  double qflx_snwcp_ice = 0.;
  double qflx_snow_grnd_patch = 0.;
  double qflx_rain_grnd = 0.;
  double total_precip_loop = 0.;
  
  Kokkos::deep_copy( elai, h_elai);
  Kokkos::deep_copy( esai, h_esai);
  Kokkos::deep_copy( forc_rain, h_forc_rain);
  Kokkos::deep_copy( forc_snow, h_forc_snow);
  Kokkos::deep_copy( forc_air_temp, h_forc_air_temp);

  std::cout << "Timestep, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr, total_precip_loop" << std::endl;
  for(size_t itime = 0; itime < n_times; itime += 1) { //Kokkos::parallel_for(n_times, KOKKOS_LAMBDA (const int itime) { 
    // note this call puts all precip as rain for testing
      Kokkos::parallel_reduce(n_times, KOKKOS_LAMBDA ( const int itime, double &update ) {
      update += forc_rain(itime,0) + forc_snow(itime,0);
      }, total_precip_loop );
      
    double total_precip = forc_rain(itime,0) + forc_snow(itime,0); 
    ELM::CanopyHydrology_Interception(dtime, total_precip, 0., 0.,
            ltype, ctype, urbpoi, do_capsnow,
            h_elai(5,7), h_esai(5,7), dewmx, frac_veg_nosno,
            h2ocan, n_irrig_steps_left,
            qflx_prec_intr, qflx_irrig, qflx_prec_grnd,
            qflx_snwcp_liq, qflx_snwcp_ice,
            qflx_snow_grnd_patch, qflx_rain_grnd); 
		
    std::cout << std::setprecision(16) << itime+1 << "\t" << total_precip << "\t" << h2ocan<< "\t" << qflx_prec_grnd << "\t" << qflx_prec_intr << "\t" << total_precip_loop  << std::endl; 
  }//);
	}
  Kokkos::finalize();
  return 0;
}
