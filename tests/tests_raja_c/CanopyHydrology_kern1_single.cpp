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

using namespace RAJA;
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

int main(int RAJA_UNUSED_ARG(argc), char **RAJA_UNUSED_ARG(argv[]))
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

  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  double* elai = memoryManager::allocate<double>(n_months * n_pfts);
  double* esai = memoryManager::allocate<double>(n_months * n_pfts);

  const int DIM = 2;
  RAJA::View<double, RAJA::Layout<DIM> > h_elai( elai, n_months, n_pfts );
  RAJA::View<double, RAJA::Layout<DIM> > h_esai( esai, n_months, n_pfts );
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
    

  // forcing state
  // ELM::Utils::MatrixForc forc_rain;
  // ELM::Utils::MatrixForc forc_snow;
  // ELM::Utils::MatrixForc forc_air_temp;
  double* forc_rain = memoryManager::allocate<double>( n_max_times * 1 );
  double* forc_snow = memoryManager::allocate<double>( n_max_times * 1 );
  double* forc_air_temp = memoryManager::allocate<double>( n_max_times * 1 );

  RAJA::View<double, RAJA::Layout<DIM> > h_forc_rain( forc_rain, n_max_times,1 );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_snow( forc_snow, n_max_times,1 );
  RAJA::View<double, RAJA::Layout<DIM> > h_forc_air_temp( forc_air_temp, n_max_times,1 );
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 6, 1, h_forc_rain, h_forc_snow, h_forc_air_temp);

  


  double h2ocan = 0.0;
  double qflx_prec_intr = 0.;
  double qflx_irrig = 0.;
  double qflx_prec_grnd = 0.;
  double qflx_snwcp_liq = 0.;
  double qflx_snwcp_ice = 0.;
  double qflx_snow_grnd_patch = 0.;
  double qflx_rain_grnd = 0.;
  double total_precip = 0.;
  
  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_single.soln");
  std::cout << "Timestep, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr, total_precip_loop" << std::endl;
  soln_file << "Timestep, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr, total_precip_loop" << std::endl;
  for(size_t itime = 0; itime < n_times; itime += 1) { 
    // note this call puts all precip as rain for testing
      
    total_precip = h_forc_rain(itime,0) + h_forc_snow(itime,0); 
    ELM::CanopyHydrology_Interception(dtime, total_precip, 0., 0.,
            ltype, ctype, urbpoi, do_capsnow,
            h_elai(5,7), h_esai(5,7), dewmx, frac_veg_nosno,
            h2ocan, n_irrig_steps_left,
            qflx_prec_intr, qflx_irrig, qflx_prec_grnd,
            qflx_snwcp_liq, qflx_snwcp_ice,
            qflx_snow_grnd_patch, qflx_rain_grnd); 
		
    soln_file << std::setprecision(16) << itime+1 << "\t" << total_precip << "\t" << h2ocan<< "\t" << qflx_prec_grnd << "\t" << qflx_prec_intr << std::endl; 
    std::cout << std::setprecision(16) << itime+1 << "\t" << total_precip << "\t" << h2ocan<< "\t" << qflx_prec_grnd << "\t" << qflx_prec_intr << std::endl; 
  }soln_file.close();

	memoryManager::deallocate(elai);
  memoryManager::deallocate(esai);
  memoryManager::deallocate(forc_rain);
  memoryManager::deallocate(forc_snow);
  memoryManager::deallocate(forc_air_temp);
  return 0;
}
