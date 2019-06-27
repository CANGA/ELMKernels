#include <netcdf.h>
#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "utils.hh"
#include "readers.hh"

#include "CanopyHydrology.hh"

using namespace std::chrono;
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

  // phenology state
  ELM::Utils::MatrixState elai;
  ELM::Utils::MatrixState esai;
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, elai, esai);

  // forcing state
  ELM::Utils::MatrixForc forc_rain;
  ELM::Utils::MatrixForc forc_snow;
  ELM::Utils::MatrixForc forc_air_temp;
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 6, 1, forc_rain, forc_snow, forc_air_temp);

  double h2ocan = 0.0;
  double qflx_prec_intr = 0.;
  double qflx_irrig = 0.;
  double qflx_prec_grnd = 0.;
  double qflx_snwcp_liq = 0.;
  double qflx_snwcp_ice = 0.;
  double qflx_snow_grnd_patch = 0.;
  double qflx_rain_grnd = 0.;
  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_single.soln");
  soln_file << "Timestep, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr" << std::endl;
  auto start = high_resolution_clock::now();
  for(size_t itime = 0; itime < n_times; itime += 1) {
    // note this call puts all precip as rain for testing
    double total_precip = forc_rain[itime][0] + forc_snow[itime][0];
    ELM::CanopyHydrology_Interception(dtime, total_precip, 0., 0.,
            ltype, ctype, urbpoi, do_capsnow,
            elai[5][7], esai[5][7], dewmx, frac_veg_nosno,
            h2ocan, n_irrig_steps_left,
            qflx_prec_intr, qflx_irrig, qflx_prec_grnd,
            qflx_snwcp_liq, qflx_snwcp_ice,
            qflx_snow_grnd_patch, qflx_rain_grnd);
		
    soln_file << std::setprecision(16) << itime+1 << "\t" << total_precip << "\t" << h2ocan<< "\t" << qflx_prec_grnd << "\t" << qflx_prec_intr << std::endl;
  }
	auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: "<< duration.count() << " microseconds" << std::endl;
  return 0;
}
