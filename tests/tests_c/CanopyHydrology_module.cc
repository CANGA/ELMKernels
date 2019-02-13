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

#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"


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
  
  // FIXME: confusing parameters that need to be fixed --etc
  const double n_melt = 0.7;
                               
  // phenology input
  ELM::Utils::MatrixStatePFT elai;
  ELM::Utils::MatrixStatePFT esai;
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, elai, esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, elai, esai);

  // forcing input
  ELM::Utils::MatrixForc forc_rain;
  ELM::Utils::MatrixForc forc_snow;
  ELM::Utils::MatrixForc forc_air_temp;
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, forc_rain, forc_snow, forc_air_temp);
  ELM::Utils::MatrixForc forc_irrig; forc_irrig = 0.;
  double qflx_floodg = 0.0;

  
  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // need them to be for these kernels. --etc
  auto z = ELM::Utils::MatrixStateSoilColumn(0.);
  auto zi = ELM::Utils::MatrixStateSoilColumn(0.);
  auto dz = ELM::Utils::MatrixStateSoilColumn(0.);

  // state variables that require ICs and evolve (in/out)
  auto h2ocan = ELM::Utils::MatrixStatePFT(0.);
  auto swe_old = ELM::Utils::MatrixStateSoilColumn(0.);
  auto h2osoi_liq = ELM::Utils::MatrixStateSoilColumn(0.);
  auto h2osoi_ice = ELM::Utils::MatrixStateSoilColumn(0.);
  auto t_soisno = ELM::Utils::MatrixStateSoilColumn(0.);
  auto frac_iceold = ELM::Utils::MatrixStateSoilColumn(0.);
  auto t_grnd = ELM::Utils::VectorColumn(0.);
  auto h2osno = ELM::Utils::VectorColumn(0.);
  auto snow_depth = ELM::Utils::VectorColumn(0.);
  auto snl = ELM::Utils::VectorColumnInt(0.); // note this tracks the snow_depth

  auto h2osfc = ELM::Utils::VectorColumn(0.);
  auto frac_h2osfc = ELM::Utils::VectorColumn(0.);

  
  // output fluxes by pft
  auto qflx_prec_intr = ELM::Utils::MatrixStatePFT();
  auto qflx_irrig = ELM::Utils::MatrixStatePFT();
  auto qflx_prec_grnd = ELM::Utils::MatrixStatePFT();
  auto qflx_snwcp_liq = ELM::Utils::MatrixStatePFT();
  auto qflx_snwcp_ice = ELM::Utils::MatrixStatePFT();
  auto qflx_snow_grnd_patch = ELM::Utils::MatrixStatePFT();
  auto qflx_rain_grnd = ELM::Utils::MatrixStatePFT();

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  auto int_snow = ELM::Utils::VectorColumn(0.);
  
  // output fluxes, state by the column
  auto qflx_snow_grnd_col = ELM::Utils::VectorColumn();
  auto qflx_snow_h2osfc = ELM::Utils::VectorColumn();
  auto qflx_h2osfc2topsoi = ELM::Utils::VectorColumn();
  auto qflx_floodc = ELM::Utils::VectorColumn();

  auto frac_sno_eff = ELM::Utils::VectorColumn();
  auto frac_sno = ELM::Utils::VectorColumn();
  

  // std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  // auto min_max = std::minmax_element(h2ocan.begin(), h2ocan.end());
  // std::cout << std::setprecision(16)
  //           << 0 << "\t" << std::accumulate(h2ocan.begin(), h2ocan.end(), 0.)
  //           << "\t" << *min_max.first
  //           << "\t" << *min_max.second << std::endl;
  
  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

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
                elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                h2ocan(g,p), n_irrig_steps_left,
                qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p));
        //printf("%i %i %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", g, p, forc_rain(t,g), forc_snow(t,g), elai(g,p), esai(g,p), h2ocan(g,p), qflx_prec_intr[g]);

        //
        // Calculate fraction of ?LAI? that is wet vs dry.
        //
        // FIXME: this currently punts on what to do with the fwet/fdry variables.
        // Surely they should be something, as such this is dead code.
        // --etc
        double fwet = 0., fdry = 0.;
        ELM::CanopyHydrology_FracWet(frac_veg_nosno, h2ocan(g,p), elai(g,p), esai(g,p), dewmx, fwet, fdry);
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
              snow_depth[g], h2osno[g], int_snow[g], swe_old[g],
              h2osoi_liq[g], h2osoi_ice[g], t_soisno[g], frac_iceold[g],
              snl[g], dz[g], z[g], zi[g], newnode,
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

    
    // auto min_max = std::minmax_element(h2ocan.begin(), h2ocan.end());
    // std::cout << std::setprecision(16)
    //           << t+1 << "\t" << std::accumulate(h2ocan.begin(), h2ocan.end(), 0.)
    //           << "\t" << *min_max.first
    //           << "\t" << *min_max.second << std::endl;

  } // end timestep loop
  return 0;
}
