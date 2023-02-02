
#include "canopy_hydrology.h"
#include "elm_constants.h"
#include "land_data.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>



/*
tests canopy_hydrology kernels: 
interception()
ground_flux()
fraction_wet()
snow_init()
fraction_h2osfc()

will need to add Irrigation() later

the following data comes from the files CanopyHydrology_IN.txt and CanopyHydrology_OUT.txt located in test/data

 ! dtime - forc_t + fwet & fdry are scalar; dz - snw_rds are arrays

CanopyHydrology_IN.txt contains
dtime
oldfflag
frac_veg_nosno
dewmx
elai
esai
h2ocan
irrig_rate ! not used
n_irrig_steps_left ! not used
qflx_irrig
qflx_snwcp_liq
qflx_snwcp_ice
qflx_rain_grnd
qflx_snow_grnd
do_capsnow
t_grnd
qflx_snow_melt
n_melt
snow_depth
h2osno
int_snow
frac_sno_eff
frac_sno
snl
micro_sigma
h2osfc
frac_h2osfc
forc_rain
forc_snow
forc_t
dz
z
zi
swe_old
h2osoi_liq
h2osoi_ice
t_soisno
frac_iceold
snw_rds

CanopyHydrology_OUT.txt contains all the variables in CanopyHydrology_IN.txt and:
fwet
fdry
*/

using namespace ELM::ELMdims;

using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files
  const std::string data_dir = TEST_DATA_DIR;
  const std::string input_file = data_dir + "CanopyHydrology_IN.txt";
  const std::string output_file = data_dir + "CanopyHydrology_OUT.txt";

  // hardwired
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;

  // temporary data to pass between functions
  double qflx_candrip;
  double qflx_through_snow;
  double qflx_through_rain;
  double fracsnow;
  double fracrain;

  // ELM state variables
  auto oldfflag = create<ArrayI1>("oldfflag", n_grid_cells);
  auto dewmx = create<ArrayD1>("dewmx", n_grid_cells);
  auto qflx_irrig = create<ArrayD1>("qflx_irrig", n_grid_cells);
  auto do_capsnow = create<ArrayB1>("do_capsnow", n_grid_cells);
  auto n_melt = create<ArrayD1>("n_melt", n_grid_cells);
  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto forc_rain = create<ArrayD1>("forc_rain", n_grid_cells);
  auto forc_snow = create<ArrayD1>("forc_snow", n_grid_cells);
  auto elai = create<ArrayD1>("elai", n_grid_cells);
  auto esai = create<ArrayD1>("esai", n_grid_cells);
  auto h2ocan = create<ArrayD1>("h2ocan", n_grid_cells);
  //auto qflx_prec_grnd = create<ArrayD1>("qflx_prec_grnd", n_grid_cells);
  auto qflx_snwcp_liq = create<ArrayD1>("qflx_snwcp_liq", n_grid_cells);
  auto qflx_snwcp_ice = create<ArrayD1>("qflx_snwcp_ice", n_grid_cells);
  auto qflx_snow_grnd = create<ArrayD1>("qflx_snow_grnd", n_grid_cells);
  auto qflx_rain_grnd = create<ArrayD1>("qflx_rain_grnd", n_grid_cells);
  auto forc_t = create<ArrayD1>("forc_t", n_grid_cells);
  auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
  auto qflx_snow_melt = create<ArrayD1>("qflx_snow_melt", n_grid_cells);
  auto micro_sigma = create<ArrayD1>("micro_sigma", n_grid_cells);
  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, nlevsno() + nlevgrnd());
  auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells);
  auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);
  //auto qflx_snow_h2osfc = create<ArrayD1>("qflx_snow_h2osfc", n_grid_cells);
  auto h2osfc = create<ArrayD1>("h2osfc", n_grid_cells);
  auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
  auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", n_grid_cells);
  auto swe_old = create<ArrayD2>("swe_old", n_grid_cells, nlevsno());
  auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);
  auto snw_rds = create<ArrayD2>("snw_rds", n_grid_cells, nlevsno());
  auto int_snow = create<ArrayD1>("int_snow", n_grid_cells);
  auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, nlevsno() + nlevgrnd());
  auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, nlevsno() + nlevgrnd());
  auto snl = create<ArrayI1>("snl", n_grid_cells);
  auto frac_iceold = create<ArrayD2>("frac_iceold", n_grid_cells, nlevsno() + nlevgrnd());
  auto dz = create<ArrayD2>("dz", n_grid_cells, nlevsno() + nlevgrnd());
  auto z = create<ArrayD2>("z", n_grid_cells, nlevsno() + nlevgrnd());
  auto zi = create<ArrayD2>("zi", n_grid_cells, nlevsno() + nlevgrnd() + 1);
  auto fwet = create<ArrayD1>("fwet", n_grid_cells);
  auto fdry = create<ArrayD1>("fdry", n_grid_cells);

  // input and output utility class objects
  ELM::IO::ELMtestinput in(input_file);
  ELM::IO::ELMtestinput out(output_file);

  for (std::size_t t = 1; t < 49; ++t) {

    // get input and output state for time t
    in.getState(t);
    out.getState(t);

    // parse input state and assign to variables
    in.parseState(oldfflag);
    in.parseState(frac_veg_nosno);
    in.parseState(dewmx);
    in.parseState(elai);
    in.parseState(esai);
    in.parseState(h2ocan);
    in.parseState(qflx_irrig);
    in.parseState(qflx_snwcp_liq);
    in.parseState(qflx_snwcp_ice);
    in.parseState(qflx_rain_grnd);
    //in.parseState(qflx_prec_grnd);
    in.parseState(qflx_snow_grnd);
    in.parseState(do_capsnow);
    in.parseState(t_grnd);
    in.parseState(qflx_snow_melt);
    in.parseState(n_melt);
    in.parseState(snow_depth);
    in.parseState(h2osno);
    in.parseState(int_snow);
    //in.parseState(qflx_snow_h2osfc);
    in.parseState(frac_sno_eff);
    in.parseState(frac_sno);
    in.parseState(snl);
    in.parseState(micro_sigma);
    in.parseState(h2osfc);
    in.parseState(frac_h2osfc);
    in.parseState(forc_rain);
    in.parseState(forc_snow);
    in.parseState(forc_t);
    in.parseState(dz[idx]);
    in.parseState(z[idx]);
    in.parseState(zi[idx]);
    in.parseState(swe_old[idx]);
    in.parseState(h2osoi_liq[idx]);
    in.parseState(h2osoi_ice[idx]);
    in.parseState(t_soisno[idx]);
    in.parseState(frac_iceold[idx]);
    in.parseState(snw_rds[idx]);


    // call CanopyHydrology kernels
    ELM::canopy_hydrology::interception(Land, frac_veg_nosno[idx], forc_rain[idx], forc_snow[idx], dewmx[idx], elai[idx], esai[idx], dtime,
                        h2ocan[idx], qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain);

    ELM::canopy_hydrology::ground_flux(Land, do_capsnow[idx], frac_veg_nosno[idx], forc_rain[idx], forc_snow[idx], qflx_irrig[idx],
                      qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain,
                      qflx_snwcp_liq[idx], qflx_snwcp_ice[idx], qflx_snow_grnd[idx], qflx_rain_grnd[idx]);

    ELM::canopy_hydrology::fraction_wet(Land, frac_veg_nosno[idx], dewmx[idx], elai[idx], esai[idx], h2ocan[idx], fwet[idx], fdry[idx]);

    ELM::canopy_hydrology::snow_init(Land, dtime, do_capsnow[idx], oldfflag[idx], forc_t[idx], t_grnd[idx], qflx_snow_grnd[idx], qflx_snow_melt[idx],
                    n_melt[idx], snow_depth[idx], h2osno[idx], int_snow[idx], swe_old[idx], h2osoi_liq[idx],
                    h2osoi_ice[idx], t_soisno[idx], frac_iceold[idx], snl[idx], dz[idx], z[idx], zi[idx], snw_rds[idx],
                    frac_sno_eff[idx], frac_sno[idx]);

    ELM::canopy_hydrology::fraction_h2osfc(Land, micro_sigma[idx], h2osno[idx], h2osfc[idx], h2osoi_liq[idx], frac_sno[idx], frac_sno_eff[idx],
                    frac_h2osfc[idx]);

    // compare kernel output to ELM output state
    out.compareOutput(oldfflag);
    out.compareOutput(frac_veg_nosno);
    out.compareOutput(dewmx);
    out.compareOutput(elai);
    out.compareOutput(esai);
    out.compareOutput(h2ocan);
    out.compareOutput(qflx_irrig);
    out.compareOutput(qflx_snwcp_liq);
    out.compareOutput(qflx_snwcp_ice);
    out.compareOutput(qflx_rain_grnd);
    //out.compareOutput(qflx_prec_grnd);
    out.compareOutput(qflx_snow_grnd);
    out.compareOutput(do_capsnow);
    out.compareOutput(t_grnd);
    out.compareOutput(qflx_snow_melt);
    out.compareOutput(n_melt);
    out.compareOutput(snow_depth);
    out.compareOutput(h2osno);
    out.compareOutput(int_snow);
    //out.compareOutput(qflx_snow_h2osfc);
    out.compareOutput(frac_sno_eff);
    out.compareOutput(frac_sno);
    out.compareOutput(snl);
    out.compareOutput(micro_sigma);
    out.compareOutput(h2osfc);
    out.compareOutput(frac_h2osfc);
    out.compareOutput(forc_rain);
    out.compareOutput(forc_snow);
    out.compareOutput(forc_t);
    out.compareOutput(fwet);
    out.compareOutput(fdry);

    out.compareOutput(dz[idx]);
    out.compareOutput(z[idx]);
    out.compareOutput(zi[idx]);
    out.compareOutput(swe_old[idx]);
    out.compareOutput(h2osoi_liq[idx]);
    out.compareOutput(h2osoi_ice[idx]);
    out.compareOutput(t_soisno[idx]);
    out.compareOutput(frac_iceold[idx]);
    out.compareOutput(snw_rds[idx]);
  }
  return 0;
}
