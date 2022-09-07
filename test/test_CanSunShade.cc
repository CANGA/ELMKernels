
#include "surface_radiation.h"
#include "elm_constants.h"
#include "land_data.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>

/*  Tests canopy_sunshade_fractions kernel
this kernel resided in SurfaceRadiation, but is called directly from the driver  

the following data comes from the files CanopySunShadeFractions_IN.txt and CanopySunShadeFractions_OUT.txt located in test/data

scalars:
nrad
elai
laisun
laisha

arrays:
tlai_z
fsun_z
forc_solad
forc_solai
fabd_sun_z
fabd_sha_z
fabi_sun_z
fabi_sha_z
parsun_z
parsha_z
laisun_z
laisha_z

*/

using namespace ELM::ELMdims;

using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files
  const std::string data_dir = TEST_DATA_DIR;
  const std::string input_file = data_dir + "CanopySunShadeFractions_IN.txt";
  const std::string output_file = data_dir + "CanopySunShadeFractions_OUT.txt";

  // hardwired
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;


  // ELM state variables
  auto nrad = create<ArrayI1>("nrad", n_grid_cells);
  auto elai = create<ArrayD1>("elai", n_grid_cells);
  auto laisun = create<ArrayD1>("laisun", n_grid_cells);
  auto laisha = create<ArrayD1>("laisha", n_grid_cells);
  auto tlai_z = create<ArrayD2>("tlai_z", n_grid_cells, nlevcan);
  auto fsun_z = create<ArrayD2>("fsun_z", n_grid_cells, nlevcan);
  auto forc_solad = create<ArrayD2>("forc_solad", n_grid_cells, numrad);
  auto forc_solai = create<ArrayD2>("forc_solai", n_grid_cells, numrad);
  auto fabd_sun_z = create<ArrayD2>("fabd_sun_z", n_grid_cells, nlevcan);
  auto fabd_sha_z = create<ArrayD2>("fabd_sha_z", n_grid_cells, nlevcan);
  auto fabi_sun_z = create<ArrayD2>("fabi_sun_z", n_grid_cells, nlevcan);
  auto fabi_sha_z = create<ArrayD2>("fabi_sha_z", n_grid_cells, nlevcan);
  auto parsun_z = create<ArrayD2>("parsun_z", n_grid_cells, nlevcan);
  auto parsha_z = create<ArrayD2>("parsha_z", n_grid_cells, nlevcan);
  auto laisun_z = create<ArrayD2>("laisun_z", n_grid_cells, nlevcan);
  auto laisha_z = create<ArrayD2>("laisha_z", n_grid_cells, nlevcan);


  // input and output utility class objects
  ELM::IO::ELMtestinput in(input_file);
  ELM::IO::ELMtestinput out(output_file);

  for (std::size_t t = 1; t < 49; ++t) {

    // get input and output state for time t
    in.getState(t);
    out.getState(t);

    // parse input state and assign to variables
    in.parseState(nrad);
    in.parseState(elai);
    in.parseState(laisun);
    in.parseState(laisha);
    in.parseState(tlai_z[idx]);
    in.parseState(fsun_z[idx]);
    in.parseState(forc_solad[idx]);
    in.parseState(forc_solai[idx]);
    in.parseState(fabd_sun_z[idx]);
    in.parseState(fabd_sha_z[idx]);
    in.parseState(fabi_sun_z[idx]);
    in.parseState(fabi_sha_z[idx]);
    in.parseState(parsun_z[idx]);
    in.parseState(parsha_z[idx]);
    in.parseState(laisun_z[idx]);
    in.parseState(laisha_z[idx]);

    // call canopy_sunshade_fractions kernel
    ELM::surface_radiation::canopy_sunshade_fractions(
        Land, nrad[idx], elai[idx], tlai_z[idx], fsun_z[idx],
        forc_solad[idx], forc_solai[idx],
        fabd_sun_z[idx], fabd_sha_z[idx],
        fabi_sun_z[idx], fabi_sha_z[idx],
        parsun_z[idx], parsha_z[idx],
        laisun_z[idx], laisha_z[idx], laisun[idx], laisha[idx]);

    // compare kernel output to ELM output state
    out.compareOutput(nrad);
    out.compareOutput(elai);
    out.compareOutput(laisun);
    out.compareOutput(laisha);
    out.compareOutput(tlai_z[idx]);
    out.compareOutput(fsun_z[idx]);
    out.compareOutput(forc_solad[idx]);
    out.compareOutput(forc_solai[idx]);
    out.compareOutput(fabd_sun_z[idx]);
    out.compareOutput(fabd_sha_z[idx]);
    out.compareOutput(fabi_sun_z[idx]);
    out.compareOutput(fabi_sha_z[idx]);
    out.compareOutput(parsun_z[idx]);
    out.compareOutput(parsha_z[idx]);
    out.compareOutput(laisun_z[idx]);
    out.compareOutput(laisha_z[idx]);

  }
  return 0;
}
