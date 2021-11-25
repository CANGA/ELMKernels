
#include "surface_radiation.h"
#include "elm_constants.h"
#include "landtype.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>

/*    

tests SurfaceRadiation kernels: 
SurfRadZeroFluxes()
SurfRadAbsorbed()
SurfRadLayers()
SurfRadReflected()

the following data comes from the files SurfaceRadiation_IN.txt and SurfaceRadiation_OUT.txt located in test/data

sabg_soil
sabg_snow
sabg
sabv
fsa
sabg_lyr[nlevsno+1]
snl
ftdd[numrad]
ftid[numrad]
ftii[numrad]
forc_solad[numrad]
forc_solai[numrad]
fabd[numrad]
fabi[numrad]
albsod[numrad]
albsoi[numrad]
albsnd_hst[numrad]
albsni_hst[numrad]
albgrd[numrad]
albgri[numrad]
trd[numrad]
tri[numrad]
snow_depth
flx_absdv[nlevsno+1]
flx_absdn[nlevsno+1]
flx_absiv[nlevsno+1]
flx_absin[nlevsno+1]
albd[numrad]
albi[numrad]
fsr
*/



using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files 
  const std::string data_dir("/Users/80x/Software/elm_kernels/test/data/");
  const std::string input_file = data_dir + "SurfaceRadiation_IN.txt";
  const std::string output_file = data_dir + "SurfaceRadiation_OUT.txt";

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
  auto sabg_soil = create<ArrayD1>("sabg_soil", n_grid_cells);
  auto sabg_snow = create<ArrayD1>("sabg_snow", n_grid_cells);
  auto sabg = create<ArrayD1>("sabg", n_grid_cells);
  auto sabv = create<ArrayD1>("sabv", n_grid_cells);
  auto fsa = create<ArrayD1>("fsa", n_grid_cells);
  auto snl = create<ArrayI1>("snl", n_grid_cells);
  auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells);
  auto fsr = create<ArrayD1>("fsr", n_grid_cells);
  auto sabg_lyr = create<ArrayD2>("sabg_lyr", n_grid_cells, ELM::nlevsno + 1);
  auto ftdd = create<ArrayD2>("ftdd", n_grid_cells, ELM::numrad);
  auto ftid = create<ArrayD2>("ftid", n_grid_cells, ELM::numrad);
  auto ftii = create<ArrayD2>("ftii", n_grid_cells, ELM::numrad);
  auto forc_solad = create<ArrayD2>("forc_solad", n_grid_cells, ELM::numrad);
  auto forc_solai = create<ArrayD2>("forc_solai", n_grid_cells, ELM::numrad);
  auto fabd = create<ArrayD2>("fabd", n_grid_cells, ELM::numrad);
  auto fabi = create<ArrayD2>("fabi", n_grid_cells, ELM::numrad);
  auto albsod = create<ArrayD2>("albsod", n_grid_cells, ELM::numrad);
  auto albsoi = create<ArrayD2>("albsoi", n_grid_cells, ELM::numrad);
  auto albsnd_hst = create<ArrayD2>("albsnd_hst", n_grid_cells, ELM::numrad);
  auto albsni_hst = create<ArrayD2>("albsni_hst", n_grid_cells, ELM::numrad);
  auto albgrd = create<ArrayD2>("albgrd", n_grid_cells, ELM::numrad);
  auto albgri = create<ArrayD2>("albgri", n_grid_cells, ELM::numrad);
  auto flx_absdv = create<ArrayD2>("flx_absdv", n_grid_cells, ELM::nlevsno + 1);
  auto flx_absdn = create<ArrayD2>("flx_absdn", n_grid_cells, ELM::nlevsno + 1);
  auto flx_absiv = create<ArrayD2>("flx_absiv", n_grid_cells, ELM::nlevsno + 1);
  auto flx_absin = create<ArrayD2>("flx_absin", n_grid_cells, ELM::nlevsno + 1);
  auto albd = create<ArrayD2>("albd", n_grid_cells, ELM::numrad);
  auto albi = create<ArrayD2>("albi", n_grid_cells, ELM::numrad);

  // arrays to compare output (trd & tri are double[numrad])
  auto trd_array = create<ArrayD2>("trd", n_grid_cells, ELM::numrad);
  auto tri_array = create<ArrayD2>("tri", n_grid_cells, ELM::numrad);

  // input and output utility class objects
  ELM::IO::ELMtestinput in(input_file);
  ELM::IO::ELMtestinput out(output_file);

  for (std::size_t t = 1; t < 49; ++t) {

    // get input and output state for time t
    in.getState(t);
    out.getState(t);

    // parse input state and assign to variables
    in.parseState(sabg_soil);
    in.parseState(sabg_snow);
    in.parseState(sabg);
    in.parseState(sabv);
    in.parseState(fsa);
    in.parseState(snl);
    in.parseState(snow_depth);
    in.parseState(fsr);
    in.parseState(sabg_lyr[idx]);
    in.parseState(ftdd[idx]);
    in.parseState(ftid[idx]);
    in.parseState(ftii[idx]);
    in.parseState(forc_solad[idx]);
    in.parseState(forc_solai[idx]);
    in.parseState(fabd[idx]);
    in.parseState(fabi[idx]);
    in.parseState(albsod[idx]);
    in.parseState(albsoi[idx]);
    in.parseState(albsnd_hst[idx]);
    in.parseState(albsni_hst[idx]);
    in.parseState(albgrd[idx]);
    in.parseState(albgri[idx]);
    in.parseState(flx_absdv[idx]);
    in.parseState(flx_absdn[idx]);
    in.parseState(flx_absiv[idx]);
    in.parseState(flx_absin[idx]);
    in.parseState(albd[idx]);
    in.parseState(albi[idx]);

    // local to these kernel calls
    double trd[ELM::numrad] = {0.0,0.0};
    double tri[ELM::numrad] = {0.0,0.0};


    // call SurfaceRadiation kernels
    ELM::surface_radiation::SurfRadZeroFluxes(Land, sabg_soil[idx], sabg_snow[idx], sabg[idx], sabv[idx], fsa[idx],
                          sabg_lyr[idx]);

    ELM::surface_radiation::SurfRadAbsorbed(Land, snl[idx], ftdd[idx], ftid[idx],
         ftii[idx], forc_solad[idx],
         forc_solai[idx], fabd[idx],
         fabi[idx], albsod[idx],
         albsoi[idx], albsnd_hst[idx],
         albsni_hst[idx], albgrd[idx],
         albgri[idx], sabv[idx], fsa[idx], sabg[idx], sabg_soil[idx], sabg_snow[idx],
         trd, tri);


    ELM::surface_radiation::SurfRadLayers(Land, snl[idx], sabg[idx], sabg_snow[idx], snow_depth[idx], flx_absdv[idx],
                       flx_absdn[idx], flx_absiv[idx],
                       flx_absin[idx], trd, tri, sabg_lyr[idx]);

    ELM::surface_radiation::SurfRadReflected(Land, albd[idx], albi[idx],
                          forc_solad[idx], forc_solai[idx],
                          fsr[idx]);

    // compare kernel output to ELM output state
    out.compareOutput(sabg_soil);
    out.compareOutput(sabg_snow);
    out.compareOutput(sabg);
    out.compareOutput(sabv);
    out.compareOutput(fsa);
    out.compareOutput(snl);
    out.compareOutput(snow_depth);
    out.compareOutput(fsr);
    out.compareOutput(sabg_lyr[idx]);
    out.compareOutput(ftdd[idx]);
    out.compareOutput(ftid[idx]);
    out.compareOutput(ftii[idx]);
    out.compareOutput(forc_solad[idx]);
    out.compareOutput(forc_solai[idx]);
    out.compareOutput(fabd[idx]);
    out.compareOutput(fabi[idx]);
    out.compareOutput(albsod[idx]);
    out.compareOutput(albsoi[idx]);
    out.compareOutput(albsnd_hst[idx]);
    out.compareOutput(albsni_hst[idx]);
    out.compareOutput(albgrd[idx]);
    out.compareOutput(albgri[idx]);
    out.compareOutput(flx_absdv[idx]);
    out.compareOutput(flx_absdn[idx]);
    out.compareOutput(flx_absiv[idx]);
    out.compareOutput(flx_absin[idx]);
    out.compareOutput(albd[idx]);
    out.compareOutput(albi[idx]);

    // put trd & tri into ELM::Array for ouput comparison
    for (std::size_t i = 0; i < ELM::numrad; ++i) {
        trd_array(idx,i) = trd[i];
        tri_array(idx,i) = tri[i];
    }
    out.compareOutput(trd_array[idx]);
    out.compareOutput(tri_array[idx]);
  }
  return 0;
}
