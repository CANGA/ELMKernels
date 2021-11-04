
#include "SurfaceAlbedo.h"
#include "SnowSNICAR.h"
#include "ELMConstants.h"
#include "LandType.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>

/*    
need to reconcile rhol/rhos taul/taus dimensionality
ReadPFT has 1d var[numpft] and these kernels use 2d var[numrad][numpft]

tests SurfaceAlbedo kernels: 
SurfRadZeroFluxes()
SurfRadAbsorbed()
SurfRadLayers()
SurfRadReflected()

the following data comes from the files SurfaceAlbedo_IN.txt and SurfaceAlbedo_OUT.txt located in test/data


urbpoi
coszen
snl

elai
mss_cnc_bcphi
mss_cnc_bcpho
mss_cnc_dst1
mss_cnc_dst2
mss_cnc_dst3
mss_cnc_dst4

vcmaxcintsun
vcmaxcintsha
albsod
albsoi
albgrd
albgri
albd
albi
fabd
fabd_sun
fabd_sha
fabi
fabi_sun
fabi_sha
ftdd
ftid
ftii
flx_absdv
flx_absdn
flx_absiv
flx_absin
mss_cnc_aer_in_fdb
frac_sno
albsnd
albsni

flx_absd_snw
flx_absi_snw

esai
tlai
tsai
nrad
ncan
tlai_z
tsai_z
fsun_z
fabd_sun_z
fabd_sha_z
fabi_sun_z
fabi_sha_z

t_veg
fwet
xl
rhol
rhos
taul
taus
h2osoi_vol

t_grnd
albsat
albdry
*/

using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;
using ArrayD3 = ELM::Array<double, 3>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(name, D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files 
  const std::string data_dir("/Users/80x/Software/elm_kernels/test/data/");
  const std::string input_file = data_dir + "SurfaceAlbedo_IN.txt";
  const std::string output_file = data_dir + "SurfaceAlbedo_OUT.txt";

  // hardwired
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;


// ELM state variables



//auto urbpoi = create<ArrayI1>("urbpoi", n_grid_cells);
auto coszen = create<ArrayD1>("coszen", n_grid_cells);
auto elai = create<ArrayD1>("elai", n_grid_cells);
auto snl = create<ArrayI1>("snl", n_grid_cells);


auto mss_cnc_bcphi = create<ArrayD2>("mss_cnc_bcphi", n_grid_cells, ELM::nlevsno);
auto mss_cnc_bcpho = create<ArrayD2>("mss_cnc_bcpho", n_grid_cells, ELM::nlevsno);
auto mss_cnc_dst1 = create<ArrayD2>("mss_cnc_dst1", n_grid_cells, ELM::nlevsno);
auto mss_cnc_dst2 = create<ArrayD2>("mss_cnc_dst2", n_grid_cells, ELM::nlevsno);
auto mss_cnc_dst3 = create<ArrayD2>("mss_cnc_dst3", n_grid_cells, ELM::nlevsno);
auto mss_cnc_dst4 = create<ArrayD2>("mss_cnc_dst4", n_grid_cells, ELM::nlevsno);



auto vcmaxcintsun = create<ArrayD1>("vcmaxcintsun", n_grid_cells);
auto vcmaxcintsha = create<ArrayD1>("vcmaxcintsha", n_grid_cells);
auto albsod = create<ArrayD2>("albsod", n_grid_cells, ELM::numrad);
auto albsoi = create<ArrayD2>("albsoi", n_grid_cells, ELM::numrad);
auto albgrd = create<ArrayD2>("albgrd", n_grid_cells, ELM::numrad);
auto albgri = create<ArrayD2>("albgri", n_grid_cells, ELM::numrad);
auto albd = create<ArrayD2>("albd", n_grid_cells, ELM::numrad);
auto albi = create<ArrayD2>("albi", n_grid_cells, ELM::numrad);
auto fabd = create<ArrayD2>("fabd", n_grid_cells, ELM::numrad);
auto fabd_sun = create<ArrayD2>("fabd_sun", n_grid_cells, ELM::numrad);
auto fabd_sha = create<ArrayD2>("fabd_sha", n_grid_cells, ELM::numrad);
auto fabi = create<ArrayD2>("fabi", n_grid_cells, ELM::numrad);
auto fabi_sun = create<ArrayD2>("fabi_sun", n_grid_cells, ELM::numrad);
auto fabi_sha = create<ArrayD2>("fabi_sha", n_grid_cells, ELM::numrad);
auto ftdd = create<ArrayD2>("ftdd", n_grid_cells, ELM::numrad);
auto ftid = create<ArrayD2>("ftid", n_grid_cells, ELM::numrad);
auto ftii = create<ArrayD2>("ftii", n_grid_cells, ELM::numrad);
auto flx_absdv = create<ArrayD2>("flx_absdv", n_grid_cells, ELM::nlevsno+1);
auto flx_absdn = create<ArrayD2>("flx_absdn", n_grid_cells, ELM::nlevsno+1);
auto flx_absiv = create<ArrayD2>("flx_absiv", n_grid_cells, ELM::nlevsno+1);
auto flx_absin = create<ArrayD2>("flx_absin", n_grid_cells, ELM::nlevsno+1);

auto mss_cnc_aer_in_fdb = create<ArrayD3>("mss_cnc_aer_in_fdb", n_grid_cells, ELM::nlevsno, ELM::SNICAR::sno_nbr_aer);

auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);
auto albsnd = create<ArrayD2>("albsnd", n_grid_cells, ELM::numrad);
auto albsni = create<ArrayD2>("albsni", n_grid_cells, ELM::numrad);

auto flx_absd_snw = create<ArrayD3>("flx_absd_snw", n_grid_cells, ELM::nlevsno+1, ELM::numrad);
auto flx_absi_snw = create<ArrayD3>("flx_absi_snw", n_grid_cells, ELM::nlevsno+1, ELM::numrad);


auto esai = create<ArrayD1>("esai", n_grid_cells);
auto tlai = create<ArrayD1>("tlai", n_grid_cells);
auto tsai = create<ArrayD1>("tsai", n_grid_cells);

auto nrad = create<ArrayI1>("nrad", n_grid_cells);
auto ncan = create<ArrayI1>("ncan", n_grid_cells);

auto tlai_z = create<ArrayD2>("tlai_z", n_grid_cells, ELM::nlevcan);
auto tsai_z = create<ArrayD2>("tsai_z", n_grid_cells, ELM::nlevcan);
auto fsun_z = create<ArrayD2>("fsun_z", n_grid_cells, ELM::nlevcan);
auto fabd_sun_z = create<ArrayD2>("fabd_sun_z", n_grid_cells, ELM::nlevcan);
auto fabd_sha_z = create<ArrayD2>("fabd_sha_z", n_grid_cells, ELM::nlevcan);
auto fabi_sun_z = create<ArrayD2>("fabi_sun_z", n_grid_cells, ELM::nlevcan);
auto fabi_sha_z = create<ArrayD2>("fabi_sha_z", n_grid_cells, ELM::nlevcan);



auto t_veg = create<ArrayD1>("t_veg", n_grid_cells);
auto fwet = create<ArrayD1>("fwet", n_grid_cells);

auto xl = create<ArrayD1>("xl", ELM::numpft);


auto rhol = create<ArrayD2>("rhol", ELM::numrad, ELM::numpft);
auto rhos = create<ArrayD2>("rhos", ELM::numrad, ELM::numpft);
auto taul = create<ArrayD2>("taul", ELM::numrad, ELM::numpft);
auto taus = create<ArrayD2>("taus", ELM::numrad, ELM::numpft);

auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);



auto albsat = create<ArrayD2>("albsat", n_grid_cells, ELM::numrad);
auto albdry = create<ArrayD2>("albdry", n_grid_cells, ELM::numrad);

auto h2osoi_vol = create<ArrayD2>("h2osoi_vol", n_grid_cells, ELM::nlevgrnd);



 ELM::SurfaceAlbedo::InitTimestep(Land.urbpoi, elai[idx], mss_cnc_bcphi[idx], mss_cnc_bcpho[idx], 
  mss_cnc_dst1[idx], mss_cnc_dst2[idx], mss_cnc_dst3[idx], mss_cnc_dst4[idx],
  vcmaxcintsun[idx], vcmaxcintsha[idx], albsod[idx], albsoi[idx], albgrd[idx], albgri[idx], albd[idx], 
  albi[idx], fabd[idx], fabd_sun[idx], fabd_sha[idx], fabi[idx], fabi_sun[idx], fabi_sha[idx], 
  ftdd[idx], ftid[idx], ftii[idx], flx_absdv[idx], flx_absdn[idx], flx_absiv[idx], flx_absin[idx], 
  mss_cnc_aer_in_fdb[idx]);

 ELM::SurfaceAlbedo::SoilAlbedo(
  Land, snl[idx], t_grnd[idx], coszen[idx], 
  h2osoi_vol[idx], albsat[idx], albdry[idx],
  albsod[idx], albsoi[idx]);

}






