

#include <iostream>
#include <string>


// utilities
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"

// constants
#include "elm_constants.h"

// input data readers and structs
#include "vegdata.h"
#include "atm_data_manager.h"
#include "read_soil_colors.h"
#include "snicar_data.h"
#include "satellite_phenology.h"

#include "InitSoil.hh"




#include "InitSnowLayers.hh"
#include "InitTimestep.hh"
#include "InitTopography.hh"
#include "landtype.h"
//#include "read_atmosphere.h"
//#include "ReadTestData.hh"

#include "incident_shortwave.h"

#include "canopy_hydrology.h"
#include "surface_radiation.h"
#include "canopy_temperature.h"
#include "bareground_fluxes.h"


using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;
using ArrayD3 = ELM::Array<double, 3>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(D0, D1); }
template <class Array_t, class T> Array_t create(const std::string &name, int D0, T D1) { return Array_t(D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(D0, D1, D2); }
template <class Array_t, class T> Array_t create(const std::string &name, int D0, int D1, T D2) { return Array_t(D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }



template<ELM::atm_data_manager::AtmForcType ftype>
using atm_forc_util = ELM::atm_data_manager::AtmDataManager<ArrayD1, ArrayD2, ftype>;

template <ELM::atm_data_manager::AtmForcType ftype>
atm_forc_util<ftype> create_forc_util(const std::string& filename, const ELM::Utils::Date &file_start_time, const int ntimes, const int ncells) {
  return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }


int main(int argc, char **argv) {

  int MPI_COMM_WORLD;
  const int n_procs = 1;
  const int ncells = 1;
  const int ntimes = 100;
  const int myrank = 0;
  const double dtime = 1800.0;
  const double dtime_d = 1800.0 / 86400.0;
  const auto start = ELM::Utils::Date(2012, 7, 1);
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;

  // this is what sat phenology wants
  auto vtype = create<ArrayI1>("vtype", ncells); // pft type
  assign(vtype, 12);


// hardwired grid info
double dz_hardwire[] = {
0.0, 0.0, 0.0,
0.0, 0.0, 0.017512817916255204,
0.02757896925967625, 0.0454700332424132, 0.07496741098620856,
0.12360036510228053, 0.20378255101043175, 0.33598062644843263,
0.5539384053686849, 0.9132900315890611, 1.5057607013992766,
2.482579696981332, 4.0930819526214, 6.7483512780057175,
11.12615029420442, 13.851152141963599 };

double zsoi_hardwire[] = {
0.0, 0.0, 0.0,
0.0, 0.0, 0.007100635417193535,
0.02792500041531687, 0.06225857393654604, 0.11886506690014327,
0.21219339590896316, 0.3660657971047043, 0.6197584979298266,
1.0380270500015696, 1.7276353086671965, 2.8646071131796917,
4.73915671146575, 7.829766507142356, 12.92532061670855,
21.32646906315379, 35.17762120511739 };

double zisoi_hardwire[] = {
0.0, 0.0, 0.0,
0.0, 0.0, 0.0,
0.017512817916255204, 0.04509178717593146, 0.09056182041834465, 
0.16552923140455322, 0.28912959650683373, 0.4929121475172655,
0.8288927739656982, 1.382831179334383, 2.2961212109234443,
3.8018819123227208, 6.284461609304053, 10.377543561925453,
17.12589483993117, 28.252045134135592, 42.10319727609919 };

auto dz = create<ArrayD2>("dz", ncells, ELM::nlevsno + ELM::nlevgrnd, dz_hardwire);
auto zsoi = create<ArrayD2>("zsoi", ncells, ELM::nlevsno + ELM::nlevgrnd, zsoi_hardwire);
auto zisoi = create<ArrayD2>("zisoi", ncells, ELM::nlevsno + ELM::nlevgrnd + 1, zisoi_hardwire);


  auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
  std::cout << " proc_decomp: " << proc_decomp[0] << "," << proc_decomp[1] << std::endl;

  auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { 1, 1 },
          { 0, 0 });

// subsurface params that should be defined by hydrology code:
// ELM uses soil texture to estimate soil hydraulic properties via Clapp & Hornberger
// input vars are percent sand/clay & organic matter

std::string fname_surfdata(
"/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850.nc");

auto isoicol = create<ArrayI1>("isoicol", ncells);
auto albsat = create<ArrayD2>("albsat", ncells, 2);
auto albdry = create<ArrayD2>("albdry", ncells, 2);
ELM::read_soil::read_soil_colors(dd, fname_surfdata, isoicol, albsat, albdry);

auto pct_sand = create<ArrayD2>("pct_sand", ncells, ELM::nlevsoi);
auto pct_clay = create<ArrayD2>("pct_clay", ncells, ELM::nlevsoi);
auto organic = create<ArrayD2>("organic", ncells, ELM::nlevsoi);
ELM::read_soil::read_soil_texture(dd, fname_surfdata, pct_sand, pct_clay, organic);


std::string fname_snicar(
  "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
ELM::snicar_data::SnicarData *snicar_data = new ELM::snicar_data::SnicarData();
ELM::snicar_data::read_snicar_data(dd.comm, fname_snicar, snicar_data);


int nsteps = ntimes + 1; // need n+1 forc data for an n length period of record
std::string fname_forc(
  "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015/all_halfhourly.nc");
auto test_TBOT = create_forc_util<ELM::atm_data_manager::AtmForcType::TBOT>(fname_forc, start, nsteps, ncells);
auto test_PBOT = create_forc_util<ELM::atm_data_manager::AtmForcType::PBOT>(fname_forc, start, nsteps, ncells);
auto test_QBOT = create_forc_util<ELM::atm_data_manager::AtmForcType::RH>(fname_forc, start, nsteps, ncells);
auto test_FLDS = create_forc_util<ELM::atm_data_manager::AtmForcType::FLDS>(fname_forc, start, nsteps, ncells);
auto test_FSDS = create_forc_util<ELM::atm_data_manager::AtmForcType::FSDS>(fname_forc, start, nsteps, ncells);
auto test_PREC = create_forc_util<ELM::atm_data_manager::AtmForcType::PREC>(fname_forc, start, nsteps, ncells);
auto test_WIND = create_forc_util<ELM::atm_data_manager::AtmForcType::WIND>(fname_forc, start, nsteps, ncells);
auto test_ZBOT = create_forc_util<ELM::atm_data_manager::AtmForcType::ZBOT>(fname_forc, start, nsteps, ncells);

std::string data_dir(
  "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata");
std::string fname_pft(
  "lnd/clm2/paramdata/clm_params_c180524.nc");
ELM::VegData<ArrayD1, ArrayD2> vegdata;
  vegdata.read_veg_data(data_dir, fname_pft);


// hardwired
auto smpmin = create<ArrayD1>("smpmin", ncells, -1.0e8);
auto topo_slope = create<ArrayD1>("topo_slope", ncells, 0.070044865858546);
auto topo_std = create<ArrayD1>("topo_std", ncells, 3.96141847422387);
auto t_grnd = create<ArrayD1>("t_grnd", ncells, 275.87936473562456);
auto h2ocan = create<ArrayD1>("t_grnd", ncells, 5.264645812293883e-06); // - BeginGridWaterBalance() - other stuff, too
const double lat = 71.3225;
const double lon = 203.3741;
const double lat_r = lat * ELM::ELM_PI / 180.0;
const double lon_r = lon * ELM::ELM_PI / 180.0;
const double dewmx = 0.1;
const double irrig_rate = 0.0;
const int n_irrig_steps_left = 0;
bool do_capsnow = false;
const bool lakpoi = false;
auto veg_active = create<ArrayB1>("veg_active", ncells); // need value
assign(veg_active, true);                                      // hardwired

auto h2osno_old = create<ArrayD1>("h2osno_old", ncells);   // get rid of this - not needed
auto eflx_bot = create<ArrayD1>("eflx_bot", ncells);       // don't need this, either
auto qflx_glcice = create<ArrayD1>("qflx_glcice", ncells); // don't need this, either

auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", ncells);
auto frac_iceold = create<ArrayD2>("frac_iceold", ncells, ELM::nlevsno + ELM::nlevgrnd);
auto h2osno = create<ArrayD1>("h2osno", ncells);
auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", ncells, ELM::nlevsno + ELM::nlevgrnd);
auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", ncells, ELM::nlevsno + ELM::nlevgrnd);
auto snw_rds = create<ArrayD2>("snw_rds", ncells, ELM::nlevsno);




// forcing data
auto forc_tbot = create<ArrayD1>("forc_tbot", ncells);
auto forc_thbot = create<ArrayD1>("forc_thbot", ncells);
auto forc_pbot = create<ArrayD1>("forc_pbot", ncells);
auto forc_qbot = create<ArrayD1>("forc_qbot", ncells);
auto forc_rh = create<ArrayD1>("forc_rh", ncells);
auto forc_lwrad = create<ArrayD1>("forc_lwrad", ncells);
auto forc_solai = create<ArrayD2>("forc_solai", ncells, 2);
auto forc_solad = create<ArrayD2>("forc_solad", ncells, 2);
auto cos1 = create<ArrayD1>("coszen", ncells);
auto forc_rain = create<ArrayD1>("forc_rain", ncells);
auto forc_snow = create<ArrayD1>("forc_snow", ncells);
auto forc_u = create<ArrayD1>("forc_u", ncells);
auto forc_v = create<ArrayD1>("forc_v", ncells);
auto forc_hgt = create<ArrayD1>("forc_hgt", ncells);
auto forc_hgt_u = create<ArrayD1>("forc_hgt_u", ncells);
auto forc_hgt_t = create<ArrayD1>("forc_hgt_t", ncells);
auto forc_hgt_q = create<ArrayD1>("forc_hgt_q", ncells);

auto forc_vp = create<ArrayD1>("forc_vp", ncells);
auto forc_rho = create<ArrayD1>("forc_rho", ncells);
auto forc_po2 = create<ArrayD1>("forc_po2", ncells);
auto forc_pco2 = create<ArrayD1>("forc_pco2", ncells);



auto snow_depth = create<ArrayD1>("snow_depth", ncells, 0.0); // NEED VALUES! - probably always init at 0
auto frac_sno = create<ArrayD1>("frac_sno", ncells, 0.0);     // NEED VALUES!  \ if not glc, icemec, etc, always init these @ 0.0
auto int_snow = create<ArrayD1>("int_snow", ncells, 0.0);     // NEED VALUES!


 // allocate 2 months of LAI, SAI, HTOP, HBOT for interpolation
auto mlai2t = create<ArrayD2>("mlai2t", ncells, 2);
auto msai2t = create<ArrayD2>("msai2t", ncells, 2);
auto mhvt2t = create<ArrayD2>("mhvt2t", ncells, 2);
auto mhvb2t = create<ArrayD2>("mhvb2t", ncells, 2);

// calculate monthly weights; read 2 months of LAI, SAI, HTOP, HBOT if required
auto timwt = create<ArrayD1>("timwt", 2);            // time weights for lai/sai/hgt
const int n_pfts = 17;
std::string basename_phen("");


// prescribed sat phenology
auto tlai = create<ArrayD1>("tlai", ncells);
auto tsai = create<ArrayD1>("tsai", ncells);
auto elai = create<ArrayD1>("elai", ncells);
auto esai = create<ArrayD1>("esai", ncells);
auto htop = create<ArrayD1>("htop", ncells);
auto hbot = create<ArrayD1>("hbot", ncells);
auto frac_veg_nosno_alb = create<ArrayI1>("frac_veg_nosno_alb", ncells);




ELM::Utils::Date time_plus_half_dt(start);
time_plus_half_dt.increment_seconds(900);

test_TBOT.read_atm_forcing(dd, start, nsteps);
test_PBOT.read_atm_forcing(dd, start, nsteps);
test_QBOT.read_atm_forcing(dd, start, nsteps);
test_FLDS.read_atm_forcing(dd, start, nsteps);
test_FSDS.read_atm_forcing(dd, start, nsteps);
test_PREC.read_atm_forcing(dd, start, nsteps);
test_WIND.read_atm_forcing(dd, start, nsteps);
test_ZBOT.read_atm_forcing(dd, start, nsteps);

auto watsat = create<ArrayD2>("watsat", ncells, ELM::nlevgrnd);
auto sucsat = create<ArrayD2>("sucsat", ncells, ELM::nlevgrnd);
auto bsw = create<ArrayD2>("bsw", ncells, ELM::nlevgrnd);
auto watdry = create<ArrayD2>("watdry", ncells, ELM::nlevgrnd);
auto watopt = create<ArrayD2>("watopt", ncells, ELM::nlevgrnd);
auto watfc = create<ArrayD2>("watfc", ncells, ELM::nlevgrnd);
ELM::init_soil_hydraulics(pct_sand[0], pct_clay[0], organic[0], zsoi[0], 
  watsat[0], bsw[0], sucsat[0], watdry[0], watopt[0], watfc[0]);


auto n_melt = create<ArrayD1>("n_melt", ncells);
auto micro_sigma = create<ArrayD1>("micro_sigma", ncells);
auto snl = create<ArrayI1>("snl", ncells);

ELM::InitTopoSlope(topo_slope[0]);
ELM::InitMicroTopo(Land.ltype, topo_slope[0], topo_std[0], n_melt[0], micro_sigma[0]);
ELM::InitSnowLayers(snow_depth[0], lakpoi, snl[0], dz[0], zsoi[0], zisoi[0]);






// for Canopy hydrology
auto qflx_prec_grnd = create<ArrayD1>("qflx_prec_grnd", ncells);
auto qflx_snwcp_liq = create<ArrayD1>("qflx_snwcp_liq", ncells);
auto qflx_snwcp_ice = create<ArrayD1>("qflx_snwcp_ice", ncells);
auto qflx_snow_grnd = create<ArrayD1>("qflx_snow_grnd", ncells);
auto qflx_rain_grnd = create<ArrayD1>("qflx_rain_grnd", ncells);
auto fwet = create<ArrayD1>("fwet", ncells);
auto fdry = create<ArrayD1>("fdry", ncells);
auto qflx_snow_melt = create<ArrayD1>("qflx_snow_melt", ncells);
auto h2osfc = create<ArrayD1>("h2osfc", ncells);
auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", ncells);
auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", ncells);
auto swe_old = create<ArrayD2>("swe_old", ncells, ELM::nlevsno);
auto t_soisno = create<ArrayD2>("t_soisno", ncells, ELM::nlevsno + ELM::nlevgrnd);



// for surface rad
auto sabg_soil = create<ArrayD1>("sabg_soil", ncells);
auto sabg_snow = create<ArrayD1>("sabg_snow", ncells);
auto sabg = create<ArrayD1>("sabg", ncells);
auto sabv = create<ArrayD1>("sabv", ncells);
auto fsa = create<ArrayD1>("fsa", ncells);
auto fsr = create<ArrayD1>("fsr", ncells);
auto sabg_lyr = create<ArrayD2>("sabg_lyr", ncells, ELM::nlevsno + 1);
auto ftdd = create<ArrayD2>("ftdd", ncells, ELM::numrad);
auto ftid = create<ArrayD2>("ftid", ncells, ELM::numrad);
auto ftii = create<ArrayD2>("ftii", ncells, ELM::numrad);
auto fabd = create<ArrayD2>("fabd", ncells, ELM::numrad);
auto fabi = create<ArrayD2>("fabi", ncells, ELM::numrad);
auto albsod = create<ArrayD2>("albsod", ncells, ELM::numrad);
auto albsoi = create<ArrayD2>("albsoi", ncells, ELM::numrad);
auto albsnd_hst = create<ArrayD2>("albsnd_hst", ncells, ELM::numrad);
auto albsni_hst = create<ArrayD2>("albsni_hst", ncells, ELM::numrad);
auto albgrd = create<ArrayD2>("albgrd", ncells, ELM::numrad);
auto albgri = create<ArrayD2>("albgri", ncells, ELM::numrad);
auto flx_absdv = create<ArrayD2>("flx_absdv", ncells, ELM::nlevsno + 1);
auto flx_absdn = create<ArrayD2>("flx_absdn", ncells, ELM::nlevsno + 1);
auto flx_absiv = create<ArrayD2>("flx_absiv", ncells, ELM::nlevsno + 1);
auto flx_absin = create<ArrayD2>("flx_absin", ncells, ELM::nlevsno + 1);
auto albd = create<ArrayD2>("albd", ncells, ELM::numrad);
auto albi = create<ArrayD2>("albi", ncells, ELM::numrad);





// variables for CanopyTemperature
auto t_h2osfc = create<ArrayD1>("t_h2osfc", ncells);
auto t_h2osfc_bef = create<ArrayD1>("t_h2osfc_bef", ncells);
auto z_0_town = create<ArrayD1>("z_0_town", ncells);
auto z_d_town = create<ArrayD1>("z_d_town", ncells);
auto soilalpha = create<ArrayD1>("soilalpha", ncells);
auto soilalpha_u = create<ArrayD1>("soilalpha_u", ncells);
auto soilbeta = create<ArrayD1>("soilbeta", ncells);
auto qg_snow = create<ArrayD1>("qg_snow", ncells);
auto qg_soil = create<ArrayD1>("qg_soil", ncells);
auto qg = create<ArrayD1>("qg", ncells);
auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", ncells);
auto dqgdT = create<ArrayD1>("dqgdT", ncells);
auto htvp = create<ArrayD1>("htvp", ncells);
auto emg = create<ArrayD1>("emg", ncells);
auto emv = create<ArrayD1>("emv", ncells);
auto z0mg = create<ArrayD1>("z0mg", ncells);
auto z0hg = create<ArrayD1>("z0hg", ncells);
auto z0qg = create<ArrayD1>("z0qg", ncells);
auto z0mv = create<ArrayD1>("z0mv", ncells);
auto z0hv = create<ArrayD1>("z0hv", ncells);
auto z0qv = create<ArrayD1>("z0qv", ncells);
auto thv = create<ArrayD1>("thv", ncells);
auto z0m = create<ArrayD1>("z0m", ncells);
auto displa = create<ArrayD1>("displa", ncells);
auto thm = create<ArrayD1>("thm", ncells);
auto eflx_sh_tot = create<ArrayD1>("eflx_sh_tot", ncells);
auto eflx_sh_tot_u = create<ArrayD1>("eflx_sh_tot_u", ncells);
auto eflx_sh_tot_r = create<ArrayD1>("eflx_sh_tot_r", ncells);
auto eflx_lh_tot = create<ArrayD1>("eflx_lh_tot", ncells);
auto eflx_lh_tot_u = create<ArrayD1>("eflx_lh_tot_u", ncells);
auto eflx_lh_tot_r = create<ArrayD1>("eflx_lh_tot_r", ncells);
auto eflx_sh_veg = create<ArrayD1>("eflx_sh_veg", ncells);
auto qflx_evap_tot = create<ArrayD1>("qflx_evap_tot", ncells);
auto qflx_evap_veg = create<ArrayD1>("qflx_evap_veg", ncells);
auto qflx_tran_veg = create<ArrayD1>("qflx_tran_veg", ncells);
auto tssbef = create<ArrayD2>("tssbef", ncells, ELM::nlevgrnd + ELM::nlevsno);
//auto watsat = create<ArrayD2>("watsat", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
//auto sucsat = create<ArrayD2>("sucsat", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
//auto bsw = create<ArrayD2>("bsw", n_grid_cells, ELM::nlevgrnd);       // comes from SoilStateType.F90
//auto watdry = create<ArrayD2>("watdry", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
//auto watopt = create<ArrayD2>("watopt", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
auto rootfr_road_perv =
    create<ArrayD2>("rootfr_road_perv", ncells, ELM::nlevgrnd); // comes from SoilStateType.F90
auto rootr_road_perv =
    create<ArrayD2>("rootr_road_perv", ncells, ELM::nlevgrnd); // comes from SoilStateType.F90
//auto watfc = create<ArrayD2>("watfc", n_grid_cells, ELM::nlevgrnd);  // comes from SoilStateType.F90





// bareground fluxes
//auto forc_u = create<ArrayD1>("forc_u", ncells);
//auto forc_v = create<ArrayD1>("forc_v", ncells);
auto dlrad = create<ArrayD1>("dlrad", ncells);
auto ulrad = create<ArrayD1>("ulrad", ncells);
//auto forc_rho = create<ArrayD1>("forc_rho", ncells);
auto eflx_sh_grnd = create<ArrayD1>("eflx_sh_grnd", ncells);
auto eflx_sh_snow = create<ArrayD1>("eflx_sh_snow", ncells);
auto eflx_sh_soil = create<ArrayD1>("eflx_sh_soil", ncells);
auto eflx_sh_h2osfc = create<ArrayD1>("eflx_sh_h2osfc", ncells);
auto qflx_evap_soi = create<ArrayD1>("qflx_evap_soi", ncells);
auto qflx_ev_snow = create<ArrayD1>("qflx_ev_snow", ncells);
auto qflx_ev_soil = create<ArrayD1>("qflx_ev_soil", ncells);
auto qflx_ev_h2osfc = create<ArrayD1>("qflx_ev_h2osfc", ncells);
auto t_ref2m = create<ArrayD1>("t_ref2m", ncells);
auto t_ref2m_r = create<ArrayD1>("t_ref2m_r", ncells);
auto q_ref2m = create<ArrayD1>("q_ref2m", ncells);
auto rh_ref2m = create<ArrayD1>("rh_ref2m", ncells);
auto rh_ref2m_r = create<ArrayD1>("rh_ref2m_r", ncells);
auto cgrnds = create<ArrayD1>("cgrnds", ncells);
auto cgrndl = create<ArrayD1>("cgrndl", ncells);
auto cgrnd = create<ArrayD1>("cgrnd", ncells);







// canopy fluxes 
auto altmax_indx = create<ArrayI1>("altmax_indx", ncells);
auto altmax_lastyear_indx = create<ArrayI1>("altmax_lastyear_indx", ncells);

auto max_dayl = create<ArrayD1>("max_dayl", ncells);
auto dayl = create<ArrayD1>("dayl", ncells);
auto tc_stress = create<ArrayD1>("tc_stres", ncells);
//auto forc_lwrad = create<ArrayD1>("forc_lwrad", ncells);
auto t10 = create<ArrayD1>("t10", ncells);
auto vcmaxcintsha = create<ArrayD1>("vcmaxcintsha", ncells);
auto vcmaxcintsun = create<ArrayD1>("vcmaxcintsun", ncells);
//auto forc_pco2 = create<ArrayD1>("forc_pco2", ncells);
//auto forc_po2 = create<ArrayD1>("forc_po2", ncells);
auto btran = create<ArrayD1>("btran", ncells);
auto t_veg = create<ArrayD1>("t_veg", ncells);

auto rootfr = create<ArrayD2>("rootfr", ncells, ELM::nlevgrnd);
auto smpso = create<ArrayD2>("smpso", ncells, ELM::numpft);
auto smpsc = create<ArrayD2>("smpsc", ncells, ELM::numpft);
auto dleaf = create<ArrayD2>("dleaf", ncells, ELM::numpft);
auto rootr = create<ArrayD2>("rootr", ncells, ELM::nlevgrnd);
auto eff_porosity = create<ArrayD2>("eff_porosity", ncells, ELM::nlevgrnd);


ELM::Utils::Date current(start);

for (int t = 0; t < ntimes; ++t) {

ELM::Utils::Date time_plus_half_dt(current);
time_plus_half_dt.increment_seconds(900);


// coszen
double coszen;
{
  ELM::Utils::Date forc_dt_start{test_FSDS.get_data_start_time()};
  forc_dt_start.increment_seconds(round(test_FSDS.forc_t_idx(time_plus_half_dt, test_FSDS.get_data_start_time()) * test_FSDS.get_forc_dt_secs()));
  //std::cout << "forc_dt_start: " << forc_dt_start.doy << "  " << forc_dt_start.sec << "  " << forc_dt_start.year << std::endl;
  //std::cout << "forc_dt_start: " << test_FSDS.forc_t_idx(time_plus_half_dt, test_FSDS.get_data_start_time()) << "  " << test_FSDS.get_forc_dt_secs() << std::endl;
  //std::cout << "time_plus_half_dt: " << time_plus_half_dt.doy << "  " << time_plus_half_dt.sec << std::endl;
  //std::cout << "current: " << current.doy << "  " << current.sec << std::endl;

  int cosz_doy = current.doy + 1;
  double declin = ELM::incident_shortwave::declination_angle2(cosz_doy); // should be the same for model and forcing - don't cross day barrier in forcing timeseries
  double cosz_decday = ELM::Utils::decimal_doy(current) + 1.0;
  double cosz_forc_decday = ELM::Utils::decimal_doy(forc_dt_start) + 1.0;
  double cosz_forcdt_avg = ELM::incident_shortwave::average_cosz(lat_r, lon_r, declin, test_FSDS.get_forc_dt_secs(), cosz_forc_decday);
  double thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, cosz_decday + dtime_d / 2.0);

  double cosz_factor = (thiscosz > 0.001) ? std::min(thiscosz/cosz_forcdt_avg, 10.0) : 0.0;
  coszen = thiscosz * cosz_factor;
}
std::cout << "coszen: " << coszen << std::endl;







ELM::InitTimestep(lakpoi, h2osno[0], veg_active[0], snl[0], h2osoi_ice[0], h2osoi_liq[0],
                    frac_veg_nosno_alb[0], h2osno_old[0], do_capsnow, eflx_bot[0], qflx_glcice[0],
                    frac_veg_nosno[0], frac_iceold[0]);

ELM::interp_monthly_veg<ArrayD3>(current, dtime, n_pfts, fname_surfdata, basename_phen, dd, vtype, timwt, mlai2t, msai2t, mhvt2t,
                      mhvb2t);

test_TBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_thbot);
test_PBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot);
test_QBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_pbot, forc_qbot, forc_rh);
test_FLDS.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot, forc_qbot, forc_tbot, forc_lwrad);
test_FSDS.get_atm_forcing(dtime_d, time_plus_half_dt, cos1, forc_solai, forc_solad);
test_PREC.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_rain, forc_snow);
test_WIND.get_atm_forcing(dtime_d, time_plus_half_dt, forc_u, forc_v);
test_ZBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_hgt, forc_hgt_u, forc_hgt_t,  forc_hgt_q);

ELM::atm_forcing_physics::ConstitutiveAirProperties air_constitutive(forc_qbot, forc_pbot, 
                             forc_tbot, forc_vp,
                             forc_rho, forc_po2, forc_pco2);
air_constitutive(0);
//for (auto i : forc_rho) std::cout << "forc_rho: " << i << std::endl;
//for (auto i : forc_pco2) std::cout << "forc_pco2: " << i << std::endl;
//for (auto i : forc_po2) std::cout << "forc_po2: " << i << std::endl;  

ELM::satellite_phenology(mlai2t[0], msai2t[0], mhvt2t[0], mhvb2t[0], timwt, Land.vtype, snow_depth[0],
                          frac_sno[0], tlai[0], tsai[0], htop[0], hbot[0], elai[0], esai[0],
                          frac_veg_nosno_alb[0]);





int idx = 0; 
{                                  // canopyhydrology calls
    // local vars - these need to be thread local in parallel runs
    static const double dewmx = 0.1; // hardwired
    // static const bool do_capsnow = false; // hardwired
    static const int oldfflag = 1; // hardwired
    double qflx_candrip;
    double qflx_through_snow;
    double qflx_through_rain;
    double fracsnow;
    double fracrain;

    double qflx_irrig[1] = {0.0}; // hardwired here

    ELM::canopy_hydrology::interception(Land, frac_veg_nosno[idx], forc_rain[idx], forc_snow[idx], dewmx, elai[idx], esai[idx], dtime,
                      h2ocan[idx], qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain);

    ELM::canopy_hydrology::ground_flux(Land, do_capsnow, frac_veg_nosno[idx], forc_rain[idx], forc_snow[idx], qflx_irrig[idx],
                    qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain, qflx_prec_grnd[idx],
                    qflx_snwcp_liq[idx], qflx_snwcp_ice[idx], qflx_snow_grnd[idx], qflx_rain_grnd[idx]);

    ELM::canopy_hydrology::fraction_wet(Land, frac_veg_nosno[idx], dewmx, elai[idx], esai[idx], h2ocan[idx], fwet[idx], fdry[idx]);

    ELM::canopy_hydrology::snow_init(Land, dtime, do_capsnow, oldfflag, forc_tbot[idx], t_grnd[idx], qflx_snow_grnd[idx], qflx_snow_melt[idx],
                  n_melt[idx], snow_depth[idx], h2osno[idx], int_snow[idx], swe_old[idx], h2osoi_liq[idx],
                  h2osoi_ice[idx], t_soisno[idx], frac_iceold[idx], snl[idx], dz[idx], zsoi[idx], zisoi[idx], snw_rds[idx],
                  frac_sno_eff[idx], frac_sno[idx]);

    ELM::canopy_hydrology::fraction_h2osfc(Land, micro_sigma[idx], h2osno[idx], h2osfc[idx], h2osoi_liq[idx], frac_sno[idx], frac_sno_eff[idx],
                    frac_h2osfc[idx]);

  }

  //for (auto i : h2ocan) std::cout << "h2ocan: " << i << std::endl;
  //for (auto i : fwet) std::cout << "fwet: " << i << std::endl;
  //for (auto i : fdry) std::cout << "fdry: " << i << std::endl;
  //for (auto i : qflx_prec_grnd) std::cout << "qflx_prec_grnd: " << i << std::endl;
  //for (auto i : qflx_snow_grnd) std::cout << "qflx_snow_grnd: " << i << std::endl;
  //for (auto i : qflx_rain_grnd) std::cout << "qflx_rain_grnd: " << i << std::endl;




{
  // local to these kernel calls
    double trd[ELM::numrad] = {0.0,0.0};
    double tri[ELM::numrad] = {0.0,0.0};

    // call SurfaceRadiation kernels
    ELM::surface_radiation::initialize_flux(Land, sabg_soil[idx], sabg_snow[idx], sabg[idx], sabv[idx], fsa[idx],
                          sabg_lyr[idx]);

    ELM::surface_radiation::total_absorbed_radiation(Land, snl[idx], ftdd[idx], ftid[idx],
         ftii[idx], forc_solad[idx],
         forc_solai[idx], fabd[idx],
         fabi[idx], albsod[idx],
         albsoi[idx], albsnd_hst[idx],
         albsni_hst[idx], albgrd[idx],
         albgri[idx], sabv[idx], fsa[idx], sabg[idx], sabg_soil[idx], sabg_snow[idx],
         trd, tri);


    ELM::surface_radiation::layer_absorbed_radiation(Land, snl[idx], sabg[idx], sabg_snow[idx], snow_depth[idx], flx_absdv[idx],
                       flx_absdn[idx], flx_absiv[idx],
                       flx_absin[idx], trd, tri, sabg_lyr[idx]);

    ELM::surface_radiation::reflected_radiation(Land, albd[idx], albi[idx],
                          forc_solad[idx], forc_solai[idx],
                          fsr[idx]);

}





{
  // call CanopyTemperature kernels
  double qred; // soil surface relative humidity
  double hr;   // relative humidity

    ELM::canopy_temperature::old_ground_temp(Land, t_h2osfc[idx], t_soisno[idx], t_h2osfc_bef[idx], tssbef[idx]);

    ELM::canopy_temperature::ground_temp(Land, snl[idx], frac_sno_eff[idx], frac_h2osfc[idx], t_h2osfc[idx], t_soisno[idx],
                             t_grnd[idx]);

    ELM::canopy_temperature::calc_soilalpha(Land, frac_sno[idx], frac_h2osfc[idx], smpmin[idx], h2osoi_liq[idx], h2osoi_ice[idx],
                            dz[idx], t_soisno[idx], watsat[idx], sucsat[idx], bsw[idx], watdry[idx], watopt[idx],
                            rootfr_road_perv[idx], rootr_road_perv[idx], qred, hr, soilalpha[idx], soilalpha_u[idx]);

    ELM::canopy_temperature::calc_soilbeta(Land, frac_sno[idx], frac_h2osfc[idx], watsat[idx], watfc[idx], h2osoi_liq[idx],
                           h2osoi_ice[idx], dz[idx], soilbeta[idx]);

    ELM::canopy_temperature::humidities(Land, snl[idx], forc_qbot[idx], forc_pbot[idx], t_h2osfc[idx], t_grnd[idx], frac_sno[idx],
                             frac_sno_eff[idx], frac_h2osfc[idx], qred, hr, t_soisno[idx], qg_snow[idx], qg_soil[idx],
                             qg[idx], qg_h2osfc[idx], dqgdT[idx]);

    ELM::canopy_temperature::ground_properties(Land, snl[idx], frac_sno[idx], forc_thbot[idx], forc_qbot[idx], elai[idx], esai[idx], htop[idx],
                           vegdata.displar, vegdata.z0mr,
                           h2osoi_liq[idx], h2osoi_ice[idx],
                           emg[idx], emv[idx], htvp[idx], z0mg[idx], z0hg[idx], z0qg[idx], z0mv[idx], z0hv[idx], z0qv[idx],
                           thv[idx], z0m[idx], displa[idx]);

    ELM::canopy_temperature::forcing_height(Land, veg_active[idx], frac_veg_nosno[idx], forc_hgt_u[idx], forc_hgt_t[idx], forc_hgt_q[idx],
                                 z0m[idx], z0mg[idx], z_0_town[idx], z_d_town[idx], forc_tbot[idx], displa[idx], forc_hgt_u[idx],
                                 forc_hgt_t[idx], forc_hgt_q[idx], thm[idx]);

    ELM::canopy_temperature::init_energy_fluxes(Land, eflx_sh_tot[idx], eflx_sh_tot_u[idx], eflx_sh_tot_r[idx], eflx_lh_tot[idx],
                                eflx_lh_tot_u[idx], eflx_lh_tot_r[idx], eflx_sh_veg[idx], qflx_evap_tot[idx], qflx_evap_veg[idx],
                                qflx_tran_veg[idx]);
}






  

{ // call bareground_fluxes kernels
  int fake_frac_veg_nosno = 0;

  // temporary data to pass between functions
  double zldis;   // reference height "minus" zero displacement height [m]
  double displa;  // displacement height [m]
  double dth;     // diff of virtual temp. between ref. height and surface
  double dqh;     // diff of humidity between ref. height and surface
  double obu;     // Monin-Obukhov length (m)
  double ur;      // wind speed at reference height [m/s]
  double um;      // wind speed including the stablity effect [m/s]
  double temp1;   // relation for potential temperature profile
  double temp2;   // relation for specific humidity profile
  double temp12m; // relation for potential temperature profile applied at 2-m
  double temp22m; // relation for specific humidity profile applied at 2-m
  double ustar;   // friction velocity [m/s]

  ELM::bareground_fluxes::initialize_flux(Land, fake_frac_veg_nosno, forc_u[idx], forc_v[idx], forc_qbot[idx], forc_thbot[idx], forc_hgt_u[idx], thm[idx],
                      thv[idx], t_grnd[idx], qg[idx], z0mg[idx], dlrad[idx], ulrad[idx], zldis, displa, dth, dqh, obu, ur, um);
    
    ELM::bareground_fluxes::stability_iteration(Land, fake_frac_veg_nosno, forc_hgt_t[idx], forc_hgt_u[idx], forc_hgt_q[idx],
                          z0mg[idx], zldis, displa, dth, dqh, ur, forc_qbot[idx], forc_thbot[idx], thv[idx], z0hg[idx],
                          z0qg[idx], obu, um, temp1, temp2, temp12m, temp22m, ustar);
    
    ELM::bareground_fluxes::compute_flux(Land, fake_frac_veg_nosno, snl[idx], forc_rho[idx], soilbeta[idx], dqgdT[idx], htvp[idx], t_h2osfc[idx],
                   qg_snow[idx], qg_soil[idx], qg_h2osfc[idx], t_soisno[idx], forc_pbot[idx], dth,
                   dqh, temp1, temp2, temp12m, temp22m, ustar, forc_qbot[idx], thm[idx], cgrnds[idx], cgrndl[idx], cgrnd[idx],
                   eflx_sh_grnd[idx], eflx_sh_tot[idx], eflx_sh_snow[idx], eflx_sh_soil[idx], eflx_sh_h2osfc[idx],
                   qflx_evap_soi[idx], qflx_evap_tot[idx], qflx_ev_snow[idx], qflx_ev_soil[idx], qflx_ev_h2osfc[idx], t_ref2m[idx],
                   t_ref2m_r[idx], q_ref2m[idx], rh_ref2m[idx], rh_ref2m_r[idx]);
}







current.increment_seconds(1800);

}



}

