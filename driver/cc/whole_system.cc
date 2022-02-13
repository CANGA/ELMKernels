

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
#include "pft_data.h"
#include "atm_data.h"
#include "soil_data.h"
#include "snicar_data.h"
//#include "satellite_phenology.h"

#include "InitSoil.hh"




#include "InitSnowLayers.hh"
#include "InitTimestep.hh"
#include "InitTopography.hh"
#include "land_data.h"
//#include "read_atmosphere.h"
//#include "ReadTestData.hh"

#include "incident_shortwave.h"
#include "day_length.h"

#include "canopy_hydrology.h"
#include "surface_radiation.h"
#include "canopy_temperature.h"
#include "bareground_fluxes.h"
#include "canopy_fluxes.h"
#include "aerosol_data.h"
#include "aerosol_physics.h"
#include "phenology_data.h"
#include "surface_albedo.h"
#include "snow_snicar.h"
#include "RootBioPhys.hh"
#include "surface_fluxes.h"



using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayI2 = ELM::Array<int, 2>;
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


using AtmForcType = ELM::ELMconstants::AtmForcType;

template<AtmForcType ftype>
using atm_forc_util = ELM::AtmDataManager<ArrayD1, ArrayD2, ftype>;

template <AtmForcType ftype>
atm_forc_util<ftype> create_forc_util(const std::string& filename, const ELM::Utils::Date &file_start_time, const int ntimes, const int ncells) {
  return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }


int main(int argc, char **argv) {

int MPI_COMM_WORLD;
const int n_procs = 1;
const int ncells = 1;
const int ntimes = 1008;
const int myrank = 0;
const double dtime = 1800.0;
const double dtime_d = 1800.0 / 86400.0;
const auto start = ELM::Utils::Date(2014, 7, 15);
ELM::LandType Land;
Land.ltype = 1;
Land.ctype = 1;
Land.vtype = 12;

auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
std::cout << " proc_decomp: " << proc_decomp[0] << "," << proc_decomp[1] << std::endl;

auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
        { 1, 1 },
        { 0, 0 });

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


// subsurface params that should be defined by hydrology code:
// ELM uses soil texture to estimate soil hydraulic properties via Clapp & Hornberger
// input vars are percent sand/clay & organic matter

// hardwired
//auto smpmin = create<ArrayD1>("smpmin", ncells, -1.0e8);
auto topo_slope = create<ArrayD1>("topo_slope", ncells, 0.070044865858546);
auto topo_std = create<ArrayD1>("topo_std", ncells, 3.96141847422387);
//auto t_grnd = create<ArrayD1>("t_grnd", ncells, 275.87936473562456);
//auto h2ocan = create<ArrayD1>("h2ocan", ncells, 5.264645812293883e-06); // - BeginGridWaterBalance() - other stuff, too
auto t_grnd = create<ArrayD1>("t_grnd", ncells);
auto h2ocan = create<ArrayD1>("h2ocan", ncells); // - BeginGridWaterBalance() - other stuff, too
const double lat = 71.323;
const double lon = 203.3886;
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
auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", ncells, ELM::nlevsno + ELM::nlevgrnd, 0.0);
auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", ncells, ELM::nlevsno + ELM::nlevgrnd, 0.0);
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

// prescribed sat phenology
auto tlai = create<ArrayD1>("tlai", ncells);
auto tsai = create<ArrayD1>("tsai", ncells);
auto elai = create<ArrayD1>("elai", ncells);
auto esai = create<ArrayD1>("esai", ncells);
auto htop = create<ArrayD1>("htop", ncells);
auto hbot = create<ArrayD1>("hbot", ncells);
auto frac_veg_nosno_alb = create<ArrayI1>("frac_veg_nosno_alb", ncells);

// soil hydraulics
auto watsat = create<ArrayD2>("watsat", ncells, ELM::nlevgrnd);
auto sucsat = create<ArrayD2>("sucsat", ncells, ELM::nlevgrnd);
auto bsw = create<ArrayD2>("bsw", ncells, ELM::nlevgrnd);
auto watdry = create<ArrayD2>("watdry", ncells, ELM::nlevgrnd);
auto watopt = create<ArrayD2>("watopt", ncells, ELM::nlevgrnd);
auto watfc = create<ArrayD2>("watfc", ncells, ELM::nlevgrnd);

// topo, microtopography
auto n_melt = create<ArrayD1>("n_melt", ncells);
auto micro_sigma = create<ArrayD1>("micro_sigma", ncells);
auto snl = create<ArrayI1>("snl", ncells, 0);

// for Canopy hydrology
auto qflx_prec_grnd = create<ArrayD1>("qflx_prec_grnd", ncells);
auto qflx_snwcp_liq = create<ArrayD1>("qflx_snwcp_liq", ncells);
auto qflx_snwcp_ice = create<ArrayD1>("qflx_snwcp_ice", ncells, 0.0);
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
double tsoi[] = {0.0, 0.0, 0.0, 0.0, 0.0, 278.3081064745931, 276.1568781897738,
  275.55803480737063, 275.2677090940866, 274.7286996980052, 273.15, 272.4187794248787, 270.65049816473027,
  267.8224112387398, 265.7450135695632, 264.49481140089864, 264.14163363048056, 264.3351872934207, 264.1163763444719, 263.88852987294865};
auto tsoi_test = create<ArrayD2>("tsoi_test", ncells, ELM::nlevsno + ELM::nlevgrnd, tsoi);
//t_soisno = tsoi_test;

// for can_sun_shade
auto nrad = create<ArrayI1>("nrad", ncells);
auto laisun = create<ArrayD1>("laisun", ncells);
auto laisha = create<ArrayD1>("laisha", ncells);
auto tlai_z = create<ArrayD2>("tlai_z", ncells, ELM::nlevcan);
auto fsun_z = create<ArrayD2>("fsun_z", ncells, ELM::nlevcan);
auto fabd_sun_z = create<ArrayD2>("fabd_sun_z", ncells, ELM::nlevcan);
auto fabd_sha_z = create<ArrayD2>("fabd_sha_z", ncells, ELM::nlevcan);
auto fabi_sun_z = create<ArrayD2>("fabi_sun_z", ncells, ELM::nlevcan);
auto fabi_sha_z = create<ArrayD2>("fabi_sha_z", ncells, ELM::nlevcan);
auto parsun_z = create<ArrayD2>("parsun_z", ncells, ELM::nlevcan);
auto parsha_z = create<ArrayD2>("parsha_z", ncells, ELM::nlevcan);
auto laisun_z = create<ArrayD2>("laisun_z", ncells, ELM::nlevcan);
auto laisha_z = create<ArrayD2>("laisha_z", ncells, ELM::nlevcan);

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
auto t_h2osfc = create<ArrayD1>("t_h2osfc", ncells, 274.0);
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
auto eflx_sh_grnd = create<ArrayD1>("eflx_sh_grnd", ncells, 0.0);
auto eflx_sh_snow = create<ArrayD1>("eflx_sh_snow", ncells, 0.0);
auto eflx_sh_soil = create<ArrayD1>("eflx_sh_soil", ncells, 0.0);
auto eflx_sh_h2osfc = create<ArrayD1>("eflx_sh_h2osfc", ncells, 0.0);
auto qflx_evap_soi = create<ArrayD1>("qflx_evap_soi", ncells, 0.0);
auto qflx_ev_snow = create<ArrayD1>("qflx_ev_snow", ncells, 0.0);
auto qflx_ev_soil = create<ArrayD1>("qflx_ev_soil", ncells, 0.0);
auto qflx_ev_h2osfc = create<ArrayD1>("qflx_ev_h2osfc", ncells, 0.0);
auto t_ref2m = create<ArrayD1>("t_ref2m", ncells);
auto t_ref2m_r = create<ArrayD1>("t_ref2m_r", ncells);
auto q_ref2m = create<ArrayD1>("q_ref2m", ncells);
auto rh_ref2m = create<ArrayD1>("rh_ref2m", ncells);
auto rh_ref2m_r = create<ArrayD1>("rh_ref2m_r", ncells);
auto cgrnds = create<ArrayD1>("cgrnds", ncells);
auto cgrndl = create<ArrayD1>("cgrndl", ncells);
auto cgrnd = create<ArrayD1>("cgrnd", ncells);

// canopy fluxes
//auto max_dayl = create<ArrayD1>("max_dayl", ncells);
//auto dayl = create<ArrayD1>("dayl", ncells);
auto altmax_indx = create<ArrayI1>("altmax_indx", ncells, 5);
auto altmax_lastyear_indx = create<ArrayI1>("altmax_lastyear_indx", ncells, 0);
//auto forc_lwrad = create<ArrayD1>("forc_lwrad", ncells);
auto t10 = create<ArrayD1>("t10", ncells, 276.0);
auto vcmaxcintsha = create<ArrayD1>("vcmaxcintsha", ncells);
auto vcmaxcintsun = create<ArrayD1>("vcmaxcintsun", ncells);
//auto forc_pco2 = create<ArrayD1>("forc_pco2", ncells);
//auto forc_po2 = create<ArrayD1>("forc_po2", ncells);
auto btran = create<ArrayD1>("btran", ncells);
auto t_veg = create<ArrayD1>("t_veg", ncells, 283.0);
auto rootfr = create<ArrayD2>("rootfr", ncells, ELM::nlevgrnd);
auto rootr = create<ArrayD2>("rootr", ncells, ELM::nlevgrnd);
auto eff_porosity = create<ArrayD2>("eff_porosity", ncells, ELM::nlevgrnd);

// surface albedo and snicar
// required for SurfaceAlbedo kernels
// I1
auto snl_top = create<ArrayI1>("snl_top", ncells);
auto snl_btm = create<ArrayI1>("snl_btm", ncells);
//auto snl = create<ArrayI1>("snl", n_grid_cells);
//auto nrad = create<ArrayI1>("nrad", n_grid_cells);
auto ncan = create<ArrayI1>("ncan", ncells);
auto flg_nosnl = create<ArrayI1>("flg_nosnl", ncells);
// I2
auto snw_rds_lcl = create<ArrayI2>("snw_rds_lcl", ncells, ELM::nlevsno);
// D1
//auto coszen = create<ArrayD1>("coszen", n_grid_cells);
//auto elai = create<ArrayD1>("elai", n_grid_cells);
//auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);
//auto esai = create<ArrayD1>("esai", n_grid_cells);
//auto tlai = create<ArrayD1>("tlai", n_grid_cells);
//auto tsai = create<ArrayD1>("tsai", n_grid_cells);
//auto vcmaxcintsun = create<ArrayD1>("vcmaxcintsun", n_grid_cells);
//auto vcmaxcintsha = create<ArrayD1>("vcmaxcintsha", n_grid_cells);
//auto t_veg = create<ArrayD1>("t_veg", n_grid_cells);
//auto fwet = create<ArrayD1>("fwet", n_grid_cells);
//auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
auto mu_not = create<ArrayD1>("mu_not", ncells);
// D2
//auto snw_rds = create<ArrayD2>("snw_rds", n_grid_cells, ELM::nlevsno);
//auto mss_cnc_bcphi = create<ArrayD2>("mss_cnc_bcphi", n_grid_cells, ELM::nlevsno);
//auto mss_cnc_bcpho = create<ArrayD2>("mss_cnc_bcpho", n_grid_cells, ELM::nlevsno);
//auto mss_cnc_dst1 = create<ArrayD2>("mss_cnc_dst1", n_grid_cells, ELM::nlevsno);
//auto mss_cnc_dst2 = create<ArrayD2>("mss_cnc_dst2", n_grid_cells, ELM::nlevsno);
//auto mss_cnc_dst3 = create<ArrayD2>("mss_cnc_dst3", n_grid_cells, ELM::nlevsno);
//auto mss_cnc_dst4 = create<ArrayD2>("mss_cnc_dst4", n_grid_cells, ELM::nlevsno);
//auto albsod = create<ArrayD2>("albsod", n_grid_cells, ELM::numrad);
//auto albsoi = create<ArrayD2>("albsoi", n_grid_cells, ELM::numrad);
//auto albgrd = create<ArrayD2>("albgrd", n_grid_cells, ELM::numrad);
//auto albgri = create<ArrayD2>("albgri", n_grid_cells, ELM::numrad);
//auto albd = create<ArrayD2>("albd", n_grid_cells, ELM::numrad);
//auto albi = create<ArrayD2>("albi", n_grid_cells, ELM::numrad);
//auto fabd = create<ArrayD2>("fabd", n_grid_cells, ELM::numrad);
auto fabd_sun = create<ArrayD2>("fabd_sun", ncells, ELM::numrad);
auto fabd_sha = create<ArrayD2>("fabd_sha", ncells, ELM::numrad);
//auto fabi = create<ArrayD2>("fabi", n_grid_cells, ELM::numrad);
auto fabi_sun = create<ArrayD2>("fabi_sun", ncells, ELM::numrad);
auto fabi_sha = create<ArrayD2>("fabi_sha", ncells, ELM::numrad);
//auto ftdd = create<ArrayD2>("ftdd", n_grid_cells, ELM::numrad);
//auto ftid = create<ArrayD2>("ftid", n_grid_cells, ELM::numrad);
//auto ftii = create<ArrayD2>("ftii", n_grid_cells, ELM::numrad);
//auto flx_absdv = create<ArrayD2>("flx_absdv", n_grid_cells, ELM::nlevsno+1);
//auto flx_absdn = create<ArrayD2>("flx_absdn", n_grid_cells, ELM::nlevsno+1);
//auto flx_absiv = create<ArrayD2>("flx_absiv", n_grid_cells, ELM::nlevsno+1);
//auto flx_absin = create<ArrayD2>("flx_absin", n_grid_cells, ELM::nlevsno+1);
auto albsnd = create<ArrayD2>("albsnd", ncells, ELM::numrad);
auto albsni = create<ArrayD2>("albsni", ncells, ELM::numrad);
//auto tlai_z = create<ArrayD2>("tlai_z", n_grid_cells, ELM::nlevcan);
auto tsai_z = create<ArrayD2>("tsai_z", ncells, ELM::nlevcan);
//auto fsun_z = create<ArrayD2>("fsun_z", n_grid_cells, ELM::nlevcan);
//auto fabd_sun_z = create<ArrayD2>("fabd_sun_z", n_grid_cells, ELM::nlevcan);
//auto fabd_sha_z = create<ArrayD2>("fabd_sha_z", n_grid_cells, ELM::nlevcan);
//auto fabi_sun_z = create<ArrayD2>("fabi_sun_z", n_grid_cells, ELM::nlevcan);
//auto fabi_sha_z = create<ArrayD2>("fabi_sha_z", n_grid_cells, ELM::nlevcan);
//auto albsat = create<ArrayD2>("albsat", n_grid_cells, ELM::numrad);
//auto albdry = create<ArrayD2>("albdry", n_grid_cells, ELM::numrad);
auto h2osoi_vol = create<ArrayD2>("h2osoi_vol", ncells, ELM::nlevgrnd);
// D3
auto mss_cnc_aer_in_fdb = create<ArrayD3>("mss_cnc_aer_in_fdb", ncells, ELM::nlevsno, ELM::snow_snicar::sno_nbr_aer);
auto flx_absd_snw = create<ArrayD3>("flx_absd_snw", ncells, ELM::nlevsno+1, ELM::numrad);
auto flx_absi_snw = create<ArrayD3>("flx_absi_snw", ncells, ELM::nlevsno+1, ELM::numrad);
auto flx_abs_lcl = create<ArrayD3>("flx_abs_lcl", ncells, ELM::nlevsno+1, ELM::snow_snicar::numrad_snw);
//auto flx_abs = create<ArrayD3>("flx_abs", n_grid_cells, ELM::nlevsno+1, ELM::numrad);
// remaining variables required for SNICAR kernels
// D1
//auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);
//auto ss_alb_oc1 = create<ArrayD1>("ss_alb_oc1", ELM::snow_snicar::numrad_snw);
//auto asm_prm_oc1 = create<ArrayD1>("asm_prm_oc1", ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_oc1 = create<ArrayD1>("ext_cff_mss_oc1", ELM::snow_snicar::numrad_snw);
//auto ss_alb_oc2 = create<ArrayD1>("ss_alb_oc2", ELM::snow_snicar::numrad_snw);
//auto asm_prm_oc2 = create<ArrayD1>("asm_prm_oc2", ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_oc2 = create<ArrayD1>("ext_cff_mss_oc2", ELM::snow_snicar::numrad_snw);
//auto ss_alb_dst1 = create<ArrayD1>("ss_alb_dst1", ELM::snow_snicar::numrad_snw);
//auto asm_prm_dst1 = create<ArrayD1>("asm_prm_dst1", ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_dst1 = create<ArrayD1>("ext_cff_mss_dst1", ELM::snow_snicar::numrad_snw);
//auto ss_alb_dst2 = create<ArrayD1>("ss_alb_dst2", ELM::snow_snicar::numrad_snw);
//auto asm_prm_dst2 = create<ArrayD1>("asm_prm_dst2", ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_dst2 = create<ArrayD1>("ext_cff_mss_dst2", ELM::snow_snicar::numrad_snw);
//auto ss_alb_dst3 = create<ArrayD1>("ss_alb_dst3", ELM::snow_snicar::numrad_snw);
//auto asm_prm_dst3 = create<ArrayD1>("asm_prm_dst3", ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_dst3 = create<ArrayD1>("ext_cff_mss_dst3", ELM::snow_snicar::numrad_snw);
//auto ss_alb_dst4 = create<ArrayD1>("ss_alb_dst4", ELM::snow_snicar::numrad_snw);
//auto asm_prm_dst4 = create<ArrayD1>("asm_prm_dst4", ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_dst4 = create<ArrayD1>("ext_cff_mss_dst4", ELM::snow_snicar::numrad_snw);
// D2
//auto ss_alb_snw_drc = create<ArrayD2>("ss_alb_snw_drc", ELM::snow_snicar::numrad_snw, ELM::snow_snicar::idx_Mie_snw_mx);
//auto asm_prm_snw_drc = create<ArrayD2>("asm_prm_snw_drc", ELM::snow_snicar::numrad_snw, ELM::snow_snicar::idx_Mie_snw_mx);
//auto ext_cff_mss_snw_drc = create<ArrayD2>("ext_cff_mss_snw_drc", ELM::snow_snicar::numrad_snw, ELM::snow_snicar::idx_Mie_snw_mx);
//auto ss_alb_snw_dfs = create<ArrayD2>("ss_alb_snw_dfs", ELM::snow_snicar::numrad_snw, ELM::snow_snicar::idx_Mie_snw_mx);
//auto asm_prm_snw_dfs = create<ArrayD2>("asm_prm_snw_dfs", ELM::snow_snicar::numrad_snw, ELM::snow_snicar::idx_Mie_snw_mx);
//auto ext_cff_mss_snw_dfs = create<ArrayD2>("ext_cff_mss_snw_dfs", ELM::snow_snicar::numrad_snw, ELM::snow_snicar::idx_Mie_snw_mx);
//auto ss_alb_bc1 = create<ArrayD2>("ss_alb_bc1", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
//auto asm_prm_bc1 = create<ArrayD2>("asm_prm_bc1", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_bc1 = create<ArrayD2>("ext_cff_mss_bc1", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
//auto ss_alb_bc2 = create<ArrayD2>("ss_alb_bc2", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
//auto asm_prm_bc2 = create<ArrayD2>("asm_prm_bc2", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
//auto ext_cff_mss_bc2 = create<ArrayD2>("ext_cff_mss_bc2", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
auto albout_lcl = create<ArrayD2>("albout_lcl", ncells, ELM::snow_snicar::numrad_snw);
auto flx_slrd_lcl = create<ArrayD2>("flx_slrd_lcl", ncells, ELM::snow_snicar::numrad_snw);
auto flx_slri_lcl = create<ArrayD2>("flx_slri_lcl", ncells, ELM::snow_snicar::numrad_snw);
//auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
//auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto h2osoi_ice_lcl = create<ArrayD2>("h2osoi_ice_lcl", ncells, ELM::nlevsno);
auto h2osoi_liq_lcl = create<ArrayD2>("h2osoi_liq_lcl", ncells, ELM::nlevsno);
// D3
//auto bcenh = create<ArrayD3>("bcenh", ELM::snow_snicar::idx_bcint_icerds_max+1, ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
auto g_star = create<ArrayD3>("g_star", ncells, ELM::snow_snicar::numrad_snw, ELM::nlevsno);
auto omega_star = create<ArrayD3>("omega_star", ncells, ELM::snow_snicar::numrad_snw, ELM::nlevsno);
auto tau_star = create<ArrayD3>("tau_star", ncells, ELM::snow_snicar::numrad_snw, ELM::nlevsno);



// soil fluxes (outputs)
auto eflx_soil_grnd = create<ArrayD1>("eflx_soil_grnd", ncells);
auto qflx_evap_grnd = create<ArrayD1>("qflx_evap_grnd", ncells);
auto qflx_sub_snow = create<ArrayD1>("qflx_sub_snow", ncells);
auto qflx_dew_snow = create<ArrayD1>("qflx_dew_snow", ncells);
auto qflx_dew_grnd = create<ArrayD1>("qflx_dew_grnd", ncells);


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/*                          input files                                                                */
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
   // std::string fname_surfdata(
   // "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850.nc");
   // std::string fname_snicar(
   //   "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
   // std::string fname_forc(
   //   "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc");
   // std::string data_dir(
   //   "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata");
   // std::string fname_pft(
   //   "lnd/clm2/paramdata/clm_params_c180524.nc");
   // std::string fname_aerosol(
   //   "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc");


// for testing
std::string fname_surfdata(
"/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc");
std::string fname_snicar(
  "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
std::string fname_forc(
  "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc");
std::string data_dir(
  "/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata");
std::string fname_pft(
  "lnd/clm2/paramdata/clm_params_c180524.nc");
std::string fname_aerosol(
  "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc");


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// init data containers and read time-invariant data from files
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    //int nsteps = ntimes + 1; // need n+1 forc data for an n length period of record
    int nsteps = 543120; // for testing
    const auto fstart = ELM::Utils::Date(1985, 1, 1);
    auto test_TBOT = create_forc_util<AtmForcType::TBOT>(fname_forc, fstart, nsteps, ncells);
    auto test_PBOT = create_forc_util<AtmForcType::PBOT>(fname_forc, fstart, nsteps, ncells);
    auto test_QBOT = create_forc_util<AtmForcType::RH>(fname_forc, fstart, nsteps, ncells);
    auto test_FLDS = create_forc_util<AtmForcType::FLDS>(fname_forc, fstart, nsteps, ncells);
    auto test_FSDS = create_forc_util<AtmForcType::FSDS>(fname_forc, fstart, nsteps, ncells);
    auto test_PREC = create_forc_util<AtmForcType::PREC>(fname_forc, fstart, nsteps, ncells);
    auto test_WIND = create_forc_util<AtmForcType::WIND>(fname_forc, fstart, nsteps, ncells);
    auto test_ZBOT = create_forc_util<AtmForcType::ZBOT>(fname_forc, fstart, nsteps, ncells);
    
    ELM::AerosolMasses aerosol_masses(ncells);
    ELM::AerosolConcentrations aerosol_concentrations(ncells);
    ELM::PhenologyDataManager phen_data(ncells, 17);
    
    auto isoicol = create<ArrayI1>("isoicol", ncells);
    auto albsat = create<ArrayD2>("albsat", ncells, 2);
    auto albdry = create<ArrayD2>("albdry", ncells, 2);
    ELM::read_soil::read_soil_colors(dd, fname_surfdata, isoicol, albsat, albdry);
    
    auto pct_sand = create<ArrayD2>("pct_sand", ncells, ELM::nlevgrnd);
    auto pct_clay = create<ArrayD2>("pct_clay", ncells, ELM::nlevgrnd);
    auto organic = create<ArrayD2>("organic", ncells, ELM::nlevgrnd);
    ELM::read_soil::read_soil_texture(dd, fname_surfdata, pct_sand, pct_clay, organic);
    
    ELM::SnicarData *snicar_data = new ELM::SnicarData();
    ELM::read_snicar_data(dd.comm, fname_snicar, snicar_data);
    
    ELM::VegData<ArrayD1, ArrayD2> vegdata;
    vegdata.read_veg_data(data_dir, fname_pft);

    //ELM::AerosolDataManager<ArrayD1> aerosoldata;
    ELM::AerosolDataManager aerosol_data;
    aerosol_data.read_data(dd.comm, fname_aerosol, lon, lat);



/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// some free functions that calculate initial values
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    
    ELM::InitTopoSlope(topo_slope[0]);
    ELM::InitMicroTopo(Land.ltype, topo_slope[0], topo_std[0], n_melt[0], micro_sigma[0]);
    ELM::InitSnowLayers(snow_depth[0], lakpoi, snl[0], dz[0], zsoi[0], zisoi[0]);

    ELM::init_soil_hydraulics(pct_sand[0], pct_clay[0], organic[0], zsoi[0], 
      watsat[0], bsw[0], sucsat[0], watdry[0], watopt[0], watfc[0]);
    ELM::init_vegrootfr(vtype[0], vegdata.roota_par[vtype[0]], vegdata.rootb_par[vtype[0]], zisoi[0], rootfr[0]);
    ELM::init_soil_temp(Land, snl[0], t_soisno[0], t_grnd[0]);
    
    // !!!!! testing !!!!!!!!!!
    for (int q = 0; q < 20; ++q) {
      t_soisno(0, q) = tsoi_test(0, q); // for testing
    }
    t_grnd[0] = t_soisno(0, ELM::nlevsno - snl(0));


    ELM::init_snow_state(Land.urbpoi, snl[0], h2osno[0], int_snow[0], snow_depth[0], h2osfc[0], h2ocan[0], frac_h2osfc[0], fwet[0], fdry[0], frac_sno[0], snw_rds[0]);
    ELM::init_soilh2o_state(Land, snl[0], watsat[0], t_soisno[0], dz[0], h2osoi_vol[0], h2osoi_liq[0], h2osoi_ice[0]);

double h2osoi_ice_test[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 51.095355179469955, 131.99213225849098, 17.829256395227745, 95.72899575304584, 155.31526899797177, 0.01, 0.01, 0.01, 0.01, 0.01};
double h2osoi_liq_test[] = {0.0, 0.0, 0.0, 0.0, 0.0, 7.045411435071487, 14.353496179256807, 36.308518784697064, 62.46145027256513, 97.14000248023912, 97.47148319510016, 78.52160092062527, 65.63904088905001, 41.25305599181871, 70.8566046019581, 0.01, 0.01, 0.01, 0.01, 0.01};
double h2osoi_vol_test[] = {0.4016484663460637, 0.5196481455614503, 0.7967166638201649, 0.8331813710901114, 0.7859200286330449, 0.7517405589446893, 0.6621235242027332, 0.1535948180493002, 0.15947477948341815, 0.15954052527228618, 8.420726808634413e-06, 5.107428986500891e-06, 3.0978122726178113e-06, 1.8789181213767733e-06, 1.5092697845407248e-06};
for (int q = 0; q < ELM::nlevgrnd+ELM::nlevsno; ++q) {
  h2osoi_ice(0,q) = h2osoi_ice_test[q];
  h2osoi_liq(0,q) = h2osoi_liq_test[q];
}

for (int q = 0; q < ELM::nlevgrnd; ++q) {
  h2osoi_vol(0,q) = h2osoi_vol_test[q];
}


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/*                          TIME LOOP                                                                  */
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// these calls should be inside time loop
test_TBOT.read_atm_forcing(dd, fstart, nsteps);
test_PBOT.read_atm_forcing(dd, fstart, nsteps);
test_QBOT.read_atm_forcing(dd, fstart, nsteps);
test_FLDS.read_atm_forcing(dd, fstart, nsteps);
test_FSDS.read_atm_forcing(dd, fstart, nsteps);
test_PREC.read_atm_forcing(dd, fstart, nsteps);
test_WIND.read_atm_forcing(dd, fstart, nsteps);
test_ZBOT.read_atm_forcing(dd, fstart, nsteps);

auto coszen = create<ArrayD1>("coszen", ncells);
auto cosz_factor = create<ArrayD1>("cosz_factor", ncells);

ELM::Utils::Date current(start);
ELM::Utils::Date big(start);
int idx = 0; // hardwire for ncells = 1
for (int t = 0; t < ntimes; ++t) {

ELM::Utils::Date time_plus_half_dt(current);
time_plus_half_dt.increment_seconds(900);


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// get coszen - should make this a function
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    double max_dayl;
    double dayl;
    {
      ELM::Utils::Date forc_dt_start{test_FSDS.get_data_start_time()};
      forc_dt_start.increment_seconds(round(test_FSDS.forc_t_idx(time_plus_half_dt, test_FSDS.get_data_start_time()) * test_FSDS.get_forc_dt_secs()));
      //std::cout << "forc_dt_start: " << forc_dt_start.doy << "  " << forc_dt_start.sec << "  " << forc_dt_start.year << std::endl;
      //std::cout << "forc_dt_start: " << test_FSDS.forc_t_idx(time_plus_half_dt, test_FSDS.get_data_start_time()) << "  " << test_FSDS.get_forc_dt_secs() << std::endl;
      //std::cout << "time_plus_half_dt: " << time_plus_half_dt.doy << "  " << time_plus_half_dt.sec << std::endl;
      //std::cout << "current: " << current.doy << "  " << current.sec << std::endl;
    
      int cosz_doy = current.doy + 1;
      //std::cout << "cosz_doy::  " << cosz_doy << std::endl;
      double declin = ELM::incident_shortwave::declination_angle2(cosz_doy); // should be the same for model and forcing - don't cross day barrier in forcing timeseries
      double cosz_decday = ELM::Utils::decimal_doy(current) + 1.0;
      //std::cout << std::setprecision(15) << "cosz_decday::  " << cosz_decday << std::endl;
      double cosz_forc_decday = ELM::Utils::decimal_doy(forc_dt_start) + 1.0;
      //std::cout << "cosz_forc_decday::  " << cosz_forc_decday << std::endl;
      double cosz_forcdt_avg = ELM::incident_shortwave::average_cosz(lat_r, lon_r, declin, test_FSDS.get_forc_dt_secs(), cosz_forc_decday);
      //std::cout << "cosz_forcdt_avg::  " << cosz_forcdt_avg << std::endl;
      double thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, cosz_decday + dtime_d/2.0);
      //std::cout << "thiscosz::  " << thiscosz << std::endl;

      //std::cout << "test_FSDS.get_forc_dt_secs()::  " << test_FSDS.get_forc_dt_secs() << std::endl;
    
      //double cosz_factor = (thiscosz > 0.001) ? std::min(thiscosz/cosz_forcdt_avg, 10.0) : 0.0;
      cosz_factor[0] = (thiscosz > 0.001) ? 1.0 : 0.0;
      //std::cout << "cosz_factor::  " << cosz_factor << std::endl;
      //coszen[0] = thiscosz * cosz_factor;
      coszen[0] = thiscosz;
      //std::cout << "coszen::  " << coszen[0] << std::endl;
      max_dayl = ELM::max_daylength(lat_r);
      dayl = ELM::daylength(lat_r, declin);
    }

/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// read time-variable data
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    phen_data.read_data(dd, fname_surfdata, current, vtype); // if needed
    phen_data.get_data(current, snow_depth, frac_sno, vtype, elai, esai,
    htop, hbot, tlai, tsai, frac_veg_nosno_alb);
    // should read forcing data here

/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// timestep init functions
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    test_TBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_thbot);
    test_PBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot);
    test_QBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_pbot, forc_qbot, forc_rh);
    test_FLDS.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot, forc_qbot, forc_tbot, forc_lwrad);
    test_FSDS.get_atm_forcing(dtime_d, time_plus_half_dt, cosz_factor, forc_solai, forc_solad);
    test_PREC.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_rain, forc_snow);
    test_WIND.get_atm_forcing(dtime_d, time_plus_half_dt, forc_u, forc_v);
    test_ZBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_hgt, forc_hgt_u, forc_hgt_t,  forc_hgt_q);

    // calc constitutive air props - call properly later
    ELM::atm_forcing_physics::ConstitutiveAirProperties air_constitutive(forc_qbot, forc_pbot, 
                             forc_tbot, forc_vp,
                             forc_rho, forc_po2, forc_pco2);
    air_constitutive(0); // - call properly

    // get aerosol mss and cnc
    aerosol_data.invoke_aerosol_source(time_plus_half_dt, dtime, snl, aerosol_masses);
    ELM::aerosols::invoke_aerosol_concen_and_mass(do_capsnow, dtime, snl, h2osoi_liq,
    h2osoi_ice, snw_rds, qflx_snwcp_ice, aerosol_masses, aerosol_concentrations);

    ELM::InitTimestep(lakpoi, h2osno[0], veg_active[0], snl[0], h2osoi_ice[0], h2osoi_liq[0],
                    frac_veg_nosno_alb[0], h2osno_old[0], do_capsnow, eflx_bot[0], qflx_glcice[0],
                    frac_veg_nosno[0], frac_iceold[0]);


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call surface albedo kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{

//for (auto i : snl) std::cout << "snl: " << i << std::endl;
//for (auto i : elai) std::cout << "elai: " << i << std::endl;
//for (auto i : esai) std::cout << "esai: " << i << std::endl;
//for (auto i : tlai) std::cout << "tlai: " << i << std::endl;
//for (auto i : tsai) std::cout << "tsai: " << i << std::endl;
//for (auto i : coszen) std::cout << "coszen: " << i << std::endl;
//for (auto i : frac_sno) std::cout << "frac_sno: " << i << std::endl;
//for (auto i : h2osno) std::cout << "h2osno: " << i << std::endl;
//for (auto i : t_veg) std::cout << "t_veg: " << i << std::endl;
//for (auto i : fwet) std::cout << "fwet: " << i << std::endl;
////for (auto i : mss_cnc_bcphi) std::cout << "mss_cnc_bcphi: " << i << std::endl;
////for (auto i : mss_cnc_bcpho) std::cout << "mss_cnc_bcpho: " << i << std::endl;
////for (auto i : mss_cnc_dst1) std::cout << "mss_cnc_dst1: " << i << std::endl;
////for (auto i : mss_cnc_dst2) std::cout << "mss_cnc_dst2: " << i << std::endl;
////for (auto i : mss_cnc_dst3) std::cout << "mss_cnc_dst3: " << i << std::endl;
////for (auto i : mss_cnc_dst4) std::cout << "mss_cnc_dst4: " << i << std::endl;
//for (auto i : albsod) std::cout << "albsod: " << i << std::endl;
//for (auto i : albsoi) std::cout << "albsoi: " << i << std::endl;
//for (auto i : albsnd) std::cout << "albsnd: " << i << std::endl;
//for (auto i : albsni) std::cout << "albsni: " << i << std::endl;
//for (auto i : h2osoi_liq) std::cout << "h2osoi_liq: " << i << std::endl;
//for (auto i : h2osoi_ice) std::cout << "h2osoi_ice: " << i << std::endl;
//for (auto i : snw_rds) std::cout << "snw_rds: " << i << std::endl;
//for (auto i : h2osoi_vol) std::cout << "h2osoi_vol: " << i << std::endl;
//for (auto i : albsat) std::cout << "albsat: " << i << std::endl;
//for (auto i : albdry) std::cout << "albdry: " << i << std::endl;
//for (auto i : albgrd) std::cout << "albgrd: " << i << std::endl;
//for (auto i : albgri) std::cout << "albgri: " << i << std::endl;
//for (auto i : flx_absdv) std::cout << "flx_absdv: " << i << std::endl;
//for (auto i : flx_absdn) std::cout << "flx_absdn: " << i << std::endl;
//for (auto i : flx_absiv) std::cout << "flx_absiv: " << i << std::endl;
//for (auto i : flx_absin) std::cout << "flx_absin: " << i << std::endl;
//for (auto i : tlai_z) std::cout << "tlai_z: " << i << std::endl;
//for (auto i : tsai_z) std::cout << "tsai_z: " << i << std::endl;
//for (auto i : fsun_z) std::cout << "fsun_z: " << i << std::endl;
//for (auto i : fabd_sun_z) std::cout << "fabd_sun_z: " << i << std::endl;
//for (auto i : fabd_sha_z) std::cout << "fabd_sha_z: " << i << std::endl;
//for (auto i : fabi_sun_z) std::cout << "fabi_sun_z: " << i << std::endl;
//for (auto i : fabi_sha_z) std::cout << "fabi_sha_z: " << i << std::endl;
//for (auto i : albd) std::cout << "albd: " << i << std::endl;
//for (auto i : ftid) std::cout << "ftid: " << i << std::endl;
//for (auto i : ftdd) std::cout << "ftdd: " << i << std::endl;
//for (auto i : fabd) std::cout << "fabd: " << i << std::endl;
//for (auto i : fabd_sun) std::cout << "fabd_sun: " << i << std::endl;
//for (auto i : fabd_sha) std::cout << "fabd_sha: " << i << std::endl;
//for (auto i : albi) std::cout << "albi: " << i << std::endl;
//for (auto i : ftii) std::cout << "ftii: " << i << std::endl;
//for (auto i : fabi) std::cout << "fabi: " << i << std::endl;
//for (auto i : fabi_sun) std::cout << "fabi_sun: " << i << std::endl;
//for (auto i : fabi_sha) std::cout << "fabi_sha: " << i << std::endl;
//for (auto i : t_grnd) std::cout << "t_grnd: " << i << std::endl;


  // parse pft data for Land.vtype
  ELM::AlbedoVegData albveg = vegdata.get_pft_albveg(Land.vtype);

  ELM::surface_albedo::init_timestep(Land.urbpoi, elai[idx], aerosol_concentrations.mss_cnc_bcphi[idx], aerosol_concentrations.mss_cnc_bcpho[idx], 
    aerosol_concentrations.mss_cnc_dst1[idx], aerosol_concentrations.mss_cnc_dst2[idx], aerosol_concentrations.mss_cnc_dst3[idx], aerosol_concentrations.mss_cnc_dst4[idx], vcmaxcintsun[idx], 
    vcmaxcintsha[idx], albsod[idx], albsoi[idx], albgrd[idx], albgri[idx], albd[idx], albi[idx], 
    fabd[idx], fabd_sun[idx], fabd_sha[idx], fabi[idx], fabi_sun[idx], fabi_sha[idx], ftdd[idx], 
    ftid[idx], ftii[idx], flx_absdv[idx], flx_absdn[idx], flx_absiv[idx], flx_absin[idx], 
    mss_cnc_aer_in_fdb[idx]);

  ELM::surface_albedo::soil_albedo(
    Land, snl[idx], t_grnd[idx], coszen[idx], h2osoi_vol[idx], 
    albsat[idx], albdry[idx], albsod[idx], albsoi[idx]);

  {
    int flg_slr_in = 1; // direct-beam

    ELM::snow_snicar::init_timestep (Land.urbpoi, flg_slr_in, coszen[idx], h2osno[idx], snl[idx], 
      h2osoi_liq[idx], h2osoi_ice[idx], snw_rds[idx], snl_top[idx], snl_btm[idx], 
      flx_abs_lcl[idx], flx_absd_snw[idx], flg_nosnl[idx], h2osoi_ice_lcl[idx], 
      h2osoi_liq_lcl[idx], snw_rds_lcl[idx], mu_not[idx], flx_slrd_lcl[idx], 
      flx_slri_lcl[idx]);


    ELM::snow_snicar::snow_aerosol_mie_params(Land.urbpoi, flg_slr_in, snl_top[idx], snl_btm[idx], coszen[idx], 
      h2osno[idx], snw_rds_lcl[idx], h2osoi_ice_lcl[idx], h2osoi_liq_lcl[idx], snicar_data->ss_alb_oc1, snicar_data->asm_prm_oc1, 
      snicar_data->ext_cff_mss_oc1, snicar_data->ss_alb_oc2, snicar_data->asm_prm_oc2, snicar_data->ext_cff_mss_oc2, snicar_data->ss_alb_dst1, snicar_data->asm_prm_dst1, 
      snicar_data->ext_cff_mss_dst1, snicar_data->ss_alb_dst2, snicar_data->asm_prm_dst2, snicar_data->ext_cff_mss_dst2, snicar_data->ss_alb_dst3, snicar_data->asm_prm_dst3, 
      snicar_data->ext_cff_mss_dst3, snicar_data->ss_alb_dst4, snicar_data->asm_prm_dst4, snicar_data->ext_cff_mss_dst4, snicar_data->ss_alb_snw_drc, 
      snicar_data->asm_prm_snw_drc, snicar_data->ext_cff_mss_snw_drc, snicar_data->ss_alb_snw_dfs, snicar_data->asm_prm_snw_dfs, snicar_data->ext_cff_mss_snw_dfs, 
      snicar_data->ss_alb_bc1, snicar_data->asm_prm_bc1, snicar_data->ext_cff_mss_bc1, snicar_data->ss_alb_bc2, snicar_data->asm_prm_bc2, snicar_data->ext_cff_mss_bc2, 
      snicar_data->bcenh, mss_cnc_aer_in_fdb[idx], g_star[idx], omega_star[idx], tau_star[idx]);


    ELM::snow_snicar::snow_radiative_transfer_solver(Land.urbpoi, flg_slr_in,flg_nosnl[idx], snl_top[idx], 
      snl_btm[idx], coszen[idx], h2osno[idx], mu_not[idx], flx_slrd_lcl[idx], flx_slri_lcl[idx], 
      albsoi[idx], g_star[idx], omega_star[idx], tau_star[idx], albout_lcl[idx], flx_abs_lcl[idx]);


    ELM::snow_snicar::snow_albedo_radiation_factor(Land.urbpoi, flg_slr_in, snl_top[idx], coszen[idx], mu_not[idx], 
      h2osno[idx], snw_rds_lcl[idx], albsoi[idx],albout_lcl[idx], flx_abs_lcl[idx], albsnd[idx], 
      flx_absd_snw[idx]);
  }


  {
    int flg_slr_in = 2; // diffuse

    ELM::snow_snicar::init_timestep (Land.urbpoi,
      flg_slr_in, coszen[idx], h2osno[idx], snl[idx], h2osoi_liq[idx], h2osoi_ice[idx],
      snw_rds[idx], snl_top[idx], snl_btm[idx], flx_abs_lcl[idx], flx_absi_snw[idx],
      flg_nosnl[idx], h2osoi_ice_lcl[idx], h2osoi_liq_lcl[idx], snw_rds_lcl[idx], mu_not[idx], 
      flx_slrd_lcl[idx], flx_slri_lcl[idx]);


    ELM::snow_snicar::snow_aerosol_mie_params(Land.urbpoi, flg_slr_in, snl_top[idx], snl_btm[idx], coszen[idx], 
      h2osno[idx], snw_rds_lcl[idx], h2osoi_ice_lcl[idx], h2osoi_liq_lcl[idx], snicar_data->ss_alb_oc1, snicar_data->asm_prm_oc1, 
      snicar_data->ext_cff_mss_oc1, snicar_data->ss_alb_oc2, snicar_data->asm_prm_oc2, snicar_data->ext_cff_mss_oc2, snicar_data->ss_alb_dst1, snicar_data->asm_prm_dst1, 
      snicar_data->ext_cff_mss_dst1, snicar_data->ss_alb_dst2, snicar_data->asm_prm_dst2, snicar_data->ext_cff_mss_dst2, snicar_data->ss_alb_dst3, snicar_data->asm_prm_dst3, 
      snicar_data->ext_cff_mss_dst3, snicar_data->ss_alb_dst4, snicar_data->asm_prm_dst4, snicar_data->ext_cff_mss_dst4, snicar_data->ss_alb_snw_drc, 
      snicar_data->asm_prm_snw_drc, snicar_data->ext_cff_mss_snw_drc, snicar_data->ss_alb_snw_dfs, snicar_data->asm_prm_snw_dfs, snicar_data->ext_cff_mss_snw_dfs, 
      snicar_data->ss_alb_bc1, snicar_data->asm_prm_bc1, snicar_data->ext_cff_mss_bc1, snicar_data->ss_alb_bc2, snicar_data->asm_prm_bc2, snicar_data->ext_cff_mss_bc2, 
      snicar_data->bcenh, mss_cnc_aer_in_fdb[idx], 
      g_star[idx], omega_star[idx], tau_star[idx]);


    ELM::snow_snicar::snow_radiative_transfer_solver(Land.urbpoi, flg_slr_in,flg_nosnl[idx], snl_top[idx], snl_btm[idx], 
      coszen[idx], h2osno[idx], mu_not[idx], flx_slrd_lcl[idx], flx_slri_lcl[idx], albsoi[idx], g_star[idx], 
      omega_star[idx], tau_star[idx], albout_lcl[idx], flx_abs_lcl[idx]);


    ELM::snow_snicar::snow_albedo_radiation_factor(Land.urbpoi, flg_slr_in, snl_top[idx], coszen[idx], mu_not[idx], 
      h2osno[idx], snw_rds_lcl[idx], albsoi[idx], albout_lcl[idx], flx_abs_lcl[idx], albsni[idx], 
      flx_absi_snw[idx]);
  }


  ELM::surface_albedo::ground_albedo(Land.urbpoi, coszen[idx], frac_sno[idx], albsod[idx], 
    albsoi[idx], albsnd[idx], albsni[idx], albgrd[idx], albgri[idx]);


  ELM::surface_albedo::flux_absorption_factor(Land, coszen[idx], frac_sno[idx], albsod[idx], albsoi[idx], 
    albsnd[idx], albsni[idx], flx_absd_snw[idx], flx_absi_snw[idx], flx_absdv[idx], flx_absdn[idx], 
    flx_absiv[idx], flx_absin[idx]);


  ELM::surface_albedo::canopy_layer_lai(Land.urbpoi, elai[idx], esai[idx], tlai[idx], tsai[idx], 
    nrad[idx], ncan[idx], tlai_z[idx], tsai_z[idx], fsun_z[idx], fabd_sun_z[idx], fabd_sha_z[idx], 
    fabi_sun_z[idx], fabi_sha_z[idx]);


  ELM::surface_albedo::two_stream_solver(Land, nrad[idx], coszen[idx], t_veg[idx], fwet[idx], elai[idx], 
    esai[idx], tlai_z[idx], tsai_z[idx], albgrd[idx], albgri[idx], albveg, vcmaxcintsun[idx], 
    vcmaxcintsha[idx], albd[idx], ftid[idx], ftdd[idx], fabd[idx], fabd_sun[idx], fabd_sha[idx], 
    albi[idx], ftii[idx], fabi[idx], fabi_sun[idx], fabi_sha[idx], fsun_z[idx], fabd_sun_z[idx], 
    fabd_sha_z[idx], fabi_sun_z[idx], fabi_sha_z[idx]);

}


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call canopy_hydrology kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{
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
//for (auto i : h2ocan) std::cout << "h2ocan:  " << i << std::endl;
//for (auto i : qflx_prec_grnd) std::cout << "qflx_prec_grnd:  " << i << std::endl;
//for (auto i : qflx_snwcp_liq) std::cout << "qflx_snwcp_liq:  " << i << std::endl;
//for (auto i : qflx_snwcp_ice) std::cout << "qflx_snwcp_ice:  " << i << std::endl;
//for (auto i : qflx_snow_grnd) std::cout << "qflx_snow_grnd:  " << i << std::endl;
//for (auto i : qflx_rain_grnd) std::cout << "qflx_rain_grnd:  " << i << std::endl;
//for (auto i : fwet) std::cout << "fwet:  " << i << std::endl;
//for (auto i : fdry) std::cout << "fdry:  " << i << std::endl;
//for (auto i : snow_depth) std::cout << "snow_depth:  " << i << std::endl;
//for (auto i : h2osno) std::cout << "h2osno:  " << i << std::endl;
//for (auto i : int_snow) std::cout << "int_snow:  " << i << std::endl;
//for (auto i : swe_old) std::cout << "swe_old:  " << i << std::endl;
//for (auto i : h2osoi_liq) std::cout << "h2osoi_liq:  " << i << std::endl;
//for (auto i : h2osoi_ice) std::cout << "h2osoi_ice:  " << i << std::endl;
//for (auto i : t_soisno) std::cout << "t_soisno:  " << i << std::endl;
//for (auto i : frac_iceold) std::cout << "frac_iceold:  " << i << std::endl;
//for (auto i : snl) std::cout << "snl:  " << i << std::endl;
//for (auto i : dz) std::cout << "dz:  " << i << std::endl;
//for (auto i : zsoi) std::cout << "zsoi:  " << i << std::endl;
//for (auto i : zisoi) std::cout << "zisoi:  " << i << std::endl;
//for (auto i : snw_rds) std::cout << "snw_rds:  " << i << std::endl;
//for (auto i : frac_sno_eff) std::cout << "frac_sno_eff:  " << i << std::endl;
//for (auto i : frac_sno) std::cout << "frac_sno:  " << i << std::endl;
//for (auto i : h2osfc) std::cout << "h2osfc:  " << i << std::endl;
//for (auto i : frac_h2osfc) std::cout << "frac_h2osfc:  " << i << std::endl;
  //for (auto i : h2ocan) std::cout << "h2ocan: " << i << std::endl;
  //for (auto i : fwet) std::cout << "fwet: " << i << std::endl;
  //for (auto i : fdry) std::cout << "fdry: " << i << std::endl;
  //for (auto i : qflx_prec_grnd) std::cout << "qflx_prec_grnd: " << i << std::endl;
  //for (auto i : qflx_snow_grnd) std::cout << "qflx_snow_grnd: " << i << std::endl;
  //for (auto i : qflx_rain_grnd) std::cout << "qflx_rain_grnd: " << i << std::endl;



/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call surface_radiation kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{
  // local to these kernel calls
    double trd[ELM::numrad] = {0.0,0.0};
    double tri[ELM::numrad] = {0.0,0.0};

    // call canopy_sunshade_fractions kernel
    ELM::surface_radiation::canopy_sunshade_fractions(
        Land, nrad[idx], elai[idx], tlai_z[idx], fsun_z[idx],
        forc_solad[idx], forc_solai[idx],
        fabd_sun_z[idx], fabd_sha_z[idx],
        fabi_sun_z[idx], fabi_sha_z[idx],
        parsun_z[idx], parsha_z[idx],
        laisun_z[idx], laisha_z[idx], laisun[idx], laisha[idx]);

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




/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call canopy_temperature kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{
  double qred; // soil surface relative humidity
  double hr;   // relative humidity

    ELM::canopy_temperature::old_ground_temp(Land, t_h2osfc[idx], t_soisno[idx], t_h2osfc_bef[idx], tssbef[idx]);

    ELM::canopy_temperature::ground_temp(Land, snl[idx], frac_sno_eff[idx], frac_h2osfc[idx], t_h2osfc[idx], t_soisno[idx],
                             t_grnd[idx]);

    ELM::canopy_temperature::calc_soilalpha(Land, frac_sno[idx], frac_h2osfc[idx], ELM::smpmin, h2osoi_liq[idx], h2osoi_ice[idx],
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






  
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call bareground_fluxes kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{
  int fake_frac_veg_nosno = frac_veg_nosno(0);

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


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call canopy_fluxes kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{
  



//for (auto i : t_soisno) std::cout << "TEST:: t_soisno: " << i << std::endl;
//for (auto i : tsoi_test) std::cout << "TEST:: tsoi_test: " << i << std::endl;

  const ELM::PSNVegData psnveg = vegdata.get_pft_psnveg(Land.vtype);
  // temporary data to pass between functions
  double wtg = 0.0;         // heat conductance for ground [m/s]
  double wtgq = 0.0;        // latent heat conductance for ground [m/s]
  double wtalq = 0.0;       // normalized latent heat cond. for air and leaf [-]
  double wtlq0 = 0.0;       // normalized latent heat conductance for leaf [-]
  double wtaq0 = 0.0;       // normalized latent heat conductance for air [-]
  double wtl0 = 0.0;        // normalized heat conductance for leaf [-]
  double wta0 = 0.0;        // normalized heat conductance for air [-]
  double wtal = 0.0;        // normalized heat conductance for air and leaf [-]
  double dayl_factor = 0.0; // scalar (0-1) for daylength effect on Vcmax
  double air = 0.0;         // atmos. radiation temporay set
  double bir = 0.0;         // atmos. radiation temporay set
  double cir = 0.0;         // atmos. radiation temporay set
  double el = 0.0;          // vapor pressure on leaf surface [pa]
  double qsatl = 0.0;       // leaf specific humidity [kg/kg]
  double qsatldT = 0.0;     // derivative of "qsatl" on "t_veg"
  double taf = 0.0;         // air temperature within canopy space [K]
  double qaf = 0.0;         // humidity of canopy air [kg/kg]
  double um = 0.0;          // wind speed including the stablity effect [m/s]
  double ur = 0.0;          // wind speed at reference height [m/s]
  double dth = 0.0;         // diff of virtual temp. between ref. height and surface
  double dqh = 0.0;         // diff of humidity between ref. height and surface
  double obu = 0.0;         // Monin-Obukhov length (m)
  double zldis = 0.0;       // reference height "minus" zero displacement height [m]
  double temp1 = 0.0;       // relation for potential temperature profile
  double temp2 = 0.0;       // relation for specific humidity profile
  double temp12m = 0.0;     // relation for potential temperature profile applied at 2-m
  double temp22m = 0.0;     // relation for specific humidity profile applied at 2-m
  double tlbef = 0.0;       // leaf temperature from previous iteration [K]
  double delq = 0.0;        // temporary
  double dt_veg = 0.0;      // change in t_veg, last iteration (Kelvin)

    ELM::canopy_fluxes::initialize_flux(
        Land, snl[idx], frac_veg_nosno[idx], frac_sno[idx], forc_hgt_u[idx],
        thm[idx], thv[idx], max_dayl, dayl, altmax_indx[idx], altmax_lastyear_indx[idx], 
        t_soisno[idx], h2osoi_ice[idx], h2osoi_liq[idx], dz[idx], rootfr[idx], psnveg.tc_stress, 
        sucsat[idx], watsat[idx], bsw[idx], psnveg.smpso, psnveg.smpsc, elai[idx], esai[idx], 
        emv[idx], emg[idx], qg[idx], t_grnd[idx], forc_tbot[idx], forc_pbot[idx], forc_lwrad[idx], 
        forc_u[idx], forc_v[idx], forc_qbot[idx], forc_thbot[idx], z0mg[idx], btran[idx], displa[idx], 
        z0mv[idx], z0hv[idx], z0qv[idx], rootr[idx], eff_porosity[idx], dayl_factor, air, bir, 
        cir, el, qsatl, qsatldT, taf, qaf, um, ur, obu, zldis, delq, t_veg[idx]);

    ELM::canopy_fluxes::stability_iteration(
        Land, dtime, snl[idx], frac_veg_nosno[idx], frac_sno[idx], forc_hgt_u[idx], 
        forc_hgt_t[idx], forc_hgt_q[idx], fwet[idx], fdry[idx], laisun[idx], 
        laisha[idx], forc_rho[idx], snow_depth[idx], soilbeta[idx], frac_h2osfc[idx], 
        t_h2osfc[idx], sabv[idx], h2ocan[idx], htop[idx],t_soisno[idx], air, bir, cir, 
        ur, zldis, displa[idx], elai[idx], esai[idx], t_grnd[idx],forc_pbot[idx], 
        forc_qbot[idx], forc_thbot[idx], z0mg[idx], z0mv[idx], z0hv[idx], z0qv[idx], thm[idx], 
        thv[idx], qg[idx], psnveg, nrad[idx], t10[idx], tlai_z[idx], vcmaxcintsha[idx], 
        vcmaxcintsun[idx], parsha_z[idx], parsun_z[idx], laisha_z[idx], laisun_z[idx], 
        forc_pco2[idx], forc_po2[idx], dayl_factor, btran[idx], qflx_tran_veg[idx], 
        qflx_evap_veg[idx], eflx_sh_veg[idx], wtg, wtl0, wta0, wtal, el, qsatl, qsatldT, 
        taf, qaf, um, dth, dqh, obu, temp1, temp2, temp12m, temp22m, tlbef, delq, dt_veg, 
        t_veg[idx], wtgq, wtalq, wtlq0, wtaq0);

    ELM::canopy_fluxes::compute_flux(
        Land, dtime, snl[idx], frac_veg_nosno[idx], frac_sno[idx], t_soisno[idx], frac_h2osfc[idx], 
        t_h2osfc[idx], sabv[idx], qg_snow[idx], qg_soil[idx], qg_h2osfc[idx], dqgdT[idx], htvp[idx], 
        wtg, wtl0, wta0, wtal, air, bir, cir, qsatl, qsatldT, dth, dqh, temp1, temp2, temp12m, 
        temp22m, tlbef, delq, dt_veg, t_veg[idx], t_grnd[idx], forc_pbot[idx], qflx_tran_veg[idx], 
        qflx_evap_veg[idx], eflx_sh_veg[idx], forc_qbot[idx], forc_rho[idx], thm[idx], emv[idx], 
        emg[idx], forc_lwrad[idx], wtgq, wtalq, wtlq0, wtaq0, h2ocan[idx], eflx_sh_grnd[idx], 
        eflx_sh_snow[idx], eflx_sh_soil[idx], eflx_sh_h2osfc[idx], qflx_evap_soi[idx], 
        qflx_ev_snow[idx], qflx_ev_soil[idx], qflx_ev_h2osfc[idx], dlrad[idx], ulrad[idx], 
        cgrnds[idx], cgrndl[idx], cgrnd[idx], t_ref2m[idx], t_ref2m_r[idx], q_ref2m[idx], 
        rh_ref2m[idx], rh_ref2m_r[idx]);
}



/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// call surface_fluxes kernels
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
{


ELM::surface_fluxes::initial_flux_calc(Land.urbpoi, snl[idx], frac_sno_eff[idx], frac_h2osfc[idx], t_h2osfc_bef[idx], tssbef(idx, ELM::nlevsno-snl(idx)), tssbef(idx, ELM::nlevsno),
  t_grnd[idx], cgrnds[idx], cgrndl[idx], eflx_sh_grnd[idx], qflx_evap_soi[idx], qflx_ev_snow[idx], qflx_ev_soil[idx],
  qflx_ev_h2osfc[idx]);



ELM::surface_fluxes::update_surface_fluxes(Land.urbpoi, do_capsnow, snl[idx], dtime, t_grnd[idx], htvp[idx], frac_sno_eff[idx], frac_h2osfc[idx], t_h2osfc_bef[idx], sabg_soil[idx],
  sabg_snow[idx], dlrad[idx], frac_veg_nosno[idx], emg[idx], forc_lwrad[idx], tssbef(idx, ELM::nlevsno-snl(idx)), tssbef(idx, ELM::nlevsno),
  h2osoi_ice(idx, ELM::nlevsno-snl(idx)), h2osoi_liq(idx, ELM::nlevsno), 
eflx_sh_veg[idx],
qflx_evap_veg[idx], qflx_evap_soi[idx], eflx_sh_grnd[idx], qflx_ev_snow[idx], qflx_ev_soil[idx], qflx_ev_h2osfc[idx],
  eflx_soil_grnd[idx], eflx_sh_tot[idx], qflx_evap_tot[idx], eflx_lh_tot[idx], qflx_evap_grnd[idx], qflx_sub_snow[idx], qflx_dew_snow[idx], qflx_dew_grnd[idx], qflx_snwcp_liq[idx],
  qflx_snwcp_ice[idx]);


}


//for (auto i : eflx_sh_snow) std::cout << "eflx_sh_snow: " << i << std::endl;
//for (auto i : eflx_sh_soil) std::cout << "eflx_sh_soil: " << i << std::endl;
//for (auto i : qflx_ev_snow) std::cout << "qflx_ev_snow: " << i << std::endl;
//for (auto i : qflx_ev_soil) std::cout << "qflx_ev_soil: " << i << std::endl;
//for (auto i : eflx_sh_grnd) std::cout << "eflx_sh_grnd: " << i << std::endl;
//for (auto i : qflx_evap_soi) std::cout << "qflx_evap_soi: " << i << std::endl;
//for (auto i : forc_rho) std::cout << "forc_rho: " << i << std::endl;
//  for (auto i : forc_qbot) std::cout << "forc_qbot: " << i << std::endl;


//std::cout << "forc_tbot: " << forc_tbot(0) << std::endl;
//std::cout << "forc_thbot: " << forc_thbot(0) << std::endl;
//std::cout << "forc_pbot: " << forc_pbot(0) << std::endl;
//std::cout << "forc_qbot: " << forc_qbot(0) << std::endl;
//std::cout << "forc_rh: " << forc_rh(0) << std::endl;
//std::cout << "forc_lwrad: " << forc_lwrad(0) << std::endl;
//std::cout << "forc_rain: " << forc_rain(0) << std::endl;
//std::cout << "forc_snow: " << forc_snow(0) << std::endl;
//std::cout << "forc_u: " << forc_u(0) << std::endl;
//std::cout << "forc_v: " << forc_v(0) << std::endl;
//std::cout << "forc_hgt: " << forc_hgt(0) << std::endl;
//std::cout << "forc_hgt_u: " << forc_hgt_u(0) << std::endl;
//std::cout << "forc_hgt_t: " << forc_hgt_t(0) << std::endl;
//std::cout << "forc_hgt_q: " << forc_hgt_q(0) << std::endl;
//std::cout << "forc_vp: " << forc_vp(0) << std::endl;
//std::cout << "forc_rho: " << forc_rho(0) << std::endl;
//std::cout << "forc_po2: " << forc_po2(0) << std::endl;
//std::cout << "forc_pco2: " << forc_pco2(0) << std::endl;

//for (auto i : eflx_sh_veg) std::cout << "eflx_sh_veg: " << i << std::endl;
//for (auto i : eflx_sh_snow) std::cout << "eflx_sh_snow: " << i << std::endl;
//for (auto i : eflx_sh_soil) std::cout << "eflx_sh_soil: " << i << std::endl;
//for (auto i : qflx_ev_snow) std::cout << "qflx_ev_snow: " << i << std::endl;
//for (auto i : qflx_ev_soil) std::cout << "qflx_ev_soil: " << i << std::endl;
//for (auto i : eflx_sh_grnd) std::cout << "eflx_sh_grnd: " << i << std::endl;
//for (auto i : qflx_evap_soi) std::cout << "qflx_evap_soi: " << i << std::endl;
//for (auto i : qflx_evap_veg) std::cout << "qflx_evap_veg: " << i << std::endl;
//for (auto i : qflx_tran_veg) std::cout << "qflx_tran_veg: " << i << std::endl;
//for (auto i : btran) std::cout << "btran: " << i << std::endl;
//for (auto i : t_soisno) std::cout << "t_soisno:  " << i << std::endl;
//for (auto i : rootfr) std::cout << "rootfr:  " << i << std::endl;


for (auto i : snl) std::cout << std::setprecision(10) << "snl:  " << i << std::endl;
for (auto i : frac_veg_nosno) std::cout << "frac_veg_nosno:  " << i << std::endl;
for (auto i : nrad) std::cout << "nrad:  " << i << std::endl;
for (auto i : altmax_indx) std::cout << "altmax_indx:  " << i << std::endl;
for (auto i : altmax_lastyear_indx) std::cout << "altmax_lastyear_indx:  " << i << std::endl;
for (auto i : frac_sno) std::cout << "frac_sno:  " << i << std::endl;
for (auto i : forc_hgt_u) std::cout << "forc_hgt_u:  " << i << std::endl;
for (auto i : thm) std::cout << "thm:  " << i << std::endl;
for (auto i : thv) std::cout << "thv:  " << i << std::endl;
std::cout << "max_dayl:  " << max_dayl << std::endl;
std::cout << "dayl:  " << dayl << std::endl;
for (auto i : elai) std::cout << "elai:  " << i << std::endl;
for (auto i : esai) std::cout << "esai:  " << i << std::endl;
for (auto i : emv) std::cout << "emv:  " << i << std::endl;
for (auto i : emg) std::cout << "emg:  " << i << std::endl;
for (auto i : qg) std::cout << "qg:  " << i << std::endl;
for (auto i : t_grnd) std::cout << "t_grnd:  " << i << std::endl;
for (auto i : forc_tbot) std::cout << "forc_t:  " << i << std::endl;
for (auto i : forc_pbot) std::cout << "forc_pbot:  " << i << std::endl;
for (auto i : forc_lwrad) std::cout << "forc_lwrad:  " << i << std::endl;
for (auto i : forc_u) std::cout << "forc_u:  " << i << std::endl;
for (auto i : forc_v) std::cout << "forc_v:  " << i << std::endl;
for (auto i : forc_qbot) std::cout << "forc_q:  " << i << std::endl;
for (auto i : forc_thbot) std::cout << "forc_th:  " << i << std::endl;
for (auto i : z0mg) std::cout << "z0mg:  " << i << std::endl;
for (auto i : btran) std::cout << "btran:  " << i << std::endl;
for (auto i : displa) std::cout << "displa:  " << i << std::endl;
for (auto i : z0mv) std::cout << "z0mv:  " << i << std::endl;
for (auto i : z0hv) std::cout << "z0hv:  " << i << std::endl;
for (auto i : z0qv) std::cout << "z0qv:  " << i << std::endl;
for (auto i : t_veg) std::cout << "t_veg:  " << i << std::endl;
for (auto i : forc_hgt_t) std::cout << "forc_hgt_t:  " << i << std::endl;
for (auto i : forc_hgt_q) std::cout << "forc_hgt_q:  " << i << std::endl;
for (auto i : fwet) std::cout << std::setprecision(10) << "fwet:  " << i << std::endl;
for (auto i : fdry) std::cout << "fdry:  " << i << std::endl;
for (auto i : laisun) std::cout << "laisun:  " << i << std::endl;
for (auto i : laisha) std::cout << "laisha:  " << i << std::endl;
for (auto i : forc_rho) std::cout << "forc_rho:  " << i << std::endl;
for (auto i : snow_depth) std::cout << "snow_depth:  " << i << std::endl;
for (auto i : soilbeta) std::cout << "soilbeta:  " << i << std::endl;
for (auto i : frac_h2osfc) std::cout << "frac_h2osfc:  " << i << std::endl;
for (auto i : t_h2osfc) std::cout << "t_h2osfc:  " << i << std::endl;
for (auto i : sabv) std::cout << "sabv:  " << i << std::endl;
for (auto i : h2ocan) std::cout << "h2ocan:  " << i << std::endl;
for (auto i : htop) std::cout << "htop:  " << i << std::endl;
for (auto i : t10) std::cout << "t10:  " << i << std::endl;
for (auto i : vcmaxcintsha) std::cout << "vcmaxcintsha:  " << i << std::endl;
for (auto i : vcmaxcintsun) std::cout << "vcmaxcintsun:  " << i << std::endl;
for (auto i : forc_pco2) std::cout << "forc_pco2:  " << i << std::endl;
for (auto i : forc_po2) std::cout << "forc_po2:  " << i << std::endl;
for (auto i : qflx_tran_veg) std::cout << "qflx_tran_veg:  " << i << std::endl;
for (auto i : qflx_evap_veg) std::cout << "qflx_evap_veg:  " << i << std::endl;
for (auto i : eflx_sh_veg) std::cout << "eflx_sh_veg:  " << i << std::endl;
for (auto i : qg_snow) std::cout << "qg_snow:  " << i << std::endl;
for (auto i : qg_soil) std::cout << "qg_soil:  " << i << std::endl;
for (auto i : qg_h2osfc) std::cout << "qg_h2osfc:  " << i << std::endl;
for (auto i : dqgdT) std::cout << "dqgdT:  " << i << std::endl;
for (auto i : htvp) std::cout << "htvp:  " << i << std::endl;
for (auto i : eflx_sh_grnd) std::cout << "eflx_sh_grnd:  " << i << std::endl;
for (auto i : eflx_sh_snow) std::cout << "eflx_sh_snow:  " << i << std::endl;
for (auto i : eflx_sh_soil) std::cout << "eflx_sh_soil:  " << i << std::endl;
for (auto i : eflx_sh_h2osfc) std::cout << "eflx_sh_h2osfc:  " << i << std::endl;
for (auto i : qflx_evap_soi) std::cout << "qflx_evap_soi:  " << i << std::endl;
for (auto i : qflx_ev_snow) std::cout << "qflx_ev_snow:  " << i << std::endl;
for (auto i : qflx_ev_soil) std::cout << "qflx_ev_soil:  " << i << std::endl;
for (auto i : qflx_ev_h2osfc) std::cout << "qflx_ev_h2osfc:  " << i << std::endl;
for (auto i : dlrad) std::cout << "dlrad:  " << i << std::endl;
for (auto i : ulrad) std::cout << "ulrad:  " << i << std::endl;
for (auto i : cgrnds) std::cout << "cgrnds:  " << i << std::endl;
for (auto i : cgrndl) std::cout << "cgrndl:  " << i << std::endl;
for (auto i : cgrnd) std::cout << "cgrnd:  " << i << std::endl;
for (auto i : t_ref2m) std::cout << "t_ref2m:  " << i << std::endl;
for (auto i : t_ref2m_r) std::cout << "t_ref2m_r:  " << i << std::endl;
for (auto i : q_ref2m) std::cout << "q_ref2m:  " << i << std::endl;
for (auto i : rh_ref2m) std::cout << "rh_ref2m:  " << i << std::endl;
for (auto i : rh_ref2m_r) std::cout << "rh_ref2m_r:  " << i << std::endl;
for (auto i : rootr) std::cout << "rootr:  " << i << std::endl;
for (auto i : eff_porosity) std::cout << "eff_porosity:  " << i << std::endl;
for (auto i : tlai_z) std::cout << "tlai_z:  " << i << std::endl;
for (auto i : parsha_z) std::cout << "parsha_z:  " << i << std::endl;
for (auto i : parsun_z) std::cout << "parsun_z:  " << i << std::endl;
for (auto i : laisha_z) std::cout << "laisha_z:  " << i << std::endl;
for (auto i : laisun_z) std::cout << "laisun_z:  " << i << std::endl;
for (auto i : t_soisno) std::cout << "t_soisno:  " << i << std::endl;
for (auto i : h2osoi_ice) std::cout << "h2osoi_ice:  " << i << std::endl;
for (auto i : h2osoi_liq) std::cout << "h2osoi_liq:  " << i << std::endl;
for (auto i : dz) std::cout << "dz:  " << i << std::endl;
for (auto i : rootfr) std::cout << "rootfr:  " << i << std::endl;
for (auto i : sucsat) std::cout << "sucsat:  " << i << std::endl;
for (auto i : watsat) std::cout << "watsat:  " << i << std::endl;
for (auto i : bsw) std::cout << "bsw:  " << i << std::endl;

for (auto i : laisun) std::cout << "laisun:  " << i << std::endl;
for (auto i : laisha) std::cout << "laisha:  " << i << std::endl;
for (auto i : tlai_z) std::cout << "tlai_z:  " << i << std::endl;
for (auto i : fsun_z) std::cout << "fsun_z:  " << i << std::endl;
for (auto i : fabd_sun_z) std::cout << "fabd_sun_z:  " << i << std::endl;
for (auto i : fabd_sha_z) std::cout << "fabd_sha_z:  " << i << std::endl;
for (auto i : fabi_sun_z) std::cout << "fabi_sun_z:  " << i << std::endl;
for (auto i : fabi_sha_z) std::cout << "fabi_sha_z:  " << i << std::endl;
for (auto i : parsun_z) std::cout << "parsun_z:  " << i << std::endl;
for (auto i : parsha_z) std::cout << "parsha_z:  " << i << std::endl;
for (auto i : laisun_z) std::cout << "laisun_z:  " << i << std::endl;
for (auto i : laisha_z) std::cout << "laisha_z:  " << i << std::endl;
for (auto i : forc_solad) std::cout << "forc_solad:  " << i << std::endl;
for (auto i : forc_solai) std::cout << "forc_solai:  " << i << std::endl;



std::cout << "END STEP   " << t << "  " << std::get<1>(current.date()) << 
"  " << std::get<2>(current.date()) << "  " << current.sec << std::endl;
current.increment_seconds(1800);

}



}

