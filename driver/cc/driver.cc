
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"
#include "vegproperties.h"
#include <iostream>
#include <string>

#include "BareGroundFluxes.h"
#include "ELMConstants.h"
#include "InitSnowLayers.hh"
#include "InitTimestep.hh"
#include "InitTopography.hh"
#include "LandType.h"
#include "ReadAtmosphere.hh"
#include "ReadPFTConstants.hh"
#include "ReadTestData.hh"
#include "SatellitePhenology.hh"

#include "CanopyHydrology.h"
#include "CanopyTemperature.h"
#include "SurfaceRadiation.h"

using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  int MPI_COMM_WORLD;
  const int n_procs = 1;
  const int myrank = 0;
  auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
  std::cout << " proc_decomp: " << proc_decomp[0] << "," << proc_decomp[1] << std::endl;

  const auto start = ELM::Utils::Date(2012, 7, 1);
  const int n_months = 1;
  const int ticks_per_day = 48; // 0.5-hourly data
  const int write_interval = 1;
  int dtime = 1800;

  std::cout << start.sec << "," << start.doy << std::endl;
  ELM::Utils::Date current_start(start);

  // input data
  const std::string data_dir("/home/jbeisman/Software/elm_test_input");
  const std::string basename_forc("");
  const std::string basename_phen("surfdata_1x1pt_US-Brw_simyr1850_c360x720_c20190221.nc");
  const std::string basename_init("ForATS_AK-BEOG_ICB1850CNPRDCTCBC.elm.r.0609-01-01-00000.nc");
  const std::string basename_pfts("clm_params.nc");

  const auto forc_dims =
      ELM::IO::get_forcing_dimensions(MPI_COMM_WORLD, data_dir, basename_forc, "PRECTmms", start, n_months);
  std::cout << " forc_dims: " << forc_dims[0] << "," << forc_dims[1] << "," << forc_dims[2] << std::endl;

  const auto phen_dims =
      ELM::IO::get_phenology_dimensions(MPI_COMM_WORLD, data_dir, basename_phen, "MONTHLY_LAI", start, n_months);
  std::cout << " phen_dims: " << phen_dims[0] << "," << phen_dims[1] << "," << phen_dims[2] << "," << phen_dims[3]
            << std::endl;

  auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp, {forc_dims[1], forc_dims[2]},
                                                       {myrank / proc_decomp[1], myrank % proc_decomp[1]});

  auto n_grid_cells = dd.n_local[0] * dd.n_local[1];
  auto ntimes = forc_dims[0];
  auto n_pfts = phen_dims[1];

  // the grid -- need to init from input or driver data; so far, only above ground snow layers get vals
  auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto z = create<ArrayD2>("z", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto zi = create<ArrayD2>("zi", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd + 1);

  auto timwt = create<ArrayD1>("timwt", 2);            // time weights for lai/sai/hgt
  static const double dewmx = 0.1;                     // hardwired
  auto vtype = create<ArrayI1>("vtype", n_grid_cells); // pft type
  assign(vtype, 2);                                    // dummy variable for now
  assign(dz, 1.0);                                     // hardwired for simple testing

  std::cout << vtype[0] << std::endl;

  // increment seconds based timestep:
  // ELM::Utils::Date current_test(current_start);
  // current_test.increment_seconds(500);
  // std::cout << current_start.sec << "," << current_test.sec << std::endl;

  if (myrank == n_procs - 1) {
    std::cout << " dimensions (rank " << myrank << "):" << std::endl
              << "   n_time = " << forc_dims[0] << std::endl
              << "   global_n_lat,lon = " << dd.n_global[0] << "," << dd.n_global[1] << std::endl
              << "   nprocs_lat,lon = " << dd.n_procs[0] << "," << dd.n_procs[1] << std::endl
              << "   local_start = " << dd.start[0] << "," << dd.start[1] << std::endl
              << "   local_n_lat,lon = " << dd.n_local[0] << "," << dd.n_local[1] << std::endl;
  }

  // PFT variables
  const auto maxpfts = ELM::IO::get_maxpfts(MPI_COMM_WORLD, data_dir, basename_pfts, "z0mr");
  std::cout << " maxpfts: " << maxpfts[0] << std::endl;
  auto pftnames = create<ArrayS1>("pftnames", maxpfts[0]);
  auto z0mr = create<ArrayD1>("pftnames", maxpfts[0]);
  auto displar = create<ArrayD1>("displar", maxpfts[0]);
  auto dleaf = create<ArrayD1>("dleaf", maxpfts[0]);
  auto c3psn = create<ArrayD1>("c3psn", maxpfts[0]);
  auto xl = create<ArrayD1>("xl", maxpfts[0]);
  auto roota_par = create<ArrayD1>("roota_par", maxpfts[0]);
  auto rootb_par = create<ArrayD1>("rootb_par", maxpfts[0]);
  auto slatop = create<ArrayD1>("slatop", maxpfts[0]);
  auto leafcn = create<ArrayD1>("leafcn", maxpfts[0]);
  auto flnr = create<ArrayD1>("flnr", maxpfts[0]);
  auto smpso = create<ArrayD1>("smpso", maxpfts[0]);
  auto smpsc = create<ArrayD1>("smpsc", maxpfts[0]);
  auto fnitr = create<ArrayD1>("fnitr", maxpfts[0]);
  auto fnr = create<ArrayD1>("fnr", maxpfts[0]);
  auto act25 = create<ArrayD1>("act25", maxpfts[0]);
  auto kcha = create<ArrayD1>("kcha", maxpfts[0]);
  auto koha = create<ArrayD1>("koha", maxpfts[0]);
  auto cpha = create<ArrayD1>("cpha", maxpfts[0]);
  auto vcmaxha = create<ArrayD1>("vcmaxha", maxpfts[0]);
  auto jmaxha = create<ArrayD1>("jmaxha", maxpfts[0]);
  auto tpuha = create<ArrayD1>("tpuha", maxpfts[0]);
  auto lmrha = create<ArrayD1>("lmrha", maxpfts[0]);
  auto vcmaxhd = create<ArrayD1>("vcmaxhd", maxpfts[0]);
  auto jmaxhd = create<ArrayD1>("jmaxhd", maxpfts[0]);
  auto tpuhd = create<ArrayD1>("tpuhd", maxpfts[0]);
  auto lmrhd = create<ArrayD1>("lmrhd", maxpfts[0]);
  auto lmrse = create<ArrayD1>("lmrse", maxpfts[0]);
  auto qe = create<ArrayD1>("qe", maxpfts[0]);
  auto theta_cj = create<ArrayD1>("theta_cj", maxpfts[0]);
  auto bbbopt = create<ArrayD1>("bbbopt", maxpfts[0]);
  auto mbbopt = create<ArrayD1>("mbbopt", maxpfts[0]);
  auto nstor = create<ArrayD1>("nstor", maxpfts[0]);
  auto br_xr = create<ArrayD1>("br_xr", maxpfts[0]);
  auto tc_stress =
      create<ArrayD1>("tc_stress", 1); // only one value - keep in container for compatibility with NetCDF reader
  auto rhol = create<ArrayD2>("rhol", 2, maxpfts[0]); // numrad
  auto rhos = create<ArrayD2>("rhol", 2, maxpfts[0]); // numrad
  auto taul = create<ArrayD2>("rhol", 2, maxpfts[0]); // numrad
  auto taus = create<ArrayD2>("rhol", 2, maxpfts[0]); // numrad

  // read pft data from clm_params.nc
  ELM::ReadPFTConstants(data_dir, basename_pfts, pftnames, z0mr, displar, dleaf, c3psn, xl, roota_par, rootb_par,
                        slatop, leafcn, flnr, smpso, smpsc, fnitr, fnr, act25, kcha, koha, cpha, vcmaxha, jmaxha, tpuha,
                        lmrha, vcmaxhd, jmaxhd, tpuhd, lmrhd, lmrse, qe, theta_cj, bbbopt, mbbopt, nstor, br_xr,
                        tc_stress, rhol[0], rhol[1], rhos[0], rhos[1], taul[0], taul[1], taus[0], taus[1]);

  for (int j = 0; j != maxpfts[0]; j++) {
    std::cout << "pftname: " << pftnames[j] << std::endl;
    std::cout << "rhol[0], rhol[1]: " << rhol[0][j] << " " << rhol[1][j] << std::endl;
    std::cout << "rhos[0], rhos[1]: " << rhos[0][j] << " " << rhos[1][j] << std::endl;
    std::cout << "taul[0], taul[1]: " << taul[0][j] << " " << taul[1][j] << std::endl;
    std::cout << "taus[0], taus[1]: " << taus[0][j] << " " << taus[1][j] << std::endl;
  }

  auto snl = create<ArrayI1>("snl", n_grid_cells);
  auto topo_slope = create<ArrayD1>("topo_slope", n_grid_cells);
  auto topo_std = create<ArrayD1>("topo_std", n_grid_cells);
  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells); // NEED VALUES! -- use init nc file for now
  auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);     // NEED VALUES! -- use init nc file for now

  auto h2ocan = create<ArrayD1>("h2ocan", n_grid_cells);
  auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd); // need value
  auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd); // need value
  auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);                                       // need value
  auto snw_rds = create<ArrayD2>("snw_rds", n_grid_cells, ELM::nlevsno);
  auto int_snow = create<ArrayD1>("int_snow", n_grid_cells); // need value

  auto sand3d = create<ArrayD2>("sand3d", n_grid_cells, ELM::nlevsoi); // should change name from sand3d, clay3d
  auto clay3d = create<ArrayD2>("clay3d", n_grid_cells, ELM::nlevsoi);
  ELM::ReadLandData(data_dir, basename_phen, dd, topo_slope, topo_std, sand3d, clay3d);
  std::cout << topo_slope[0] << " " << topo_std[0] << std::endl;
  for (int i = 0; i < ELM::nlevsoi; ++i) {
    std::cout << "sand3d & clay3d: " << i << " " << sand3d[0][i] << " " << clay3d[0][i] << std::endl;
  }

  ELM::ReadTestInitData(data_dir, basename_init, dd, t_soisno, snl, snow_depth, h2ocan, h2osoi_liq, h2osoi_ice, h2osno,
                        snw_rds, int_snow);

  for (int i = 0; i < 1; i++) {
    std::cout << snl[i] << std::endl;
    std::cout << snow_depth[i] << std::endl;
    std::cout << i << " " << h2ocan[i] << std::endl;
  }

  // allocate 2 months of LAI, SAI, HTOP, HBOT for interpolation
  auto mlai2t = create<ArrayD2>("mlai2t", n_grid_cells, 2);
  auto msai2t = create<ArrayD2>("msai2t", n_grid_cells, 2);
  auto mhvt2t = create<ArrayD2>("mhvt2t", n_grid_cells, 2);
  auto mhvb2t = create<ArrayD2>("mhvb2t", n_grid_cells, 2);

  // calculate monthly weights; read 2 months of LAI, SAI, HTOP, HBOT if required
  ELM::InterpMonthlyVeg(start, dtime, n_pfts, data_dir, basename_phen, dd, vtype, timwt, mlai2t, msai2t, mhvt2t,
                        mhvb2t);

  // prescribed sat phenology
  auto tlai = create<ArrayD1>("tlai", n_grid_cells);
  auto tsai = create<ArrayD1>("tsai", n_grid_cells);
  auto elai = create<ArrayD1>("elai", n_grid_cells);
  auto esai = create<ArrayD1>("esai", n_grid_cells);
  auto htop = create<ArrayD1>("htop", n_grid_cells);
  auto hbot = create<ArrayD1>("hbot", n_grid_cells);
  auto frac_veg_nosno_alb = create<ArrayI1>("frac_veg_nosno_alb", n_grid_cells);

  int idx = 0;
  ELM::SatellitePhenology(mlai2t[idx], msai2t[idx], mhvt2t[idx], mhvb2t[idx], timwt, vtype[idx], snow_depth[idx],
                          frac_sno[idx], tlai[idx], tsai[idx], htop[idx], hbot[idx], elai[idx], esai[idx],
                          frac_veg_nosno_alb[idx]);

  for (int mon = 0; mon < 2; ++mon) {
    for (int i = 0; i != n_grid_cells; i++) {
      std::cout << "mlai: " << i << " " << mon << " " << mlai2t(i, mon) << std::endl;
      std::cout << "msai: " << i << " " << mon << " " << msai2t(i, mon) << std::endl;
      std::cout << "mhgtt: " << i << " " << mon << " " << mhvt2t(i, mon) << std::endl;
      std::cout << "mhgtb: " << i << " " << mon << " " << mhvb2t(i, mon) << std::endl;
    }
  }

  std::cout << "::: " << tlai[idx] << " " << tsai[idx] << " " << htop[idx] << " " << hbot[idx] << " " << elai[idx]
            << " " << esai[idx] << std::endl;

  // raw ATM forcing data
  auto atm_zbot = create<ArrayD2>("atm_zbot", ntimes, n_grid_cells);
  auto atm_tbot = create<ArrayD2>("atm_tbot", ntimes, n_grid_cells);
  auto atm_rh = create<ArrayD2>("atm_rh", ntimes, n_grid_cells);
  auto atm_wind = create<ArrayD2>("atm_wind", ntimes, n_grid_cells);
  auto atm_fsds = create<ArrayD2>("atm_fsds", ntimes, n_grid_cells);
  auto atm_flds = create<ArrayD2>("atm_flds", ntimes, n_grid_cells);
  auto atm_psrf = create<ArrayD2>("atm_psrf", ntimes, n_grid_cells);
  auto atm_prec = create<ArrayD2>("atm_prec", ntimes, n_grid_cells);

  ELM::ReadAtmForcing(data_dir, basename_forc, start, dd, n_months, atm_zbot, atm_tbot, atm_rh, atm_wind, atm_fsds,
                      atm_flds, atm_psrf, atm_prec);

  //  for (int t = 0; t < 1488; ++t) {
  //    std::cout << "::: " << atm_zbot(t, 0) << " " << atm_tbot(t, 0) << " " << atm_rh(t, 0) << " " << atm_wind(t, 0)
  //              << " " << atm_fsds(t, 0) << " " << atm_flds(t, 0) << " " << atm_psrf(t, 0) << " " << atm_prec(t, 0)
  //              << std::endl;
  //  }

  // processed ATM forcing data
  auto forc_t = create<ArrayD1>("forc_t", n_grid_cells);
  auto forc_th = create<ArrayD1>("forc_th", n_grid_cells);
  auto forc_pbot = create<ArrayD1>("forc_pbot", n_grid_cells);
  auto forc_q = create<ArrayD1>("forc_q", n_grid_cells);
  auto forc_lwrad = create<ArrayD1>("forc_lwrad", n_grid_cells);
  auto forc_rain = create<ArrayD1>("forc_rain", n_grid_cells);
  auto forc_snow = create<ArrayD1>("forc_snow", n_grid_cells);
  auto forc_u = create<ArrayD1>("forc_u", n_grid_cells);
  auto forc_v = create<ArrayD1>("forc_v", n_grid_cells);
  auto forc_rh = create<ArrayD1>("forc_rh", n_grid_cells);
  auto forc_rho = create<ArrayD1>("forc_rho", n_grid_cells);
  auto forc_po2 = create<ArrayD1>("forc_po2", n_grid_cells);
  auto forc_pco2 = create<ArrayD1>("forc_pco2", n_grid_cells);
  auto forc_hgt_u = create<ArrayD1>("forc_hgt_u", n_grid_cells);
  auto forc_hgt_t = create<ArrayD1>("forc_hgt_t", n_grid_cells);
  auto forc_hgt_q = create<ArrayD1>("forc_hgt_q", n_grid_cells);
  auto forc_solad = create<ArrayD2>("forc_solad", n_grid_cells, 2);
  auto forc_solai = create<ArrayD2>("forc_solai", n_grid_cells, 2);

  for (int t = 0; t < 1488; ++t) {
    ELM::GetAtmTimestep(atm_tbot(t, idx), atm_psrf(t, idx), atm_rh(t, idx), atm_flds(t, idx), atm_fsds(t, idx),
                        atm_prec(t, idx), atm_wind(t, idx), forc_t[0], forc_th[0], forc_pbot[0], forc_q[0],
                        forc_lwrad[idx], forc_rain[idx], forc_snow[idx], forc_u[idx], forc_v[idx], forc_rh[idx],
                        forc_rho[idx], forc_po2[idx], forc_pco2[idx], forc_hgt_u[idx], forc_hgt_t[idx], forc_hgt_q[idx],
                        forc_solad[idx], forc_solai[idx]);

    //  std::cout << t << " ::: " << forc_lwrad[idx] << " " << forc_u[idx] << " " << forc_rh[idx] << " " << forc_t[idx]
    //  << " "
    //    << forc_th[idx] << " " << forc_rho[idx] << "forc_solad0: " << forc_solad[idx][0] << "forc_solad1: " <<
    //    forc_solad[idx][1] << "forc_solai0: "
    //    << forc_solai[idx][0] << "forc_solai1: " << forc_solai[idx][1] << " " << forc_rain[idx] << " " <<
    //    forc_snow[idx]
    //    << std::endl;
    // std::cout << "forc_rho: " << forc_rho[idx] << " " << forc_pbot[idx] << " " << forc_q[idx] << " " << forc_rh[idx]
    // << std::endl;
  }

  const bool lakpoi(false);
  bool do_capsnow(false);
  auto veg_active = create<ArrayB1>("veg_active", n_grid_cells); // need value
  assign(veg_active, true);                                      // hardwired
  auto h2osno_old = create<ArrayD1>("h2osno_old", n_grid_cells);
  auto eflx_bot = create<ArrayD1>("eflx_bot", n_grid_cells);
  auto qflx_glcice = create<ArrayD1>("qflx_glcice", n_grid_cells);

  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto frac_iceold = create<ArrayD2>("frac_iceold", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);

  ELM::InitSnowLayers(snow_depth[idx], lakpoi, snl[idx], dz[idx], z[idx], zi[idx]);
  ELM::InitTimestep(lakpoi, h2osno[idx], veg_active[idx], snl[idx], h2osoi_ice[idx], h2osoi_liq[idx],
                    frac_veg_nosno_alb[idx], h2osno_old[idx], do_capsnow, eflx_bot[idx], qflx_glcice[idx],
                    frac_veg_nosno[idx], frac_iceold[idx]);

  std::cout << snl[0] << std::endl;
  for (int i = 0; i < ELM::nlevsno + ELM::nlevgrnd; i++) {
    std::cout << dz[0][i] << std::endl;
  }

  // instantiate data
  ELM::LandType Land;
  ELM::VegProperties Veg;

  auto n_irrig_steps_left = create<ArrayI1>("n_irrig_steps_left", n_grid_cells);
  auto irrig_rate = create<ArrayD1>("irrig_rate", n_grid_cells);
  auto qflx_irrig = create<ArrayD1>("qflx_irrig", n_grid_cells);
  auto qflx_prec_grnd = create<ArrayD1>("qflx_prec_grnd", n_grid_cells);
  auto qflx_snwcp_liq = create<ArrayD1>("qflx_snwcp_liq", n_grid_cells);
  auto qflx_snwcp_ice = create<ArrayD1>("qflx_snwcp_ice", n_grid_cells);
  auto qflx_snow_grnd = create<ArrayD1>("qflx_snow_grnd", n_grid_cells);
  auto qflx_rain_grnd = create<ArrayD1>("qflx_rain_grnd", n_grid_cells);
  auto fwet = create<ArrayD1>("fwet", n_grid_cells);
  auto fdry = create<ArrayD1>("fdry", n_grid_cells);
  auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
  auto qflx_snow_melt = create<ArrayD1>("qflx_snow_melt", n_grid_cells);

  auto qflx_snow_h2osfc = create<ArrayD1>("qflx_snow_h2osfc", n_grid_cells);
  auto h2osfc = create<ArrayD1>("h2osfc", n_grid_cells);
  auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
  auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", n_grid_cells);
  auto swe_old = create<ArrayD2>("swe_old", n_grid_cells, ELM::nlevsno);

  ELM::InitTopoSlope(topo_slope[idx]);

  auto n_melt = create<ArrayD1>("n_melt", n_grid_cells);
  auto micro_sigma = create<ArrayD1>("micro_sigma", n_grid_cells);
  ELM::InitMicroTopo(Land.ltype, topo_slope[idx], topo_std[idx], n_melt[idx], micro_sigma[idx]);
  std::cout << "n_melt: " << n_melt[idx] << std::endl;

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

    ELM::Interception(Land, frac_veg_nosno[idx], forc_rain[idx], forc_snow[idx], dewmx, elai[idx], esai[idx], dtime,
                      h2ocan[idx], qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain);
    std::cout << frac_veg_nosno[idx] << " " << forc_rain[idx] << " " << forc_snow[idx] << " " << dewmx << " "
              << elai[idx] << " " << esai[idx] << " " << dtime << " " << h2ocan[idx] << " " << qflx_candrip << " "
              << qflx_through_snow << " " << qflx_through_rain << " " << fracsnow << " " << fracrain << std::endl;

    ELM::GroundFlux(Land, do_capsnow, frac_veg_nosno[idx], forc_rain[idx], forc_snow[idx], qflx_irrig[idx],
                    qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain, qflx_prec_grnd[idx],
                    qflx_snwcp_liq[idx], qflx_snwcp_ice[idx], qflx_snow_grnd[idx], qflx_rain_grnd[idx]);
    std::cout << qflx_prec_grnd[idx] << " " << qflx_snwcp_liq[idx] << " " << qflx_snwcp_ice[idx] << " "
              << qflx_snow_grnd[idx] << " " << qflx_rain_grnd[idx] << std::endl;

    ELM::FracWet(Land, frac_veg_nosno[idx], dewmx, elai[idx], esai[idx], h2ocan[idx], fwet[idx], fdry[idx]);
    std::cout << fwet[idx] << " " << fdry[idx] << std::endl;

    ELM::SnowInit(Land, dtime, do_capsnow, oldfflag, forc_t[idx], t_grnd[idx], qflx_snow_grnd[idx], qflx_snow_melt[idx],
                  n_melt[idx], snow_depth[idx], h2osno[idx], int_snow[idx], swe_old[idx], h2osoi_liq[idx],
                  h2osoi_ice[idx], t_soisno[idx], frac_iceold[idx], snl[idx], dz[idx], z[idx], zi[idx], snw_rds[idx],
                  qflx_snow_h2osfc[idx], frac_sno_eff[idx], frac_sno[idx]);

    std::cout << snow_depth[idx] << " " << h2osno[idx] << " " << int_snow[idx] << " " << snl[idx] << " "
              << qflx_snow_h2osfc[idx] << " " << frac_sno_eff[idx] << " " << frac_sno[idx] << std::endl;

    for (int i = ELM::nlevsno - snl[idx]; i < ELM::nlevsno; i++) {
      std::cout << i << " " << swe_old[idx][i] << " " << h2osoi_liq[idx][i] << " " << h2osoi_ice[idx][i] << " "
                << t_soisno[idx][i] << " " << frac_iceold[idx][i] << " " << dz[idx][i] << " " << z[idx][i] << " "
                << zi[idx][i] << " " << snw_rds[idx][i] << std::endl;
    }

    ELM::FracH2OSfc(Land, micro_sigma[idx], h2osno[idx], h2osfc[idx], h2osoi_liq[idx], frac_sno[idx], frac_sno_eff[idx],
                    frac_h2osfc[idx]);

    std::cout << micro_sigma[idx] << " " << h2osno[idx] << " " << h2osfc[idx] << " " << frac_sno[idx] << " "
              << frac_sno_eff[idx] << " " << frac_h2osfc[idx] << std::endl;
  }


  // surface radiation variables
  // these will all be calculated eventually
  // for now, just manually hardwire in data values 
  auto nrad = create<ArrayI1>("nrad", n_grid_cells);
  assign(nrad, 1); // hardwired for nlevcan == 1; gets initialized in SurfaceAlbedoType.F90 & recalculated in SurfaceAlbedoMod.F90
  auto fsr = create<ArrayD1>("fsr", n_grid_cells); // output
  auto laisun = create<ArrayD1>("laisun", n_grid_cells); // output
  auto laisha = create<ArrayD1>("laisha", n_grid_cells); // output
  auto sabg_soil = create<ArrayD1>("sabg_soil", n_grid_cells); // output
  auto sabg_snow = create<ArrayD1>("sabg_snow", n_grid_cells); // output
  auto sabg = create<ArrayD1>("sabg", n_grid_cells); // output
  auto sabv = create<ArrayD1>("sabv", n_grid_cells); // output
  auto fsa = create<ArrayD1>("fsa", n_grid_cells); // output
  auto tlai_z = create<ArrayD2>("tlai_z", n_grid_cells, ELM::nlevcan);
  assign(tlai_z, 0.303834600280739);
  auto fsun_z = create<ArrayD2>("fsun_z", n_grid_cells, ELM::nlevcan);
  //assign(fsun_z, 0.0);
  assign(fsun_z,  0.554991699703339);
  //auto forc_solad = create<ArrayD2>("forc_solad", n_grid_cells, ELM::numrad);
  //auto forc_solai = create<ArrayD2>("forc_solai", n_grid_cells, ELM::numrad);
  auto fabd_sun_z = create<ArrayD2>("fabd_sun_z", n_grid_cells, ELM::nlevcan);
  //assign(fabd_sun_z, 0.0);
  assign(fabd_sun_z, 0.652531958249706);
  auto fabd_sha_z = create<ArrayD2>("fabd_sha_z", n_grid_cells, ELM::nlevcan);
  //assign(fabd_sha_z, 0.0);
  assign(fabd_sha_z, 0.0626477498888005);
  auto fabi_sun_z = create<ArrayD2>("fabi_sun_z", n_grid_cells, ELM::nlevcan);
  //assign(fabi_sun_z, 0.0);
  assign(fabi_sun_z, 0.477194066186709);
  auto fabi_sha_z = create<ArrayD2>("fabi_sha_z", n_grid_cells, ELM::nlevcan);
  //assign(fabi_sha_z, 0.0);
  assign(fabi_sha_z, 0.341279365789308);
  auto parsun_z = create<ArrayD2>("parsun_z", n_grid_cells, ELM::nlevcan); // output
  auto parsha_z = create<ArrayD2>("parsha_z", n_grid_cells, ELM::nlevcan); // output
  auto laisun_z = create<ArrayD2>("laisun_z", n_grid_cells, ELM::nlevcan); // output
  auto laisha_z = create<ArrayD2>("laisha_z", n_grid_cells, ELM::nlevcan); // output
  auto sabg_lyr = create<ArrayD2>("sabg_lyr", n_grid_cells, ELM::nlevsno + 1);
  auto ftdd = create<ArrayD2>("ftdd", n_grid_cells, ELM::numrad);
  //assign(ftdd, 0.0);
  ftdd[0][0] = 0.266879726153718;
  ftdd[0][1] = 0.266879726153718;
  auto ftid = create<ArrayD2>("ftid", n_grid_cells, ELM::numrad);
  //assign(ftid, 0.0);
  ftid[0][0] = 0.048196184434411;
  ftid[0][1] = 0.228995902545678;
  auto ftii = create<ArrayD2>("ftii", n_grid_cells, ELM::numrad);
  //assign(ftii, 0.0);
  ftii[0][0] = 0.247622223576231;
  ftii[0][1] = 0.420713437119143;
  auto fabd = create<ArrayD2>("fabd", n_grid_cells, ELM::numrad);
  //assign(fabd, 0.0);
  fabd[0][0] = 0.666823073554822;
  fabd[0][1] = 0.374008118438534;
  auto fabi = create<ArrayD2>("fabi", n_grid_cells, ELM::numrad);
  //assign(fabi, 0.0);
  fabi[0][0] = 0.712441211538226;
  fabi[0][1] = 0.391710336687883;
  auto albsod = create<ArrayD2>("albsod", n_grid_cells, ELM::numrad);
  assign(albsod, 0.2);
  auto albsoi = create<ArrayD2>("albsoi", n_grid_cells, ELM::numrad);
  assign(albsoi, 0.2);
  auto albsnd_hst = create<ArrayD2>("albsnd_hst", n_grid_cells, ELM::numrad);
  assign(albsnd_hst, 0.6);
  auto albsni_hst = create<ArrayD2>("albsni_hst", n_grid_cells, ELM::numrad);
  assign(albsni_hst, 0.6);
  auto albgrd = create<ArrayD2>("albgrd", n_grid_cells, ELM::numrad);
  assign(albgrd, 0.0);
  auto albgri = create<ArrayD2>("albgri", n_grid_cells, ELM::numrad);
  assign(albgri, 0.0);
  auto flx_absdv = create<ArrayD2>("flx_absdv", n_grid_cells, ELM::nlevsno + 1);
  assign(flx_absdv, 0.0);
  auto flx_absdn = create<ArrayD2>("flx_absdn", n_grid_cells, ELM::nlevsno + 1);
  assign(flx_absdn, 0.0);
  auto flx_absiv = create<ArrayD2>("flx_absiv", n_grid_cells, ELM::nlevsno + 1);
  assign(flx_absiv, 0.0);
  auto flx_absin = create<ArrayD2>("flx_absin", n_grid_cells, ELM::nlevsno + 1);
  assign(flx_absin, 0.0);
  auto albd = create<ArrayD2>("albd", n_grid_cells, ELM::numrad);
  albd[0][0] = 0.13;
  albd[0][0] = 0.26;
  auto albi = create<ArrayD2>("albi", n_grid_cells, ELM::numrad);
  albi[0][0] = 0.13;
  albi[0][0] = 0.26;



  ELM::CanopySunShadeFractions(
        Land, nrad[idx], elai[idx], tlai_z[idx], fsun_z[idx],
        forc_solad[idx], forc_solai[idx],
        fabd_sun_z[idx], fabd_sha_z[idx],
        fabi_sun_z[idx], fabi_sha_z[idx],
        parsun_z[idx], parsha_z[idx],
        laisun_z[idx], laisha_z[idx], laisun[idx], laisha[idx]);
  std::cout << laisun[idx] << " " << laisha[idx] << std::endl;
    std::cout << "CanopySunShadeFractions: " << parsun_z[0][0] << " " <<parsha_z[0][0] << " " << 
      laisun_z[0][0] << " " << laisha_z[0][0] << std::endl;


   ELM::SurfRadZeroFluxes(Land, sabg_soil[idx], sabg_snow[idx], sabg[idx], sabv[idx], fsa[idx],
                          sabg_lyr[idx]);

   {                          // scope for some local vars
     double trd[ELM::numrad]; // transmitted solar radiation: direct (W/m**2)
     double tri[ELM::numrad]; // transmitted solar radiation: diffuse (W/m**2)
      ELM::SurfRadAbsorbed(Land, snl[idx], ftdd[idx], ftid[idx],
         ftii[idx], forc_solad[idx],
         forc_solai[idx], fabd[idx],
         fabi[idx], albsod[idx],
         albsoi[idx], albsnd_hst[idx],
         albsni_hst[idx], albgrd[idx],
         albgri[idx], sabv[idx], fsa[idx], sabg[idx], sabg_soil[idx], sabg_snow[idx],
         trd, tri);

      std::cout << "SurfRadAbsorbed: " << sabv[idx] << " " << fsa[idx] << " " << sabg[idx] << " " << sabg_soil[idx] << " " << sabg_snow[idx]  << " " <<
      trd[0]  << " " << trd[1]  << " " << tri[0]  << " " << tri[1] << std::endl;

     ELM::SurfRadLayers(Land, snl[idx], sabg[idx], sabg_snow[idx], snow_depth[idx], flx_absdv[idx],
                       flx_absdn[idx], flx_absiv[idx],
                       flx_absin[idx], trd, tri, sabg_lyr[idx]);

     std::cout << "SurfRadLayers, with snl =: " << snl[0] << std::endl;
     for (int i = 0; i < ELM::nlevsno + 1; ++i)
      std::cout << sabg_lyr[0][i] << std::endl;

    ELM::SurfRadReflected(Land, albd[idx], albi[idx],
                          forc_solad[idx], forc_solai[idx],
                          fsr[idx]);
    std::cout << "SurfRadReflected: " << fsr[0] <<  std::endl;
   }

   // variables for CanopyTemperature
   auto t_h2osfc = create<ArrayD1>("t_h2osfc", n_grid_cells);
   assign(t_h2osfc, 255.315686526596);
   auto t_h2osfc_bef = create<ArrayD1>("t_h2osfc_bef", n_grid_cells);
   auto smpmin = create<ArrayD1>("smpmin", n_grid_cells);
   assign(smpmin, -1.0e8); // gets init to this value in SoilStateType.F90
   auto z_0_town = create<ArrayD1>("z_0_town", n_grid_cells);
   auto z_d_town = create<ArrayD1>("z_d_town", n_grid_cells);
   auto soilalpha = create<ArrayD1>("soilalpha", n_grid_cells);
   auto soilalpha_u = create<ArrayD1>("soilalpha_u", n_grid_cells);
   auto soilbeta = create<ArrayD1>("soilbeta", n_grid_cells);
   auto qg_snow = create<ArrayD1>("qg_snow", n_grid_cells);
   auto qg_soil = create<ArrayD1>("qg_soil", n_grid_cells);
   auto qg = create<ArrayD1>("qg", n_grid_cells);
   auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", n_grid_cells);
   auto dqgdT = create<ArrayD1>("dqgdT", n_grid_cells);
   auto htvp = create<ArrayD1>("htvp", n_grid_cells);
   auto emg = create<ArrayD1>("emg", n_grid_cells);
   auto emv = create<ArrayD1>("emv", n_grid_cells);
   auto z0mg = create<ArrayD1>("z0mg", n_grid_cells);
   auto z0hg = create<ArrayD1>("z0hg", n_grid_cells);
   auto z0qg = create<ArrayD1>("z0qg", n_grid_cells);
   auto z0mv = create<ArrayD1>("z0mv", n_grid_cells);
   auto z0hv = create<ArrayD1>("z0hv", n_grid_cells);
   auto z0qv = create<ArrayD1>("z0qv", n_grid_cells);
   //auto beta = create<ArrayD1>("beta", n_grid_cells); - not needed - always assigned 1 - move out of CanTemp and just assign as 1.0 in *Fluxes
   //auto zii = create<ArrayD1>("zii", n_grid_cells);  - not needed - always assigned 1000 - move out of CanTemp and just assign as 1.0 in *Fluxes
   auto thv = create<ArrayD1>("thv", n_grid_cells);
   auto z0m = create<ArrayD1>("z0m", n_grid_cells);
   auto displa = create<ArrayD1>("displa", n_grid_cells);
   //auto cgrnds = create<ArrayD1>("cgrnds", n_grid_cells);
   //auto cgrndl = create<ArrayD1>("cgrndl", n_grid_cells);
   //auto cgrnd = create<ArrayD1>("cgrnd", n_grid_cells);
   auto forc_hgt_u_patch = create<ArrayD1>("forc_hgt_u_patch", n_grid_cells);
   auto forc_hgt_t_patch = create<ArrayD1>("forc_hgt_t_patch", n_grid_cells);
   auto forc_hgt_q_patch = create<ArrayD1>("forc_hgt_q_patch", n_grid_cells);
   auto thm = create<ArrayD1>("thm", n_grid_cells);
   auto eflx_sh_tot = create<ArrayD1>("eflx_sh_tot", n_grid_cells);
   auto eflx_sh_tot_u = create<ArrayD1>("eflx_sh_tot_u", n_grid_cells);
   auto eflx_sh_tot_r = create<ArrayD1>("eflx_sh_tot_r", n_grid_cells);
   auto eflx_lh_tot = create<ArrayD1>("eflx_lh_tot", n_grid_cells);
   auto eflx_lh_tot_u = create<ArrayD1>("eflx_lh_tot_u", n_grid_cells);
   auto eflx_lh_tot_r = create<ArrayD1>("eflx_lh_tot_r", n_grid_cells);
   auto eflx_sh_veg = create<ArrayD1>("eflx_sh_veg", n_grid_cells);
   auto qflx_evap_tot = create<ArrayD1>("qflx_evap_tot", n_grid_cells);
   auto qflx_evap_veg = create<ArrayD1>("qflx_evap_veg", n_grid_cells);
   auto qflx_tran_veg = create<ArrayD1>("qflx_tran_veg", n_grid_cells);
   auto tssbef = create<ArrayD2>("tssbef", n_grid_cells, ELM::nlevgrnd + ELM::nlevsno);
   auto watsat = create<ArrayD2>("watsat", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
   auto sucsat = create<ArrayD2>("sucsat", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
   auto bsw = create<ArrayD2>("bsw", n_grid_cells, ELM::nlevgrnd);       // comes from SoilStateType.F90
   auto watdry = create<ArrayD2>("watdry", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
   auto watopt = create<ArrayD2>("watopt", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
   auto rootfr_road_perv =
       create<ArrayD2>("rootfr_road_perv", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
   auto rootr_road_perv =
       create<ArrayD2>("rootr_road_perv", n_grid_cells, ELM::nlevgrnd); // comes from SoilStateType.F90
   auto watfc = create<ArrayD2>("watfc", n_grid_cells, ELM::nlevgrnd);  // comes from SoilStateType.F90

   ELM::SaveGroundTemp(Land, t_h2osfc[idx], t_soisno[idx], t_h2osfc_bef[idx], tssbef[idx]);

   std::cout << "SaveGroundTemp: " << t_h2osfc_bef[idx] << std::endl;
   for (int i = 0; i < ELM::nlevgrnd + ELM::nlevsno; ++i)
     std::cout << i << " " << tssbef[idx][i] << std::endl;

   ELM::CalculateGroundTemp(Land, snl[idx], frac_sno_eff[idx], frac_h2osfc[idx], t_h2osfc[idx], t_soisno[idx],
                            t_grnd[idx]);
   std::cout << "CalculateGroundTemp: " << t_grnd[idx] << std::endl;

   double qred, hr;
   // need vars from SoilStateType.F90 (InitSoil.hh) to be meaningful for next three calls
   ELM::CalculateSoilAlpha(Land, frac_sno[idx], frac_h2osfc[idx], smpmin[idx], h2osoi_liq[idx], h2osoi_ice[idx],
                           dz[idx], t_soisno[idx], watsat[idx], sucsat[idx], bsw[idx], watdry[idx], watopt[idx],
                           rootfr_road_perv[idx], rootr_road_perv[idx], qred, hr, soilalpha[idx], soilalpha_u[idx]);

   ELM::CalculateSoilBeta(Land, frac_sno[idx], frac_h2osfc[idx], watsat[idx], watfc[idx], h2osoi_liq[idx],
                          h2osoi_ice[idx], dz[idx], soilbeta[idx]);

   ELM::CalculateHumidities(Land, snl[idx], forc_q[idx], forc_pbot[idx], t_h2osfc[idx], t_grnd[idx], frac_sno[idx],
                            frac_sno_eff[idx], frac_h2osfc[idx], qred, hr, t_soisno[idx], qg_snow[idx], qg_soil[idx],
                            qg[idx], qg_h2osfc[idx], dqgdT[idx]);

   return 0;
}
