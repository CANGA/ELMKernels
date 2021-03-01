
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"
#include "vegproperties.h"
#include <iostream>
#include <string>

#include "BareGroundFluxes.h"
#include "InitSnowLayers.hh"
#include "InitTimestep.hh"
#include "ReadAtmosphere.hh"
#include "ReadPFTConstants.hh"
#include "ReadTestData.hh"
#include "SatellitePhenology.hh"
#include "clm_constants.h"
#include "landtype.h"

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

  // the grid
  auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto z = create<ArrayD2>("z", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto zi = create<ArrayD2>("zi", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd + 1);

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

  auto timwt = create<ArrayD1>("timwt", 2);            // time weights for lai/sai/hgt
  auto vtype = create<ArrayI1>("vtype", n_grid_cells); // pft type
  assign(vtype, 4);                                    // dummy variable for now

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

  auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells); // NEED VALUES! -- use init nc file for now
  auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);     // NEED VALUES! -- use init nc file for now

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
    //            << forc_th[idx] << " " << forc_rho[idx] << " " << forc_solad[idx][0] << " " << forc_solad[idx][1] << "
    //            "
    //            << forc_solai[idx][0] << " " << forc_solai[idx][1] << " " << forc_rain[idx] << " " << forc_snow[idx]
    //            << std::endl;
    // std::cout << "forc_rho: " << forc_rho[idx] << " " << forc_pbot[idx] << " " << forc_q[idx] << " " << forc_rh[idx]
    // << std::endl;
  }

  const bool lakpoi(false);
  bool do_capsnow(false);
  auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);         // need value
  auto veg_active = create<ArrayB1>("veg_active", n_grid_cells); // need value
  auto snl = create<ArrayI1>("snl", 17);
  auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd); // need value
  auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd); // need value

  auto h2osno_old = create<ArrayD1>("h2osno_old", n_grid_cells);
  auto eflx_bot = create<ArrayD1>("eflx_bot", n_grid_cells);
  auto qflx_glcice = create<ArrayD1>("qflx_glcice", n_grid_cells);

  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto frac_iceold = create<ArrayD2>("frac_iceold", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);

//  ELM::InitSnowLayers(snow_depth[idx], lakpoi, snl[idx], dz[idx], z[idx], zi[idx]);
//  ELM::InitTimestep(lakpoi, h2osno[idx], veg_active[idx], snl[idx], h2osoi_ice[idx], h2osoi_liq[idx],
//                    frac_veg_nosno_alb[idx], h2osno_old[idx], do_capsnow, eflx_bot[idx], qflx_glcice[idx],
//                    frac_veg_nosno[idx], frac_iceold[idx]);




  auto topo_slope = create<ArrayD1>("topo_slope", n_grid_cells);
  auto topo_std = create<ArrayD1>("topo_std", n_grid_cells);
  ELM::ReadLandData(data_dir, basename_phen, dd, topo_slope, topo_std);
  std::cout << topo_slope[0] << " " << topo_std[0] << std::endl;
  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  ELM::ReadTestInitData(data_dir, basename_init, dd, t_soisno, snl);

  for (int i = 0; i < 17; i++) {
  std::cout << snl[i] << std::endl;
}


  // instantiate data
  ELM::LandType Land;
  ELM::VegProperties Veg;

  auto n_irrig_steps_left = create<ArrayI1>("n_irrig_steps_left", n_grid_cells);
  auto h2ocan = create<ArrayD1>("h2ocan", n_grid_cells);
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
  auto n_melt = create<ArrayD1>("n_melt", n_grid_cells);
  auto micro_sigma = create<ArrayD1>("micro_sigma", n_grid_cells);

  auto int_snow = create<ArrayD1>("int_snow", n_grid_cells);
  auto qflx_snow_h2osfc = create<ArrayD1>("qflx_snow_h2osfc", n_grid_cells);
  auto h2osfc = create<ArrayD1>("h2osfc", n_grid_cells);
  auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
  auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", n_grid_cells);
  auto swe_old = create<ArrayD2>("swe_old", n_grid_cells, ELM::nlevsno);


  auto snw_rds = create<ArrayD2>("snw_rds", n_grid_cells, ELM::nlevsno);

  //  // iterate in time
  //  for (int time = 0; time != ntimes; ++time) {
  //    for (int c = 0; c != n_grid_cells; ++c) {
  //      bg_fluxes[c].InitializeFlux(land, frac_vec_nosno_in[c], forc_u[c], forc_v[c], forc_q_in[c], forc_th_in[c],
  //                                  forc_hgt_u_patch_in[c], thm_in[c], thv_in[c], t_grnd[c], qg[c], z0mg_in[c],
  //                                  dlrad[c], ulrad[c]);
  //
  //      bg_fluxes[c].StabilityIteration(land, forc_hgt_t_patch[c], forc_hgt_q_patch[c], z0mg[c], zii[c], beta[c],
  //      z0hg[c],
  //                                      z0qg[c]);
  //
  //      bg_fluxes[c].ComputeFlux(land, snl[c], forc_rho[c], soilbeta[c], dqgdT[c], htvp[c], t_h2osfc[c], qg_snow[c],
  //                               qg_soil[c], qg_h2osfc[c], t_soisno[c], forc_pbot[c], cgrnds[c], cgrndl[c], cgrnd[c],
  //                               eflx_sh_grnd[c], eflx_sh_tot[c], eflx_sh_snow[c], eflx_sh_soil[c], eflx_sh_h2osfc[c],
  //                               qflx_evap_soi[c], qflx_evap_tot[c], qflx_ev_snow[c], qflx_ev_soil[c],
  //                               qflx_ev_h2osfc[c], t_ref2m[c], t_ref2m_r[c], q_ref2m[c], rh_ref2m[c], rh_ref2m_r[c]);
  //    }
  //    std::cout << "running, t= " << time << std::endl;
  //  }

  return 0;
}
