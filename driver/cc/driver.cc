#include "BareGroundFluxes.h"
#include "SatellitePhenology.hh"
#include "array.hh"
#include "clm_constants.h"
#include "landtype.h"
#include "read_input.hh"
#include "utils.hh"
#include "vegproperties.h"
#include <iostream>
#include <string>

#include "ReadAtmInput.hh"

using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
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

  std::cout << start.sec << "," << start.doy << std::endl;
  ELM::Utils::Date current_start(start);

  // double timwt[2];
  ELM::Array<double, 1> timwt(2);
  int dt = 10000;

  // for (int i = 0; i < 220; i++) {
  //  current_start.increment_seconds(84000);
  // std::cout << i << " " << current_start.sec << " " << current_start.doy << " " << current_start.year << std::endl;
  //}
  //
  // ELM::Utils::Date current_test(current_start);
  // current_test.increment_seconds(500);
  // std::cout << current_start.sec << "," << current_test.sec << std::endl;

  const std::string data_dir("/home/jbeisman/Software/elm_test_input");
  const std::string basename_forc("");
  const std::string basename_phen("surfdata_1x1pt_US-Brw_simyr1850_c360x720_c20190221.nc");

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

  if (myrank == n_procs - 1) {
    std::cout << " dimensions (rank " << myrank << "):" << std::endl
              << "   n_time = " << forc_dims[0] << std::endl
              << "   global_n_lat,lon = " << dd.n_global[0] << "," << dd.n_global[1] << std::endl
              << "   nprocs_lat,lon = " << dd.n_procs[0] << "," << dd.n_procs[1] << std::endl
              << "   local_start = " << dd.start[0] << "," << dd.start[1] << std::endl
              << "   local_n_lat,lon = " << dd.n_local[0] << "," << dd.n_local[1] << std::endl;
  }

  // allocate 2 months of LAI, SAI, HTOP, HBOT
  // ELM::Array<double,2> mlai2t(n_grid_cells, 2);
  // ELM::Array<double,2> msai2t(n_grid_cells, 2);
  // ELM::Array<double,2> mhvt2t(n_grid_cells, 2);
  // ELM::Array<double,2> mhvb2t(n_grid_cells, 2);
  auto mlai2t = create<ArrayD2>("mlai2t", n_grid_cells, 2);
  auto msai2t = create<ArrayD2>("msai2t", n_grid_cells, 2);
  auto mhvt2t = create<ArrayD2>("mhvt2t", n_grid_cells, 2);
  auto mhvb2t = create<ArrayD2>("mhvb2t", n_grid_cells, 2);

  auto vtype = create<ArrayI1>("vtype", n_grid_cells);
  assign(vtype, 4);

  // calculate monthly weights; read 2 months of LAI, SAI, HTOP, HBOT if required
  ELM::interpMonthlyVeg(start, dt, n_pfts, data_dir, basename_phen, dd, vtype, timwt, mlai2t, msai2t, mhvt2t, mhvb2t);

  auto tlai = create<ArrayD1>("tlai", n_grid_cells);
  auto tsai = create<ArrayD1>("tsai", n_grid_cells);
  auto elai = create<ArrayD1>("elai", n_grid_cells);
  auto esai = create<ArrayD1>("esai", n_grid_cells);
  auto htop = create<ArrayD1>("htop", n_grid_cells);
  auto hbot = create<ArrayD1>("hbot", n_grid_cells);
  auto frac_veg_nosno_alb = create<ArrayI1>("frac_veg_nosno_alb", n_grid_cells);

  auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells);
  auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);

  ELM::SatellitePhenology(mlai2t[0], msai2t[0], mhvt2t[0], mhvb2t[0], timwt, vtype[0], snow_depth[0], frac_sno[0],

                          tlai[0], tsai[0], htop[0], hbot[0], elai[0], esai[0], frac_veg_nosno_alb[0]);

  for (int mon = 0; mon < 2; ++mon) {
    for (int i = 0; i != n_grid_cells; i++) {
      std::cout << "mlai: " << i << " " << mon << " " << mlai2t(i, mon) << std::endl;
    }

    for (int i = 0; i != n_grid_cells; i++) {
      std::cout << "msai: " << i << " " << mon << " " << msai2t(i, mon) << std::endl;
    }

    for (int i = 0; i != n_grid_cells; i++) {
      std::cout << "mhgtt: " << i << " " << mon << " " << mhvt2t(i, mon) << std::endl;
    }

    for (int i = 0; i != n_grid_cells; i++) {
      std::cout << "mhgtb: " << i << " " << mon << " " << mhvb2t(i, mon) << std::endl;
    }
  }

  std::cout << "::: " << tlai[0] << " " << tsai[0] << " " << htop[0] << " " << hbot[0] << " " << elai[0] << " "
            << esai[0] << std::endl;

  auto forc_zbot = create<ArrayD2>("forc_zbot", ntimes, n_grid_cells);
  auto forc_tbot = create<ArrayD2>("forc_tbot", ntimes, n_grid_cells);
  auto forc_rh = create<ArrayD2>("forc_rh", ntimes, n_grid_cells);
  auto forc_wind = create<ArrayD2>("forc_wind", ntimes, n_grid_cells);
  auto forc_fsds = create<ArrayD2>("forc_fsds", ntimes, n_grid_cells);
  auto forc_flds = create<ArrayD2>("forc_flds", ntimes, n_grid_cells);
  auto forc_psrf = create<ArrayD2>("forc_psrf", ntimes, n_grid_cells);
  auto forc_prec = create<ArrayD2>("forc_prec", ntimes, n_grid_cells);

  ReadAtmData(data_dir, basename_forc, start, dd, n_months, forc_zbot, forc_tbot, forc_rh, forc_wind, forc_fsds,
              forc_flds, forc_psrf, forc_prec);

  for (int i = 0; i < 1488; ++i) {
    std::cout << "::: " << forc_zbot(i, 0) << " " << forc_tbot(i, 0) << " " << forc_rh(i, 0) << " " << forc_wind(i, 0)
              << " " << forc_fsds(i, 0) << " " << forc_flds(i, 0) << " " << forc_psrf(i, 0) << " " << forc_prec(i, 0)
              << std::endl;
  }

  // instantiate data
  ELM::LandType Land;
  ELM::VegProperties Veg;
  double dtime = 0.5;
  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto snl = create<ArrayI1>("snl", n_grid_cells);
  assign(snl, 5);
  auto n_irrig_steps_left = create<ArrayI1>("n_irrig_steps_left", n_grid_cells);
  auto forc_rain = create<ArrayD1>("forc_rain", n_grid_cells);
  auto forc_snow = create<ArrayD1>("forc_snow", n_grid_cells);
  // auto elai = create<ArrayD1>("elai", n_grid_cells);
  // auto esai = create<ArrayD1>("esai", n_grid_cells);
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
  auto forc_t = create<ArrayD1>("forc_t", n_grid_cells);
  auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
  auto qflx_snow_melt = create<ArrayD1>("qflx_snow_melt", n_grid_cells);
  auto n_melt = create<ArrayD1>("n_melt", n_grid_cells);
  auto micro_sigma = create<ArrayD1>("micro_sigma", n_grid_cells);
  auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);
  auto int_snow = create<ArrayD1>("int_snow", n_grid_cells);
  auto qflx_snow_h2osfc = create<ArrayD1>("qflx_snow_h2osfc", n_grid_cells);
  auto h2osfc = create<ArrayD1>("h2osfc", n_grid_cells);
  auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
  auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", n_grid_cells);
  auto swe_old = create<ArrayD2>("swe_old", n_grid_cells, ELM::nlevsno);
  auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto frac_iceold = create<ArrayD2>("frac_iceold", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto z = create<ArrayD2>("z", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto zi = create<ArrayD2>("zi", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd + 1);
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
