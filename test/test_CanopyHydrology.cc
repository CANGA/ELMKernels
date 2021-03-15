
#include "ELMConstants.h"
#include "LandType.h"
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"

#include <iostream>
#include <string>

#include "CanopyHydrology.h"

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
  auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp, {1, 1},
                                                     {myrank / proc_decomp[1], myrank % proc_decomp[1]});
  auto n_grid_cells = dd.n_local[0] * dd.n_local[1];

  const auto start = ELM::Utils::Date(2014, 7, 1, 1800);
  const int n_months = 1;
  const int ticks_per_day = 48; // 0.5-hourly data
  const int write_interval = 1;
  int dtime = 1800;

  std::cout << start.sec << "," << start.doy << "," << start.year << std::endl;
  ELM::Utils::Date current(start);

  // input data
  const std::string data_dir("/home/jbeisman/Software/elm_test_input/test_CanHydro");
  const std::string basename_h1("ForATS_AK-BEOG_ICB20TRELMBC.elm.h1.");
  const std::string basename_r("ForATS_AK-BEOG_ICB20TRELMBC.elm.r.");

  auto date = current.date();
  int year = std::get<0>(date);
  int month = std::get<1>(date);
  int day = std::get<2>(date);


  std::stringstream fname_h1_start;
  std::stringstream fname_h1_stop;
  std::stringstream fname_r;

  fname_h1_start << basename_h1 << year << "-" << std::setw(2) << std::setfill('0') << month << "-" << 
  std::setw(2) << std::setfill('0') << day << "-" << std::setw(5) << std::setfill('0') << current.sec <<".nc";

  fname_r << basename_r << year << "-" << std::setw(2) << std::setfill('0') << month << "-" << 
  std::setw(2) << std::setfill('0') << day << "-" << std::setw(5) << std::setfill('0') << current.sec <<".nc";

  current.increment_seconds(dtime);
  fname_h1_stop << basename_h1 << year << "-" << std::setw(2) << std::setfill('0') << month << "-" << 
  std::setw(2) << std::setfill('0') << day << "-" << std::setw(5) << std::setfill('0') << current.sec <<".nc";


  std::cout << fname_h1_start.str() << std::endl;
  std::cout << fname_r.str() << std::endl;
  std::cout << fname_h1_stop.str() << std::endl;

  int columni = 0;
  int pfti = 12;
  int idx = 0;



  double dewmx = 0.1;
  bool do_capsnow(false);
  double qflx_irrig[1] = {0.0};
  int oldfflag = 0; // hardwired ??
  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto forc_rain = create<ArrayD1>("forc_rain", n_grid_cells);
  auto forc_snow = create<ArrayD1>("forc_snow", n_grid_cells);
  auto elai = create<ArrayD1>("elai", n_grid_cells);
  auto esai = create<ArrayD1>("esai", n_grid_cells);
  auto h2ocan = create<ArrayD1>("h2ocan", n_grid_cells);
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
  //auto n_melt = create<ArrayD1>("n_melt", n_grid_cells);
  double n_melt[1] = {20.0}; // hardwired - not in input data!! guess!!
  auto micro_sigma = create<ArrayD1>("micro_sigma", n_grid_cells);
  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells); // NEED VALUES! -- use init nc file for now
  auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);     // NEED VALUES! -- use init nc file for now
  auto qflx_snow_h2osfc = create<ArrayD1>("qflx_snow_h2osfc", n_grid_cells);
  auto h2osfc = create<ArrayD1>("h2osfc", n_grid_cells);
  auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
  auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", n_grid_cells);
  auto swe_old = create<ArrayD2>("swe_old", n_grid_cells, ELM::nlevsno);
  auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);                                       // need value
  auto snw_rds = create<ArrayD2>("snw_rds", n_grid_cells, ELM::nlevsno);
  auto int_snow = create<ArrayD1>("int_snow", n_grid_cells); // need value
  auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd); // need value
  auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd); // need value
  auto snl = create<ArrayI1>("snl", n_grid_cells);
  assign(snl, 0);
  auto frac_iceold = create<ArrayD2>("frac_iceold", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);

  // the grid -- need to init from input or driver data; so far, only above ground snow layers get vals
  auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto z = create<ArrayD2>("z", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto zi = create<ArrayD2>("zi", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd + 1);


  double qflx_candrip;
  double qflx_through_snow;
  double qflx_through_rain;
  double fracsnow;
  double fracrain;



  ELM::IO::read_distributed_scalar<int>(data_dir, fname_r.str(), "FRAC_VEG_NOSNO_ALB", pfti, frac_veg_nosno);
  std::cout << frac_veg_nosno[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "qflx_snow_melt", columni, qflx_snow_melt); // from r file??
  std::cout << qflx_snow_melt[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "SNOW_DEPTH", columni, snow_depth); // from r file??
  std::cout << snow_depth[0] << std::endl;

  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "RAIN", dd, forc_rain);
  std::cout << forc_rain[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "SNOW", dd, forc_snow);
  std::cout << forc_snow[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "TBOT", dd, forc_t);
  std::cout << forc_t[0] << std::endl;

  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "ELAI", 0, pfti, elai);
  std::cout << elai[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "ESAI", 0, pfti, esai);
  std::cout << esai[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "H2OCAN", 0, pfti, h2ocan);
  std::cout << h2ocan[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "TG", 0, columni, t_grnd);
  std::cout << t_grnd[0] << std::endl;
  ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1_start.str(), "H2OSNO", 0, columni, h2osno);
  std::cout << h2osno[0] << std::endl;

  //ELM::IO::read_bottom<double>(data_dir, fname_h1_start.str(), "SOILICE", dd, 0, h2osoi_ice);
  //ELM::IO::read_bottom<double>(data_dir, fname_h1_start.str(), "SOILLIQ", dd, 0, h2osoi_liq);

  ELM::IO::read_distributed_array<double>(data_dir, fname_r.str(), "H2OSOI_LIQ", h2osoi_liq);
  ELM::IO::read_distributed_array<double>(data_dir, fname_r.str(), "H2OSOI_ICE", h2osoi_ice);

  for (int i = 0; i < 20; i ++) {
  std::cout << "i: " << i << " " << h2osoi_ice[0][i] << " " << h2osoi_liq[0][i] << std::endl;
  }


  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;

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

    for (int i = 0; i < 20; i++) {
      std::cout << i <<  " " << h2osoi_liq[idx][i] << " " << h2osoi_ice[idx][i] << " "
                << t_soisno[idx][i] << std::endl;
    }

  ELM::FracH2OSfc(Land, micro_sigma[idx], h2osno[idx], h2osfc[idx], h2osoi_liq[idx], frac_sno[idx], frac_sno_eff[idx],
                    frac_h2osfc[idx]);

  return 0;

}

