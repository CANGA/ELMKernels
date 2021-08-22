
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

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }


template <class Array_t> Array_t createtest(const std::string &name, int D0) { return Array_t(name, D0); }

template <typename T>
bool IsAlmostEqual(const T a, const T b, double rel_tol=1e-7, double abs_tol=1e-20) {
  auto diff = std::abs(a - b);
  auto maxreldiff = std::max(std::abs(a),std::abs(b)) * rel_tol;
  if (diff > maxreldiff || diff > abs_tol) {
    std::cout << "NOT EQUAL  a: " << a << " b: " << b << std::endl;
    return false;
  } else {
    std::cout << "EQUAL" << std::endl;
    return true;
  }
}

int main(int argc, char **argv) {

  // domain decomp 
  int MPI_COMM_WORLD;
  const int n_procs = 1;
  const int myrank = 0;
  auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
  auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp, {1, 1},
                                                     {myrank / proc_decomp[1], myrank % proc_decomp[1]});
  auto n_grid_cells = dd.n_local[0] * dd.n_local[1];

  // date and timestep
  const auto start = ELM::Utils::Date(2014, 7, 15, 1800);
  int dtime = 1800;
  ELM::Utils::Date current(start);

  // put this inside loop when built
  std::cout << start.sec << "," << start.doy << "," << start.year << std::endl;


//how-to
  // get ddatetime info
  auto date = current.date();
  int year = std::get<0>(date);
  int month = std::get<1>(date);
  int day = std::get<2>(date);
  //increment seconds
  //current.increment_seconds(dtime);
  //ELM::IO::read_bottom<double>(data_dir, fname_h1.str(), "SOILICE", dd, 0, h2osoi_ice);
  //ELM::IO::read_bottom<double>(data_dir, fname_h1.str(), "SOILLIQ", dd, 0, h2osoi_liq);
  //assign(snl, 0);

  // test new array name
  std::cout << "starting array test" << std::endl;
  auto frac_veg_nosno = createtest<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto copy = frac_veg_nosno;
  std::cout << "array has been declared" << std::endl;
  std::cout << "plz??! " << frac_veg_nosno.getname() << std::endl;
  std::cout << "copy: " << copy.getname() << std::endl;

  // data files 
  const std::string data_dir("/home/jbeisman/Software/elm_test_input/test_CanHydro");
  const std::string basename_h1("ForATS_AK-BEOG_ICB20TRELMBC.elm.h1.");
  const std::string basename_r("ForATS_AK-BEOG_ICB20TRELMBC.elm.r.");
  std::stringstream fname_h1, fname_r;

  

//from here on inside run loop

  while (current.doy == start.doy) { // only run for 1 day

    fname_h1.str(std::string{});
    fname_h1 << basename_h1 << year << "-" << std::setw(2) << std::setfill('0') << month << "-" << 
    std::setw(2) << std::setfill('0') << day << "-" << std::setw(5) << std::setfill('0') << current.sec <<".nc";

    fname_r.str(std::string{});
    fname_r << basename_r << year << "-" << std::setw(2) << std::setfill('0') << month << "-" << 
    std::setw(2) << std::setfill('0') << day << "-" << std::setw(5) << std::setfill('0') << current.sec <<".nc";

    std::cout << fname_h1.str() << std::endl;
    std::cout << fname_r.str() << std::endl;

{
    // hardwired
  int columni = 0; //column
  int pfti = 12; // pft
  int idx = 0; // grid cell
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  double dewmx = 0.1;
  bool do_capsnow(false);
  double qflx_irrig[1] = {0.0};
  int oldfflag = 0; // hardwired ??
  double n_melt[1] = {20.0}; // hardwired - not in input data!! guess!!


  double qflx_candrip;
  double qflx_through_snow;
  double qflx_through_rain;
  double fracsnow;
  double fracrain;


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
  auto frac_iceold = create<ArrayD2>("frac_iceold", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);

  // the grid -- need to init from input or driver data;
  // don't need values for CanopyHydrology
  auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto z = create<ArrayD2>("z", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
  auto zi = create<ArrayD2>("zi", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd + 1);


    
    // read input vars
    // most come from r file, but a few come from h file
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "RAIN", dd, forc_rain);
    std::cout << forc_rain[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "SNOW", dd, forc_snow);
    std::cout << forc_snow[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "TBOT", dd, forc_t);
    std::cout << forc_t[0] << std::endl;

    ELM::IO::read_distributed_scalar<int>(data_dir, fname_r.str(), "FRAC_VEG_NOSNO_ALB", pfti, frac_veg_nosno);
    std::cout << frac_veg_nosno[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "qflx_snow_melt", columni, qflx_snow_melt); // from r file??
    std::cout << qflx_snow_melt[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "SNOW_DEPTH", columni, snow_depth); // from r file??
    std::cout << snow_depth[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "elai", pfti, elai);
    std::cout << elai[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "esai", pfti, esai);
    std::cout << esai[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "H2OCAN", pfti, h2ocan);
    std::cout << h2ocan[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "T_GRND", columni, t_grnd);
    std::cout << t_grnd[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "H2OSNO", columni, h2osno);
    std::cout << h2osno[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "INT_SNOW", columni, int_snow);
    std::cout << int_snow[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "SNLSNO", columni, snl);
    std::cout << "SNL:: " << snl[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "frac_sno", columni, frac_sno);
    std::cout << frac_sno[0] << std::endl;
    ELM::IO::read_distributed_scalar<double>(data_dir, fname_r.str(), "frac_sno_eff", columni, frac_sno_eff);
    std::cout << frac_sno_eff[0] << std::endl;

    // arrays
    ELM::IO::read_distributed_array<double>(data_dir, fname_r.str(), "H2OSOI_LIQ", h2osoi_liq);
    ELM::IO::read_distributed_array<double>(data_dir, fname_r.str(), "H2OSOI_ICE", h2osoi_ice);


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

     // for (int i = 0; i < 20; i++) {
     //   std::cout << i <<  " " << h2osoi_liq[idx][i] << " " << h2osoi_ice[idx][i] << " "
     //             << t_soisno[idx][i] << std::endl;
     //           }

  ELM::FracH2OSfc(Land, micro_sigma[idx], h2osno[idx], h2osfc[idx], h2osoi_liq[idx], frac_sno[idx], frac_sno_eff[idx],
                    frac_h2osfc[idx]);
  std::cout << h2osfc[idx] << " " << frac_sno[idx] << " " << frac_sno_eff[idx] << " " << frac_h2osfc[idx] << std::endl;


    {
      auto H2OCAN = create<ArrayD1>("H2OCAN", n_grid_cells);
      auto QDRIP = create<ArrayD1>("QDRIP", n_grid_cells);
      auto QSNWCPLIQ = create<ArrayD1>("QSNWCPLIQ", n_grid_cells);
      auto QSNWCPICE = create<ArrayD1>("QSNWCPICE", n_grid_cells);
      auto QFLX_SNOW_GRND = create<ArrayD1>("QFLX_SNOW_GRND", n_grid_cells);
      auto QFLX_RAIN_GRND = create<ArrayD1>("QFLX_RAIN_GRND", n_grid_cells);
      auto FWET = create<ArrayD1>("FWET", n_grid_cells);
      auto FDRY = create<ArrayD1>("FDRY", n_grid_cells);
      auto SNOW_DEPTH = create<ArrayD1>("SNOW_DEPTH", n_grid_cells);
      auto H2OSNO = create<ArrayD1>("H2OSNO", n_grid_cells);
      auto INT_SNOW = create<ArrayD1>("INT_SNOW", n_grid_cells);
      //auto snl = create<ArrayI1>("snl", n_grid_cells);
      auto H2OSFC = create<ArrayD1>("H2OSFC", n_grid_cells);
      auto FH2OSFC = create<ArrayD1>("FH2OSFC", n_grid_cells);
      auto FSNO = create<ArrayD1>("FSNO", n_grid_cells);
      auto FSNO_EFF = create<ArrayD1>("FSNO_EFF", n_grid_cells);


      auto TSOI = create<ArrayD2>("TSOI", n_grid_cells, ELM::nlevgrnd);
      auto SNO_T = create<ArrayD2>("SNO_T", n_grid_cells, ELM::nlevsno);
      auto SWE_OLD = create<ArrayD2>("SWE_OLD", n_grid_cells, ELM::nlevsno);
      auto SOILLIQ = create<ArrayD2>("SOILLIQ", n_grid_cells, ELM::nlevgrnd);
      //SNOWLIQ? 
      auto SOILICE = create<ArrayD2>("SOILICE", n_grid_cells, ELM::nlevgrnd);
      //SNOWICE??
      auto FRAC_ICEOLD = create<ArrayD2>("FRAC_ICEOLD", n_grid_cells, ELM::nlevgrnd);
      //auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
      //auto z = create<ArrayD2>("z", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
      //auto zi = create<ArrayD2>("zi", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd + 1);




      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "H2OCAN", 0, pfti, H2OCAN);

      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "QDRIP", 0, pfti, QDRIP);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "QSNWCPLIQ", 0, pfti, QSNWCPLIQ);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "QSNWCPICE", 0, 0, QSNWCPICE);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "QFLX_SNOW_GRND", 0, pfti, QFLX_SNOW_GRND);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "QFLX_RAIN_GRND", 0, pfti, QFLX_RAIN_GRND);

      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "FWET", 0, pfti, FWET);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "FDRY", 0, pfti, FDRY);

      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "SNOW_DEPTH", 0, 0, SNOW_DEPTH);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "H2OSNO", 0, 0, H2OSNO);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "INT_SNOW", 0, 0, INT_SNOW);

      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "FSNO", 0, 0, FSNO);
      ELM::IO::read_distributed_scalar<double>(data_dir, fname_h1.str(), "FSNO_EFF", 0, 0, FSNO_EFF);


      // need to read FRAC_ICEOLD, SOILLIQ, SOILICE, TSOI, SNO_T  - all are in time, levels, column format
      ELM::IO::read_distributed_array_<double>(data_dir, fname_h1.str(), "FRAC_ICEOLD", FRAC_ICEOLD);
      //ELM::IO::read_distributed_array_<double>(data_dir, fname_h1.str(), "FRAC_ICEOLD", FRAC_ICEOLD);
      //ELM::IO::read_distributed_array_<double>(data_dir, fname_h1.str(), "FRAC_ICEOLD", FRAC_ICEOLD);
      //ELM::IO::read_distributed_array_<double>(data_dir, fname_h1.str(), "FRAC_ICEOLD", FRAC_ICEOLD);
      //ELM::IO::read_distributed_array_<double>(data_dir, fname_h1.str(), "FRAC_ICEOLD", FRAC_ICEOLD);



      IsAlmostEqual(h2ocan[idx], H2OCAN[idx]);
      IsAlmostEqual(qflx_prec_grnd[idx], QDRIP[idx]);
      IsAlmostEqual(qflx_snwcp_liq[idx], QSNWCPLIQ[idx]);
      IsAlmostEqual(qflx_snwcp_ice[idx], QSNWCPICE[idx]);
      IsAlmostEqual(qflx_snow_grnd[idx], QFLX_SNOW_GRND[idx]);
      IsAlmostEqual(qflx_rain_grnd[idx], QFLX_RAIN_GRND[idx]);

      IsAlmostEqual(fwet[idx], FWET[idx]);
      IsAlmostEqual(fdry[idx], FDRY[idx]);
      IsAlmostEqual(snow_depth[idx], SNOW_DEPTH[idx]);
      IsAlmostEqual(h2osno[idx], H2OSNO[idx]);
      IsAlmostEqual(int_snow[idx], INT_SNOW[idx]);

      IsAlmostEqual(frac_sno[idx], FSNO[idx]);
      IsAlmostEqual(frac_sno_eff[idx], FSNO_EFF[idx]);


}








  }








  current.increment_seconds(dtime);

}

  return 0;

}

