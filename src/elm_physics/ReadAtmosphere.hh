#pragma once

#include "ELMConstants.h"
#include "cmath"
#include "read_input.hh"
#include "utils.hh"

namespace ELM {

// serial I/O function
template <typename Array_t>
void ReadAtmForcing(const std::string &data_dir, const std::string &basename_atm, const Utils::Date &time,
                    const Utils::DomainDecomposition<2> &dd, const int n_months, Array_t &atm_zbot, Array_t &atm_tbot,
                    Array_t &atm_rh, Array_t &atm_wind, Array_t &atm_fsds, Array_t &atm_flds, Array_t &atm_psrf,
                    Array_t &atm_prec) {

  IO::read_and_reshape_forcing(data_dir, basename_atm, "ZBOT", time, n_months, dd, atm_zbot);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "TBOT", time, n_months, dd, atm_tbot);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "RH", time, n_months, dd, atm_rh);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "WIND", time, n_months, dd, atm_wind);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "FSDS", time, n_months, dd, atm_fsds);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "FLDS", time, n_months, dd, atm_flds);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "PSRF", time, n_months, dd, atm_psrf);
  IO::read_and_reshape_forcing(data_dir, basename_atm, "PRECTmms", time, n_months, dd, atm_prec);
}

double tdc(const double &t) {
  double ret = std::min(50.0, std::max(-50.0, (t - tfrz)));
  return ret;
}

double esatw(const double &t) {
  const double a[7] = {6.107799961,     4.436518521e-01, 1.428945805e-02, 2.650648471e-04,
                       3.031240396e-06, 2.034080948e-08, 6.136820929e-11};
  double ret = 100.0 * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * (a[4] + t * (a[5] + t * a[6]))))));
  return ret;
}

double esati(const double &t) {
  const double b[7] = {6.109177956,     5.034698970e-01, 1.886013408e-02, 4.176223716e-04,
                       5.824720280e-06, 4.838803174e-08, 1.838826904e-10};
  double ret = 100.0 * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * (b[4] + t * (b[5] + t * b[6]))))));
  return ret;
}

/*

skipping time interpolation for now, will need to include later
currently assumed that model timesteps will be the same as the forcing data, with timestep centered around forcing time

need to figure out what scale_factors and add_offsets are and where they come from, haven't been able to find them in
any forcing file

need to write szenith function to calc solar zenith angle for radiation processing
forc_hgt_grc is hardwired as 30m in lnd_import_export.F90 - should it be here? what about ZBOT from forcing file (2m)?

*/
template <typename ArrayD1>
void GetAtmTimestep(const double &atm_tbot, const double &atm_psrf, const double &atm_rh, const double &atm_flds,
                    const double &atm_fsds, const double &atm_prec, const double &atm_wind, double &forc_t,
                    double &forc_th, double &forc_pbot, double &forc_q, double &forc_lwrad, double &forc_rain,
                    double &forc_snow, double &forc_u, double &forc_v, double &forc_rh, double &forc_rho,
                    double &forc_po2, double &forc_pco2, double &forc_hgt_u, double &forc_hgt_t, double &forc_hgt_q,
                    ArrayD1 forc_solad, ArrayD1 forc_solai) {

  double e;

  // forc_t, forc_th, forc_pbot, forc_q
  forc_t = std::min(atm_tbot, 323.0);
  forc_th = forc_t;
  forc_pbot = std::max(atm_psrf, 4.0e4);

  // vapor pressure
  // if given "RH" in forcing file
  forc_rh = std::max(atm_rh, 1.0e-9);
  if (forc_t > tfrz) {
    e = esatw(tdc(forc_t));
  } else {
    e = esati(tdc(forc_t));
  }
  double qsat = 0.622 * e / (forc_pbot - 0.378 * e);
  forc_q = qsat * forc_rh / 100.0;

  // if given "QBOT" in forcing file
  // forc_q = std::max(atm_rh, 1.0e-9);
  // if (forc_t > tfrz) {
  //  e = esatw(tdc(forc_t));
  //} else {
  //  e = esati(tdc(forc_t));
  //}
  // double qsat = 0.622 * e / (forc_pbot - 0.378 * e);
  // forc_rh = 100.0 * (forc_q / qsat);

  // forc_lwrad
  if (atm_flds <= 50.0 || atm_flds >= 600.0) {
    e = forc_pbot * forc_q / (0.622 + 0.378 * forc_q);
    double ea = 0.70 + 5.95e-5 * 0.01 * e * exp(1500.0 / forc_t);
    forc_lwrad = ea * sb * pow(forc_t, 4.0);
  } else {
    forc_lwrad = atm_flds;
  }

  // solar radiation - forc_solai[2], forc_solad[2] -- harwired for numrad==2
  // need to write szenith & calculate thiscosz, avgcosz (if time interpolating)
  // ultimately need wt2, which will be 1 or 0 (with current model dt == forcing dt assumption)
  // if we time interpolate wt2 can be any real number between 0.0 & 1.0
  double wt2 = 1.0;
  double swndr = std::max(atm_fsds * wt2 * 0.50, 0.0);
  double swndf = swndr;
  double swvdr = swndr;
  double swvdf = swndr;
  double ratio_rvrf = std::min(
      0.99, std::max(0.29548 + 0.00504 * swndr - 1.4957e-05 * pow(swndr, 20) + 1.4881e-08 * pow(swndr, 3.0), 0.01));
  forc_solad[1] = ratio_rvrf * swndr;
  forc_solai[1] = (1.0 - ratio_rvrf) * swndf;
  ratio_rvrf = std::min(
      0.99, std::max(0.17639 + 0.00380 * swvdr - 9.0039e-06 * pow(swvdr, 2.0) + 8.1351e-09 * pow(swvdr, 3.0), 0.01));
  forc_solad[0] = ratio_rvrf * swvdr;
  forc_solai[0] = (1.0 - ratio_rvrf) * swvdf;
  // forc_solar = forc_solad[0] + forc_solai[0] + forc_solad[1] + forc_solai[1];

  // partition precip into rain & snow
  double frac = (forc_t - tfrz) * 0.5;       // ramp near freezing
  frac = std::min(1.0, std::max(0.0, frac)); // bound in [0,1]
  // Don't interpolate rainfall data
  // do we care about convective vs large-scale precip? no? if not, we can dispense with the 10%, 90% split
  // forc_rainc = 0.1 * frac * std::max(atm_prec, 0.0);
  // forc_rainl = 0.9 * frac * std::max(atm_prec, 0.0);
  // forc_snowc = 0.1 * (1.0 - frac) * std::max(atm_prec, 0.0);
  // forc_snowl = 0.9 * (1.0 - frac) * std::max(atm_prec, 0.0);
  forc_rain = frac * std::max(atm_prec, 0.0);
  forc_snow = (1.0 - frac) * std::max(atm_prec, 0.0);

  // wind
  forc_u = atm_wind;
  forc_v = 0.0;

  // forcing height
  double forc_hgt = 30.0; // hardwired? what about zbot from forcing file?
  forc_hgt_u = forc_hgt;  // observational height of wind [m]
  forc_hgt_t = forc_hgt;  // observational height of temperature [m]
  forc_hgt_q = forc_hgt;  // observational height of humidity [m]

  // rho, pO2, pCO2
  double forc_vp = forc_q * forc_pbot / (0.622 + 0.378 * forc_q);
  forc_rho = (forc_pbot - 0.378 * forc_vp) / (rair * forc_t);
  forc_po2 = o2_molar_const * forc_pbot;
  forc_pco2 = co2_ppmv * 1.0e-6 * forc_pbot;
  // forc_rain = forc_rainc + forc_rainl;
  // forc_snow = forc_snowc + forc_snowl;
}

// double szenith() {} later?

} // namespace ELM