#pragma once

#include "read_input.hh"
#include "utils.hh"
namespace ELM {

// time interpolation? - later
// serial I/O function
template <typename Array_t>
void ReadAtmData(const std::string &data_dir, const std::string &basename_atm, const Utils::Date &time,
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

void InitAtmTimestep() {

  // Constants to compute vapor pressure
  parameter(a0 = 6.107799961_r8, a1 = 4.436518521e-01_r8, &a2 = 1.428945805e-02_r8, a3 = 2.650648471e-04_r8,
            &a4 = 3.031240396e-06_r8, a5 = 2.034080948e-08_r8, &a6 = 6.136820929e-11_r8)

      parameter(b0 = 6.109177956_r8, b1 = 5.034698970e-01_r8, &b2 = 1.886013408e-02_r8, b3 = 4.176223716e-04_r8,
                &b4 = 5.824720280e-06_r8, b5 = 4.838803174e-08_r8, &b6 = 1.838826904e-10_r8)

          tdc(t) = min(50._r8, max(-50._r8, (t - SHR_CONST_TKFRZ))) esatw(t) =
              100._r8 * (a0 + t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t * a6)))))) esati(t) =
                  100._r8 * (b0 + t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * (b5 + t * b6))))))

                                forc_t = std::min(atm_tbot, 323.0);
  forc_th = forc_t;
  forc_pbot = std::max(atm_psrf, 4.0e4);
  forc_q = std::max(atm_rh, 1.0e-9);

  if (atm_flds <= 50.0 || atm_flds >= 600.0) {
    e = forc_pbot * forc_q / (0.622 + 0.378 * forc_q);
    ea = 0.70 + 5.95e-5 * 0.01 * e * exp(1500.0 / forc_t) forc_lwrad = ea * sb * pow(forc_t, 4.0);
  } else {
    forc_lwrad = atm_flds;
  }

  // need to write szenith
  // ultimately need wt2, which will be 1 or 0
  wt2 = 1.0;

  swndr = std::max(atm_fsds * wt2 * 0.50, 0.0);
  swndf = swndr;
  swvdr = swndr;
  swvdf = swndr;
  ratio_rvrf = std::min(
      0.99, std::max(0.29548 + 0.00504 * swndr - 1.4957e-05 * pow(swndr, 20) + 1.4881e-08 * pow(swndr, 3.0), 0.01));
  forc_solad[1] = ratio_rvrf * swndr;
  forc_solai[1] = (1.0 - ratio_rvrf) * swndf;
  ratio_rvrf = std::min(
      0.99, std::max(0.17639 + 0.00380 * swvdr - 9.0039e-06 * pow(swvdr, 2.0) + 8.1351e-09 * pow(swvdr, 3.0), 0.01));
  forc_solad[0] = ratio_rvrf * swvdr;
  forc_solai[0] = (1.0 - ratio_rvrf) * swvdf;

  frac = (forc_t - tfrz) * 0.5;              // ramp near freezing
  frac = std::min(1.0, std::max(0.0, frac)); // bound in [0,1]
  // Don't interpolate rainfall data
  forc_rainc = 0.1 * frac * std::max(atm_prec, 0.0_r8);
  forc_rainl = 0.9 * frac * std::max(atm_prec, 0.0_r8);

  forc_snowc = 0.1 * (1.0 - frac) * std::max(atm_prec, 0.0_r8);
  forc_snowl = 0.9 * (1.0 - frac) * std::max(atm_prec, 0.0_r8);

  forc_u = atm_wind;
  forc_v = 0.0;

  forc_hgt = 30.0; // hardwired? what about zbot from forcing file?

  forc_hgt_u = forc_hgt; // observational height of wind [m]
  forc_hgt_t = forc_hgt; // observational height of temperature [m]
  forc_hgt_q = forc_hgt; // observational height of humidity [m]
  forc_vp = forc_q * forc_pbot / (0.622 + 0.378 * forc_q);
  forc_rho = (forc_pbot - 0.378 * forc_vp) / (rair * forc_t);
  forc_po2 = o2_molar_const * forc_pbot;
  forc_wind = forc_u;
  forc_solar = forc_solad[0] + forc_solai[0] + forc_solad[1] + forc_solai[1];
  forc_rain = forc_rainc + forc_rainl;
  forc_snow = forc_snowc + forc_snowl;
  if (forc_t > tfrz) {
    e = esatw(tdc(forc_t))
  } else {
    e = esati(tdc(forc_t))
  }
  qsat = 0.622_r8 * e / (forc_pbot - 0.378_r8 * e) atm2lnd_vars % forc_rh_grc(g) = 100.0_r8 * (forc_q / qsat)
}

double szenith() {}

} // namespace ELM