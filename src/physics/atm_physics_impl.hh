
#pragma once

namespace ELM::atm_forcing_physics {

template <typename ArrayD1>
ConstitutiveAirProperties<ArrayD1>::
ConstitutiveAirProperties(const ArrayD1 forc_qbot, const ArrayD1 forc_pbot,
                          const ArrayD1 forc_tbot,
                          ArrayD1 forc_rho, ArrayD1 forc_po2, ArrayD1 forc_pco2)
    : forc_qbot_{forc_qbot}, forc_pbot_{forc_pbot},
      forc_tbot_{forc_tbot},
      forc_rho_{forc_rho}, forc_po2_{forc_po2},
      forc_pco2_{forc_pco2}
    {}

// functor to calculate all derived forcing quantities
template <typename ArrayD1>
ACCELERATE
constexpr void ConstitutiveAirProperties<ArrayD1>::
operator()(const int i) const
{
  forc_rho_(i) = derive_forc_rho(forc_pbot_(i), forc_qbot_(i), forc_tbot_(i));
  forc_po2_(i) = derive_forc_po2(forc_pbot_(i));
  forc_pco2_(i) = derive_forc_pco2(forc_pbot_(i));
}

template <typename ArrayD1, typename ArrayD2>
ProcessTBOT<ArrayD1, ArrayD2>::
ProcessTBOT(const int& t_idx, const double& wt1, const double& wt2,
            const ArrayD2 atm_tbot, ArrayD1 forc_tbot, ArrayD1 forc_thbot)
    : t_idx_{t_idx}, wt1_{wt1}, wt2_{wt2},
      atm_tbot_{atm_tbot}, forc_tbot_{forc_tbot},
      forc_thbot_{forc_thbot}
    {}

// functor to calculate atmospheric temperature and potential temperature
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
constexpr void ProcessTBOT<ArrayD1, ArrayD2>::
operator()(const int i) const {
  forc_tbot_(i) = std::min(interp_forcing(wt1_, wt2_, atm_tbot_(t_idx_, i), atm_tbot_(t_idx_ + 1, i)), 323.0);
  forc_thbot_(i) = forc_tbot_(i);
}

template <typename ArrayD1, typename ArrayD2>
ProcessPBOT<ArrayD1, ArrayD2>::
ProcessPBOT(const int& t_idx, const double& wt1, const double& wt2,
            const ArrayD2 atm_pbot, ArrayD1 forc_pbot)
    : t_idx_{t_idx}, wt1_{wt1}, wt2_{wt2},
      atm_pbot_{atm_pbot}, forc_pbot_{forc_pbot}
    {}

// functor to calculate atmospheric pressure
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
constexpr void ProcessPBOT<ArrayD1, ArrayD2>::
operator()(const int i) const {
  forc_pbot_(i) = std::max(interp_forcing(wt1_, wt2_, atm_pbot_(t_idx_, i), atm_pbot_(t_idx_ + 1, i)), 4.0e4);
}

template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
ProcessQBOT<ArrayD1, ArrayD2, ftype>::
ProcessQBOT(const int& t_idx, const double& wt1, const double& wt2,
            const ArrayD2 atm_qbot, const ArrayD1 forc_tbot,
            const ArrayD1 forc_pbot, ArrayD1 forc_qbot)
    : t_idx_{t_idx}, wt1_{wt1}, wt2_{wt2},
      atm_qbot_{atm_qbot}, forc_tbot_{forc_tbot},
      forc_pbot_{forc_pbot}, forc_qbot_{forc_qbot}
    {}

// functor to calculate specific humidity and relative humidity
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
ACCELERATE
constexpr void ProcessQBOT<ArrayD1, ArrayD2, ftype>::
operator()(const int i) const {
  forc_qbot_(i) = std::max(interp_forcing(wt1_, wt2_, atm_qbot_(t_idx_, i), atm_qbot_(t_idx_ + 1, i)), 1.0e-9);
  if constexpr (ftype == AtmForcType::RH) {
    double e = (forc_tbot_(i) > ELMconst::TFRZ()) ? esatw(tdc(forc_tbot_(i))) : esati(tdc(forc_tbot_(i)));
    double qsat = 0.622 * e / (forc_pbot_(i) - 0.378 * e);
    forc_qbot_(i) *= qsat / 100.0;
  }
}

template <typename ArrayD1, typename ArrayD2>
ProcessFLDS<ArrayD1, ArrayD2>::
ProcessFLDS(const int& t_idx, const double& wt1, const double& wt2,
            const ArrayD2 atm_flds, const ArrayD1 forc_pbot, const ArrayD1 forc_qbot,
            const ArrayD1 forc_tbot, ArrayD1 forc_lwrad)
    : t_idx_{t_idx}, wt1_{wt1}, wt2_{wt2},
      atm_flds_{atm_flds}, forc_pbot_{forc_pbot},
      forc_qbot_{forc_qbot}, forc_tbot_{forc_tbot},
      forc_lwrad_{forc_lwrad}
    {}

// functor to calculate downward longwave radiation
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
constexpr void ProcessFLDS<ArrayD1, ArrayD2>::
operator()(const int i) const {
  const double flds = interp_forcing(wt1_, wt2_, atm_flds_(t_idx_, i), atm_flds_(t_idx_ + 1, i));
  if (flds <= 50.0 || flds >= 600.0) {
    const double e = forc_pbot_(i) * forc_qbot_(i) / (0.622 + 0.378 * forc_qbot_(i));
    const double ea = 0.70 + 5.95e-5 * 0.01 * e * exp(1500.0 / forc_tbot_(i));
    forc_lwrad_(i) = ea * ELMconst::STEBOL() * pow(forc_tbot_(i), 4.0);
  } else {
    forc_lwrad_(i) = flds;
  }
}

template <typename ArrayD1, typename ArrayD2>
ProcessFSDS<ArrayD1, ArrayD2>::
ProcessFSDS(const int& t_idx, const ArrayD2 atm_fsds,
            const ArrayD1 coszen, ArrayD2 forc_solai,
            ArrayD2 forc_solad)
    : t_idx_{t_idx}, atm_fsds_{atm_fsds},
      coszen_{coszen}, forc_solai_{forc_solai},
      forc_solad_{forc_solad}
    {}

// functor to calculate solar incident and diffuse radiation in the visible and NIR spectrums
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
constexpr void ProcessFSDS<ArrayD1, ArrayD2>::
operator()(const int i) const {
  // need to impement model for coszen factor
  // ELM uses fac = (cosz > 0.001) ? min(cosz/avg_forc_cosz, 10) : 0.0
  // ATS uses a slope based factor
  // ELM's method could probably be calculated outside the parallel region
  // ATS's method should probably be calculated inside parallel region
  const double swndr = std::max(atm_fsds_(t_idx_, i) * coszen_(i) * 0.5, 0.0);
  const double& swndf = swndr;
  const double& swvdr = swndr; // these vars are only used with a specific forcing data stream
  const double& swvdf = swndr; // maybe implement later? placeholders for now
  const double ratio_rvrf_vis = std::min(
      0.99, std::max(0.17639 + 0.00380 * swvdr - 9.0039e-06 * pow(swvdr, 2.0) + 8.1351e-09 * pow(swvdr, 3.0), 0.01));
  const double ratio_rvrf_nir = std::min(
      0.99, std::max(0.29548 + 0.00504 * swndr - 1.4957e-05 * pow(swndr, 2.0) + 1.4881e-08 * pow(swndr, 3.0), 0.01));
  forc_solad_(i, 0) = ratio_rvrf_vis * swvdr;
  forc_solad_(i, 1) = ratio_rvrf_nir * swndr;
  forc_solai_(i, 0) = (1.0 - ratio_rvrf_vis) * swvdf;
  forc_solai_(i, 1) = (1.0 - ratio_rvrf_nir) * swndf;
}

template <typename ArrayD1, typename ArrayD2>
ProcessPREC<ArrayD1, ArrayD2>::
ProcessPREC(const int& t_idx, const ArrayD2 atm_prec,
            const ArrayD1 forc_tbot, ArrayD1 forc_rain,
            ArrayD1 forc_snow)
    : t_idx_{t_idx}, atm_prec_{atm_prec},
      forc_tbot_{forc_tbot}, forc_rain_{forc_rain},
      forc_snow_{forc_snow}
    {}

// functor to calculate liquid and solid precipitation
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
constexpr void ProcessPREC<ArrayD1, ArrayD2>::
operator()(const int i) const {
  const double frac1 = (forc_tbot_(i) - ELMconst::TFRZ()) * 0.5; // ramp near freezing
  const double frac2 = std::min(1.0, std::max(0.0, frac1));        // bound in [0,1]
  forc_rain_(i) = frac2 * std::max(atm_prec_(t_idx_, i), 0.0);
  forc_snow_(i) = (1.0 - frac2) * std::max(atm_prec_(t_idx_, i), 0.0);
}

template <typename ArrayD1, typename ArrayD2>
ProcessWIND<ArrayD1, ArrayD2>::
ProcessWIND(const int& t_idx, const double& wt1, const double& wt2,
            const ArrayD2 atm_wind, ArrayD1 forc_u, ArrayD1 forc_v)
    : t_idx_{t_idx}, wt1_{wt1}, wt2_{wt2},
      atm_wind_{atm_wind}, forc_u_{forc_u},
      forc_v_{forc_v}
    {}

// functor to calculate wind speed
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
constexpr void ProcessWIND<ArrayD1, ArrayD2>::
operator()(const int i) const {
  forc_u_(i) = interp_forcing(wt1_, wt2_, atm_wind_(t_idx_, i), atm_wind_(t_idx_ + 1, i));
  forc_v_(i) = 0.0;
}

// functor to calculate forcing height
// hardwired at 30m for now
template <typename ArrayD1>
ProcessZBOT<ArrayD1>::
ProcessZBOT(ArrayD1 forc_hgt, ArrayD1 forc_hgt_u,
            ArrayD1 forc_hgt_t, ArrayD1 forc_hgt_q)
    : forc_hgt_{forc_hgt}, forc_hgt_u_{forc_hgt_u},
      forc_hgt_t_{forc_hgt_t}, forc_hgt_q_{forc_hgt_q}
    {}

// hardwired at 30m for now
template <typename ArrayD1>
ACCELERATE
constexpr void ProcessZBOT<ArrayD1>::
operator()(const int i) const {
  forc_hgt_(i) = 30.0;           // hardwired? what about zbot from forcing file?
  forc_hgt_u_(i) = forc_hgt_(i); // observational height of wind [m]
  forc_hgt_t_(i) = forc_hgt_(i); // observational height of temperature [m]
  forc_hgt_q_(i) = forc_hgt_(i); // observational height of humidity [m]
}

// calc forcing given two raw forcing inputs and corresponding weights
ACCELERATE
double interp_forcing(const double& wt1, const double& wt2, const double& forc1, const double& forc2)
{
  return forc1 * wt1 + forc2 * wt2;
}

// convert degrees K to C; bound on interval [-50,50]
ACCELERATE
double tdc(const double& t) { return std::min(50.0, std::max(-50.0, (t - ELMconst::TFRZ()))); }

// calc saturated vapor pressure as function of temp for t > freezing
// Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure.
ACCELERATE
double esatw(const double& t)
{
  static constexpr double a0{6.107799961};
  static constexpr double a1{4.436518521e-01};
  static constexpr double a2{1.428945805e-02};
  static constexpr double a3{2.650648471e-04};
  static constexpr double a4{3.031240396e-06};
  static constexpr double a5{2.034080948e-08};
  static constexpr double a6{6.136820929e-11};
  return 100.0 * (a0 + t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t * a6))))));
}

// calc saturated vapor pressure as function of temp for t <= freezing
// Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure.
ACCELERATE
double esati(const double& t)
{
  static constexpr double b0{6.109177956};
  static constexpr double b1{5.034698970e-01};
  static constexpr double b2{1.886013408e-02};
  static constexpr double b3{4.176223716e-04};
  static constexpr double b4{5.824720280e-06};
  static constexpr double b5{4.838803174e-08};
  static constexpr double b6{1.838826904e-10};
  return 100.0 * (b0 + t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * (b5 + t * b6))))));
}

// vp, rho, pO2, pCO2
// eq 26.10 in CLM tech note
// derive atmospheric vapor pressure from specific humidity and pressure
ACCELERATE
double derive_forc_vp(const double& forc_qbot, const double& forc_pbot)
{
  return forc_qbot * forc_pbot / (0.622 + 0.378 * forc_qbot);
}

// derive atmospheric density from pressure, specific humidity, and temperature
ACCELERATE
double derive_forc_rho(const double& forc_pbot, const double& forc_qbot, const double& forc_tbot)
{
  return (forc_pbot - 0.378 * derive_forc_vp(forc_qbot, forc_pbot)) / (ELMconst::RAIR() * forc_tbot);
}

// derive partial O2 pressure from atmospheric pressure
ACCELERATE
double derive_forc_po2(const double& forc_pbot) { return ELMconst::O2_MOLAR_CONST() * forc_pbot; }

// derive partial CO2 pressure from atmospheric pressure
ACCELERATE
double derive_forc_pco2(const double& forc_pbot) { return ELMconst::CO2_PPMV() * 1.0e-6 * forc_pbot; }

} // namespace ELM::atm_forcing_physics
