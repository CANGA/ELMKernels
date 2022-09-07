
#pragma once

#include "array.hh"
#include "elm_constants.h"

#include "compile_options.hh"

#include <cmath>

namespace ELM::atm_forcing_physics {

// calc forcing given two raw forcing inputs and corresponding weights
ACCELERATE
double interp_forcing(const double& wt1, const double& wt2, const double& forc1, const double& forc2);

// convert degrees K to C; bound on interval [-50,50]
ACCELERATE
double tdc(const double& t);

// calc saturated vapor pressure as function of temp for t > freezing
// Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure.
ACCELERATE
double esatw(const double& t);

// calc saturated vapor pressure as function of temp for t <= freezing
// Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure.
ACCELERATE
double esati(const double& t);

// rho, pO2, pCO2
// eq 26.10 in CLM tech note
// derive atmospheric vapor pressure from specific humidity and pressure
ACCELERATE
double derive_forc_vp(const double& forc_qbot, const double& forc_pbot);

// derive atmospheric density from pressure, vapor pressure, and temperature
ACCELERATE
double derive_forc_rho(const double& forc_pbot, const double& forc_vp, const double& forc_tbot);

// derive partial O2 pressure from atmospheric pressure
ACCELERATE
double derive_forc_po2(const double& forc_pbot);
// derive partial CO2 pressure from atmospheric pressure
ACCELERATE
double derive_forc_pco2(const double& forc_pbot);

// functor to calculate all derived forcing quantities
template <typename ArrayD1>
struct ConstitutiveAirProperties {
  ConstitutiveAirProperties(const ArrayD1 forc_qbot,
                            const ArrayD1 forc_pbot,
                            const ArrayD1 forc_tbot,
                            ArrayD1 forc_vp,
                            ArrayD1 forc_rho,
                            ArrayD1 forc_po2,
                            ArrayD1 forc_pco2);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  ArrayD1 forc_qbot_, forc_pbot_, forc_tbot_;
  ArrayD1 forc_vp_, forc_rho_, forc_po2_, forc_pco2_;
};

// functor to calculate atmospheric temperature and potential temperature
template <typename ArrayD1, typename ArrayD2>
struct ProcessTBOT {
  ProcessTBOT(const int& t_idx,
              const double& wt1,
              const double& wt2,
              const ArrayD2 atm_tbot,
              ArrayD1 forc_tbot,
              ArrayD1 forc_thbot);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_tbot_;
  ArrayD1 forc_tbot_, forc_thbot_;
};

// functor to calculate atmospheric pressure
template <typename ArrayD1, typename ArrayD2>
struct ProcessPBOT {
  ProcessPBOT(const int& t_idx,
              const double& wt1,
              const double& wt2,
              const ArrayD2 atm_pbot,
              ArrayD1 forc_pbot);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_pbot_;
  ArrayD1 forc_pbot_;
};

// functor to calculate specific humidity and relative humidty
template <typename ArrayD1, typename ArrayD2, AtmForcType type>
struct ProcessQBOT {
  ProcessQBOT(const int& t_idx,
              const double& wt1,
              const double& wt2,
              const ArrayD2 atm_qbot,
              const ArrayD1 forc_tbot,
              const ArrayD1 forc_pbot,
              ArrayD1 forc_qbot,
              ArrayD1 forc_rh);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_qbot_;
  ArrayD1 forc_tbot_, forc_pbot_, forc_qbot_, forc_rh_;
};

// functor to calculate downward longwave radiation
template <typename ArrayD1, typename ArrayD2>
struct ProcessFLDS {
  ProcessFLDS(const int& t_idx,
              const double& wt1,
              const double& wt2,
              const ArrayD2 atm_flds,
              const ArrayD1 forc_pbot,
              const ArrayD1 forc_qbot,
              const ArrayD1 forc_tbot,
              ArrayD1 forc_lwrad);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_flds_;
  ArrayD1 forc_pbot_, forc_qbot_, forc_tbot_, forc_lwrad_;
};

// functor to calculate solar incident and diffuse
// radiation in the visible and NIR spectrums
template <typename ArrayD1, typename ArrayD2>
struct ProcessFSDS {
  ProcessFSDS(const int& t_idx,
              const ArrayD2 atm_fsds,
              const ArrayD1 coszen,
              ArrayD2 forc_solai,
              ArrayD2 forc_solad);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  ArrayD2 atm_fsds_;
  ArrayD1 coszen_;
  ArrayD2 forc_solai_, forc_solad_;
};

// functor to calculate liquid and solid precipitation
template <typename ArrayD1, typename ArrayD2>
struct ProcessPREC {
  ProcessPREC(const int& t_idx,
              const ArrayD2 atm_prec,
              const ArrayD1 forc_tbot,
              ArrayD1 forc_rain,
              ArrayD1 forc_snow);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  ArrayD2 atm_prec_;
  ArrayD1 forc_tbot_, forc_rain_, forc_snow_;
};

// functor to calculate wind speed
template <typename ArrayD1, typename ArrayD2>
struct ProcessWIND {
  ProcessWIND(const int& t_idx,
              const double& wt1,
              const double& wt2,
              const ArrayD2 atm_wind,
              ArrayD1 forc_u,
              ArrayD1 forc_v);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_wind_;
  ArrayD1 forc_u_, forc_v_;
};

// functor to calculate forcing height
// hardwired at 30m for now
template <typename ArrayD1>
struct ProcessZBOT {
  ProcessZBOT(ArrayD1 forc_hgt,
              ArrayD1 forc_hgt_u,
              ArrayD1 forc_hgt_t,
              ArrayD1 forc_hgt_q);

  ACCELERATE
  constexpr void operator()(const int i) const;

private:
  ArrayD1 forc_hgt_, forc_hgt_u_, forc_hgt_t_, forc_hgt_q_;
};

} // namespace ELM::atm_forcing_physics

#include "atm_physics_impl.hh"
