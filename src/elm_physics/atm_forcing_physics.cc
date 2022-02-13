
#include "atm_forcing_physics.h"


namespace ns = ELM::atm_forcing_physics;

// calc forcing given two raw forcing inputs and corresponding weights 
double ns::interp_forcing(const double& wt1, const double& wt2, const double& forc1, 
  const double& forc2) { return forc1 * wt1 + forc2 * wt2; }

// convert degrees K to C; bound on interval [-50,50]
double ns::tdc(const double &t) {
  return std::min(50.0, std::max(-50.0, (t - ELMconstants::TFRZ))); }

// calc saturated vapor pressure as function of temp for t > freezing
// Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure. 
double ns::esatw(const double &t) {
  static constexpr double a[7] = {6.107799961, 4.436518521e-01, 1.428945805e-02, 2.650648471e-04,
                     3.031240396e-06, 2.034080948e-08, 6.136820929e-11};
  return 100.0 * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * (a[4] + t * (a[5] + t * a[6])))))); }

// calc saturated vapor pressure as function of temp for t <= freezing
// Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure.
double ns::esati(const double &t) {
  static constexpr double b[7] = {6.109177956, 5.034698970e-01, 1.886013408e-02, 4.176223716e-04,
                     5.824720280e-06, 4.838803174e-08, 1.838826904e-10};
  return 100.0 * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * (b[4] + t * (b[5] + t * b[6])))))); }

// vp, rho, pO2, pCO2
// eq 26.10 in CLM tech note 
// derive atmospheric vapor pressure from specific humidity and pressure
double ns::derive_forc_vp (const double& forc_qbot, const double& forc_pbot)
{ return forc_qbot * forc_pbot / (0.622 + 0.378 * forc_qbot); }

// derive atmospheric density from pressure, vapor pressure, and temperature
double ns::derive_forc_rho (const double& forc_pbot, const double& forc_vp, const double& forc_tbot)
{ return (forc_pbot - 0.378 * forc_vp) / (ELMconstants::RAIR * forc_tbot); }

// derive partial O2 pressure from atmospheric pressure
double ns::derive_forc_po2 (const double& forc_pbot)
{ return ELMconstants::O2_MOLAR_CONST * forc_pbot; }

// derive partial CO2 pressure from atmospheric pressure
double ns::derive_forc_pco2 (const double& forc_pbot)
{ return ELMconstants::CO2_PPMV * 1.0e-6 * forc_pbot; }

