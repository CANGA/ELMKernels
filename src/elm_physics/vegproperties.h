#pragma once

#include "ELMConstants.h"
namespace ELM {
struct VegProperties {
  double fnr[25];      //      fraction of nitrogen in RuBisCO
  double act25[25];    //    (umol/mgRubisco/min) Rubisco activity at 25 C
  double kcha[25];     //     Activation energy for kc
  double koha[25];     //     Activation energy for ko
  double cpha[25];     //     Activation energy for cp
  double vcmaxha[25];  //  Activation energy for vcmax
  double jmaxha[25];   //   Activation energy for jmax
  double tpuha[25];    //    Activation energy for tpu
  double lmrha[25];    //    Acitivation energy for lmr
  double vcmaxhd[25];  //  Deactivation energy for vcmax
  double jmaxhd[25];   //   Deactivation energy for jmax
  double tpuhd[25];    //    Deactivation energy for tpu
  double lmrhd[25];    //    Deacitivation energy for lmr
  double lmrse[25];    //    SE for lmr
  double qe[25];       //       Quantum efficiency
  double theta_cj[25]; // empirical curvature parameter for ac, aj photosynthesis co-limitation
  double bbbopt[25];   //   Ball-Berry stomatal conductance intercept
  double mbbopt[25];   //   Ball-Berry stomatal conductance slope
  double c3psn[25];    //    photosynthetic pathway: 0. = c4, 1. = c3
  double slatop[25];   //   specific leaf area at top of canopy, projected area basis [m^2/gC]
  double leafcn[25];   //   leaf C:N (gC/gN)
  double flnr[25];     //     fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
  double fnitr[25];    //    foliage nitrogen limitation factor (-)
};

} // namespace ELM
