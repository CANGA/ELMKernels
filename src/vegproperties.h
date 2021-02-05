#pragma once

#include "clm_constants.h"

struct VegProperties {
  double fnr[numpft];      //      fraction of nitrogen in RuBisCO
  double act25[numpft];    //    (umol/mgRubisco/min) Rubisco activity at 25 C
  double kcha[numpft];     //     Activation energy for kc
  double koha[numpft];     //     Activation energy for ko
  double cpha[numpft];     //     Activation energy for cp
  double vcmaxha[numpft];  //  Activation energy for vcmax
  double jmaxha[numpft];   //   Activation energy for jmax
  double tpuha[numpft];    //    Activation energy for tpu
  double lmrha[numpft];    //    Acitivation energy for lmr
  double vcmaxhd[numpft];  //  Deactivation energy for vcmax
  double jmaxhd[numpft];   //   Deactivation energy for jmax
  double tpuhd[numpft];    //    Deactivation energy for tpu
  double lmrhd[numpft];    //    Deacitivation energy for lmr
  double lmrse[numpft];    //    SE for lmr
  double qe[numpft];       //       Quantum efficiency
  double theta_cj[numpft]; // empirical curvature parameter for ac, aj photosynthesis co-limitation
  double bbbopt[numpft];   //   Ball-Berry stomatal conductance intercept
  double mbbopt[numpft];   //   Ball-Berry stomatal conductance slope
  double c3psn[numpft];    //    photosynthetic pathway: 0. = c4, 1. = c3
  double slatop[numpft];   //   specific leaf area at top of canopy, projected area basis [m^2/gC]
  double leafcn[numpft];   //   leaf C:N (gC/gN)
  double flnr[numpft];     //     fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
  double fnitr[numpft];    //    foliage nitrogen limitation factor (-)
};
