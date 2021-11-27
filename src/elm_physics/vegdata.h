// Vegetation parameters data structure
// read arrays from NetCDF file for SoA approach
// read single pft data into struct for AoS approach
#pragma once

#include "elm_constants.h"
#include "array.hh"
#include "read_input.hh"
#include "utils.hh"
#include <assert.h>


namespace ELM {

// struct for single pft vegetation parameters used in photosynthesis/canopy_flux
struct PSNVegData {
  double fnr, act25, kcha, koha, cpha, vcmaxha, jmaxha, tpuha, lmrha;
  double vcmaxhd, jmaxhd, tpuhd, lmrhd, lmrse, qe, theta_cj, bbbopt, mbbopt;
  double c3psn, slatop, leafcn, flnr, fnitr, dleaf,smpso,smpsc,tc_stress;
};

// struct for single pft vegetation parameters used in surface_albedo
struct AlbedoVegData {
  double rhol[ELM::numrad], rhos[ELM::numrad];
  double taul[ELM::numrad], taus[ELM::numrad];
  double xl;
};

// struct that stores array objects containing time-invariant vegetation data
template <typename ArrayD1, typename ArrayD2>
struct VegData {

  ArrayD1 fnr;        //  fraction of nitrogen in RuBisCO
  ArrayD1 act25;      //  Rubisco activity at 25 C (umol/mgRubisco/min)
  ArrayD1 kcha;       //  Activation energy for kc
  ArrayD1 koha;       //  Activation energy for ko
  ArrayD1 cpha;       //  Activation energy for cp
  ArrayD1 vcmaxha;    //  Activation energy for vcmax
  ArrayD1 jmaxha;     //  Activation energy for jmax
  ArrayD1 tpuha;      //  Activation energy for tpu
  ArrayD1 lmrha;      //  Acitivation energy for lmr
  ArrayD1 vcmaxhd;    //  Deactivation energy for vcmax
  ArrayD1 jmaxhd;     //  Deactivation energy for jmax
  ArrayD1 tpuhd;      //  Deactivation energy for tpu
  ArrayD1 lmrhd;      //  Deacitivation energy for lmr
  ArrayD1 lmrse;      //  SE for lmr
  ArrayD1 qe;         //  Quantum efficiency
  ArrayD1 theta_cj;   //  empirical curvature parameter for ac, aj photosynthesis co-limitation
  ArrayD1 bbbopt;     //  Ball-Berry stomatal conductance intercept
  ArrayD1 mbbopt;     //  Ball-Berry stomatal conductance slope
  ArrayD1 c3psn;      //  photosynthetic pathway: 0. = c4, 1. = c3
  ArrayD1 slatop;     //  specific leaf area at top of canopy, projected area basis [m^2/gC]
  ArrayD1 leafcn;     //  leaf C:N (gC/gN)
  ArrayD1 flnr;       //  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
  ArrayD1 fnitr;      //  foliage nitrogen limitation factor (-)
  ArrayD1 dleaf;      //  heat transfer coefficient from leaves [-]
  ArrayD1 smpso;      //  soil water potential at full stomatal opening (mm)
  ArrayD1 smpsc;      //  soil water potential at full stomatal closure (mm)
  ArrayD1 tc_stress;  //  critical soil temperature for soil water stress (C)
  ArrayD1 z0mr;       //  ratio of momentum roughness length to canopy top height (-)
  ArrayD1 displar;    //  ratio of displacement height to canopy top height (-)
  ArrayD1 xl;         //  ecophys const - leaf/stem orientation index
  ArrayD1 roota_par;  //  rooting distribution parameter [1/m]
  ArrayD1 rootb_par;  //  rooting distribution parameter [1/m]
  ArrayD1 rholvis;    //  visible leaf reflectance
  ArrayD1 rholnir;    //  nir leaf reflectance
  ArrayD1 rhosvis;    //  visible stem reflectance
  ArrayD1 rhosnir;    //  nir stem reflectance
  ArrayD1 taulvis;    //  visible leaf transmittance
  ArrayD1 taulnir;    //  nir leaf transmittance
  ArrayD1 tausvis;    //  visible stem transmittance
  ArrayD1 tausnir;    //  nir stem transmittance

  // default constructor
  VegData();
  // default destructor
  ~VegData() {};

  // Read pft time-invariant file data into member variables
  void read_veg_data(const std::string &data_dir, const std::string &basename_pfts);

  // get struct of photosynthesis variables for vegetation == vegtype
  const PSNVegData get_pft_psnveg (int vegtype) const;

  // get struct of albedo variables for vegetation == vegtype
  const AlbedoVegData get_pft_albveg(int vegtype) const;

};



} // namespace ELM

#include "vegdata_impl.hh"
