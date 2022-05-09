
#include <iostream>
#include <string>

// utilities
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"

// constants
#include "elm_constants.h"

// input data readers and structs
#include "pft_data.h"
#include "atm_data.h"
#include "soil_data.h"
#include "snicar_data.h"
//#include "satellite_phenology.h"

//#include "InitSoil.hh"
#include "init_soil_state.h"
#include "init_snow_state.h"


//#include "InitSnowLayers.hh"
//#include "InitTimestep.hh"
#include "init_timestep.h"
//#include "InitTopography.hh"
#include "init_topography.h"
#include "land_data.h"
//#include "read_atmosphere.h"
//#include "ReadTestData.hh"

#include "incident_shortwave.h"
#include "day_length.h"

#include "canopy_hydrology.h"
#include "surface_radiation.h"
#include "canopy_temperature.h"
#include "bareground_fluxes.h"
#include "canopy_fluxes.h"
#include "aerosol_data.h"
#include "aerosol_physics.h"
#include "phenology_data.h"
#include "surface_albedo.h"
#include "snow_snicar.h"
//#include "root_biophys.h"
#include "surface_fluxes.h"
#include "soil_texture_hydraulic_model.h"

#include "Kokkos_Core.hpp"

using ArrayB1 = Kokkos::View<bool *>;
using ArrayI1 = Kokkos::View<int *>;
using ArrayD1 = Kokkos::View<double *>;
using ArrayD2 = Kokkos::View<double **>;
using ArrayP1 = Kokkos::View<ELM::PSNVegData *>;

template <class Array_t> Array_t create(const std::string &name, int D0)
{ return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1)
{ return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2)
{ return Array_t(name, D0, D1, D2); }
template <class Array_t, typename T> void assign(Array_t &arr, T val)
{ Kokkos::deep_copy(arr, val); }


using AtmForcType = ELM::AtmForcType;

template<AtmForcType ftype>
using atm_forc_util = ELM::AtmDataManager<ArrayD1, ArrayD2, ftype>;

template <AtmForcType ftype>
atm_forc_util<ftype> create_forc_util(const std::string& filename,
  const ELM::Utils::Date &file_start_time, const int ntimes,
  const int ncells)
{ return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }


std::map<std::string, ArrayD2::HostMirror> get_phen_host_views(const ELM::PhenologyDataManager& phen_data)
{
  std::map<std::string, ArrayD2::HostMirror> phen_host_views;
  phen_host_views["MONTHLY_LAI"] = Kokkos::create_mirror_view(phen_data.mlai_);
  phen_host_views["MONTHLY_SAI"] = Kokkos::create_mirror_view(phen_data.msai_);
  phen_host_views["MONTHLY_HEIGHT_TOP"] = Kokkos::create_mirror_view(phen_data.mhtop_);
  phen_host_views["MONTHLY_HEIGHT_BOT"] = Kokkos::create_mirror_view(phen_data.mhbot_);
  return phen_host_views;
}

std::map<std::string, ArrayD1::HostMirror> get_aero_host_views(const ELM::AerosolDataManager& aero_data)
{
  std::map<std::string, ArrayD1::HostMirror> aero_host_views;
  aero_host_views["BCDEPWET"] = Kokkos::create_mirror_view(aero_data.bcdep_);
  aero_host_views["BCPHODRY"] = Kokkos::create_mirror_view(aero_data.bcpho_);
  aero_host_views["BCPHIDRY"] = Kokkos::create_mirror_view(aero_data.bcphi_);
  aero_host_views["DSTX01DD"] = Kokkos::create_mirror_view(aero_data.dst1_1_);
  aero_host_views["DSTX02DD"] = Kokkos::create_mirror_view(aero_data.dst1_2_);
  aero_host_views["DSTX03DD"] = Kokkos::create_mirror_view(aero_data.dst2_1_);
  aero_host_views["DSTX04DD"] = Kokkos::create_mirror_view(aero_data.dst2_2_);
  aero_host_views["DSTX01WD"] = Kokkos::create_mirror_view(aero_data.dst3_1_);
  aero_host_views["DSTX02WD"] = Kokkos::create_mirror_view(aero_data.dst3_2_);
  aero_host_views["DSTX03WD"] = Kokkos::create_mirror_view(aero_data.dst4_1_);
  aero_host_views["DSTX04WD"] = Kokkos::create_mirror_view(aero_data.dst4_2_);
  return aero_host_views;
}

void copy_aero_host_views(std::map<std::string, ArrayD1::HostMirror>& aero_host_views, ELM::AerosolDataManager& aero_data)
{
  Kokkos::deep_copy(aero_data.bcdep_, aero_host_views["BCDEPWET"]);
  Kokkos::deep_copy(aero_data.bcpho_, aero_host_views["BCPHODRY"]);
  Kokkos::deep_copy(aero_data.bcphi_, aero_host_views["BCPHIDRY"]);
  Kokkos::deep_copy(aero_data.dst1_1_, aero_host_views["DSTX01DD"]);
  Kokkos::deep_copy(aero_data.dst1_2_, aero_host_views["DSTX02DD"]);
  Kokkos::deep_copy(aero_data.dst2_1_, aero_host_views["DSTX03DD"]);
  Kokkos::deep_copy(aero_data.dst2_2_, aero_host_views["DSTX04DD"]);
  Kokkos::deep_copy(aero_data.dst3_1_, aero_host_views["DSTX01WD"]);
  Kokkos::deep_copy(aero_data.dst3_2_, aero_host_views["DSTX02WD"]);
  Kokkos::deep_copy(aero_data.dst4_1_, aero_host_views["DSTX03WD"]);
  Kokkos::deep_copy(aero_data.dst4_2_, aero_host_views["DSTX04WD"]);
}

std::map<std::string, ArrayD1::HostMirror> get_snicar_host_views_d1(const ELM::SnicarData& snicar_data)
{
  std::map<std::string, ArrayD1::HostMirror> snicar_host_views_d1;
  snicar_host_views_d1["ss_alb_ocphil"] = Kokkos::create_mirror_view(snicar_data.ss_alb_oc1);
  snicar_host_views_d1["asm_prm_ocphil"] = Kokkos::create_mirror_view(snicar_data.asm_prm_oc1);
  snicar_host_views_d1["ext_cff_mss_ocphil"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_oc1);
  snicar_host_views_d1["ss_alb_ocphob"] = Kokkos::create_mirror_view(snicar_data.ss_alb_oc2);
  snicar_host_views_d1["asm_prm_ocphob"] = Kokkos::create_mirror_view(snicar_data.asm_prm_oc2);
  snicar_host_views_d1["ext_cff_mss_ocphob"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_oc2);
  snicar_host_views_d1["ss_alb_dust01"] = Kokkos::create_mirror_view(snicar_data.ss_alb_dst1);
  snicar_host_views_d1["asm_prm_dust01"] = Kokkos::create_mirror_view(snicar_data.asm_prm_dst1);
  snicar_host_views_d1["ext_cff_mss_dust01"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_dst1);
  snicar_host_views_d1["ss_alb_dust02"] = Kokkos::create_mirror_view(snicar_data.ss_alb_dst2);
  snicar_host_views_d1["asm_prm_dust02"] = Kokkos::create_mirror_view(snicar_data.asm_prm_dst2);
  snicar_host_views_d1["ext_cff_mss_dust02"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_dst2);
  snicar_host_views_d1["ss_alb_dust03"] = Kokkos::create_mirror_view(snicar_data.ss_alb_dst3);
  snicar_host_views_d1["asm_prm_dust03"] = Kokkos::create_mirror_view(snicar_data.asm_prm_dst3);
  snicar_host_views_d1["ext_cff_mss_dust03"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_dst3);
  snicar_host_views_d1["ss_alb_dust04"] = Kokkos::create_mirror_view(snicar_data.ss_alb_dst4);
  snicar_host_views_d1["asm_prm_dust04"] = Kokkos::create_mirror_view(snicar_data.asm_prm_dst4);
  snicar_host_views_d1["ext_cff_mss_dust04"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_dst4);
  return snicar_host_views_d1;
}

std::map<std::string, ArrayD2::HostMirror> get_snicar_host_views_d2(const ELM::SnicarData& snicar_data)
{
  std::map<std::string, ArrayD2::HostMirror> snicar_host_views_d2;
  snicar_host_views_d2["ss_alb_ice_drc"] = Kokkos::create_mirror_view(snicar_data.ss_alb_snw_drc);
  snicar_host_views_d2["asm_prm_ice_drc"] = Kokkos::create_mirror_view(snicar_data.asm_prm_snw_drc);
  snicar_host_views_d2["ext_cff_mss_ice_drc"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_snw_drc);
  snicar_host_views_d2["ss_alb_ice_dfs"] = Kokkos::create_mirror_view(snicar_data.ss_alb_snw_dfs);
  snicar_host_views_d2["asm_prm_ice_dfs"] = Kokkos::create_mirror_view(snicar_data.asm_prm_snw_dfs);
  snicar_host_views_d2["ext_cff_mss_ice_dfs"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_snw_dfs);
  snicar_host_views_d2["ss_alb_bc_mam"] = Kokkos::create_mirror_view(snicar_data.ss_alb_bc1);
  snicar_host_views_d2["asm_prm_bc_mam"] = Kokkos::create_mirror_view(snicar_data.asm_prm_bc1);
  snicar_host_views_d2["ext_cff_mss_bc_mam"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_bc1);
  snicar_host_views_d2["ss_alb_bc_mam"] = Kokkos::create_mirror_view(snicar_data.ss_alb_bc2);
  snicar_host_views_d2["asm_prm_bc_mam"] = Kokkos::create_mirror_view(snicar_data.asm_prm_bc2);
  snicar_host_views_d2["ext_cff_mss_bc_mam"] = Kokkos::create_mirror_view(snicar_data.ext_cff_mss_bc2);
  return snicar_host_views_d2;
}

std::map<std::string, ArrayD3::HostMirror> get_snicar_host_views_d3(const ELM::SnicarData& snicar_data)
{
  std::map<std::string, ArrayD3::HostMirror> snicar_host_views_d3;
  snicar_host_views_d3["bcint_enh_mam"] = Kokkos::create_mirror_view(snicar_data.bcenh);
  return snicar_host_views_d3;
}




void copy_snicar_host_views_d1(std::map<std::string, ArrayD1::HostMirror>& snicar_host_views_d1, ELM::SnicarData& snicar_data)
{
  Kokkos::deep_copy(snicar_data.ss_alb_oc1, snicar_host_views_d1["ss_alb_ocphil"]);
  Kokkos::deep_copy(snicar_data.asm_prm_oc1, snicar_host_views_d1["asm_prm_ocphil"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_oc1, snicar_host_views_d1["ext_cff_mss_ocphil"]);
  Kokkos::deep_copy(snicar_data.ss_alb_oc2, snicar_host_views_d1["ss_alb_ocphob"]);
  Kokkos::deep_copy(snicar_data.asm_prm_oc2, snicar_host_views_d1["asm_prm_ocphob"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_oc2, snicar_host_views_d1["ext_cff_mss_ocphob"]);
  Kokkos::deep_copy(snicar_data.ss_alb_dst1, snicar_host_views_d1["ss_alb_dust01"]);
  Kokkos::deep_copy(snicar_data.asm_prm_dst1, snicar_host_views_d1["asm_prm_dust01"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_dst1, snicar_host_views_d1["ext_cff_mss_dust01"]);
  Kokkos::deep_copy(snicar_data.ss_alb_dst2, snicar_host_views_d1["ss_alb_dust02"]);
  Kokkos::deep_copy(snicar_data.asm_prm_dst2, snicar_host_views_d1["asm_prm_dust02"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_dst2, snicar_host_views_d1["ext_cff_mss_dust02"]);
  Kokkos::deep_copy(snicar_data.ss_alb_dst3, snicar_host_views_d1["ss_alb_dust03"]);
  Kokkos::deep_copy(snicar_data.asm_prm_dst3, snicar_host_views_d1["asm_prm_dust03"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_dst3, snicar_host_views_d1["ext_cff_mss_dust03"]);
  Kokkos::deep_copy(snicar_data.ss_alb_dst4, snicar_host_views_d1["ss_alb_dust04"]);
  Kokkos::deep_copy(snicar_data.asm_prm_dst4, snicar_host_views_d1["asm_prm_dust04"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_dst4, snicar_host_views_d1["ext_cff_mss_dust04"]);
}


void copy_snicar_host_views_d2(std::map<std::string, ArrayD2::HostMirror>& snicar_host_views_d2, ELM::SnicarData& snicar_data)
{
  Kokkos::deep_copy(snicar_data.ss_alb_snw_drc, snicar_host_views_d2["ss_alb_ice_drc"]);
  Kokkos::deep_copy(snicar_data.asm_prm_snw_drc, snicar_host_views_d2["asm_prm_ice_drc"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_snw_drc, snicar_host_views_d2["ext_cff_mss_ice_drc"]);
  Kokkos::deep_copy(snicar_data.ss_alb_snw_dfs, snicar_host_views_d2["ss_alb_ice_dfs"]);
  Kokkos::deep_copy(snicar_data.asm_prm_snw_dfs, snicar_host_views_d2["asm_prm_ice_dfs"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_snw_dfs, snicar_host_views_d2["ext_cff_mss_ice_dfs"]);
  Kokkos::deep_copy(snicar_data.ss_alb_bc1, snicar_host_views_d2["ss_alb_bc_mam"]);
  Kokkos::deep_copy(snicar_data.asm_prm_bc1, snicar_host_views_d2["asm_prm_bc_mam"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_bc1, snicar_host_views_d2["ext_cff_mss_bc_mam"]);
  Kokkos::deep_copy(snicar_data.ss_alb_bc2, snicar_host_views_d2["ss_alb_bc_mam"]);
  Kokkos::deep_copy(snicar_data.asm_prm_bc2, snicar_host_views_d2["asm_prm_bc_mam"]);
  Kokkos::deep_copy(snicar_data.ext_cff_mss_bc2, snicar_host_views_d2["ext_cff_mss_bc_mam"]);
}


void copy_snicar_host_views_d3(std::map<std::string, ArrayD3::HostMirror>& snicar_host_views_d3, ELM::SnicarData& snicar_data)
{
  Kokkos::deep_copy(snicar_data.bcenh, snicar_host_views_d3["bcint_enh_mam"]);
}







std::map<std::string, ArrayD1::HostMirror> get_veg_host_views(const ELM::VegData<ArrayD1, ArrayD2>& veg_data)
{
  std::map<std::string, ArrayD1::HostMirror> veg_host_views;
  veg_host_views["fnr"] = Kokkos::create_mirror_view(veg_data.fnr);
  veg_host_views["act25"] = Kokkos::create_mirror_view(veg_data.act25);
  veg_host_views["kcha"] = Kokkos::create_mirror_view(veg_data.kcha);
  veg_host_views["koha"] = Kokkos::create_mirror_view(veg_data.koha);
  veg_host_views["cpha"] = Kokkos::create_mirror_view(veg_data.cpha);
  veg_host_views["vcmaxha"] = Kokkos::create_mirror_view(veg_data.vcmaxha);
  veg_host_views["jmaxha"] = Kokkos::create_mirror_view(veg_data.jmaxha);
  veg_host_views["tpuha"] = Kokkos::create_mirror_view(veg_data.tpuha);
  veg_host_views["lmrha"] = Kokkos::create_mirror_view(veg_data.lmrha);
  veg_host_views["vcmaxhd"] = Kokkos::create_mirror_view(veg_data.vcmaxhd);
  veg_host_views["jmaxhd"] = Kokkos::create_mirror_view(veg_data.jmaxhd);
  veg_host_views["tpuhd"] = Kokkos::create_mirror_view(veg_data.tpuhd);
  veg_host_views["lmrhd"] = Kokkos::create_mirror_view(veg_data.lmrhd);
  veg_host_views["lmrse"] = Kokkos::create_mirror_view(veg_data.lmrse);
  veg_host_views["qe"] = Kokkos::create_mirror_view(veg_data.qe);
  veg_host_views["theta_cj"] = Kokkos::create_mirror_view(veg_data.theta_cj);
  veg_host_views["bbbopt"] = Kokkos::create_mirror_view(veg_data.bbbopt);
  veg_host_views["mbbopt"] = Kokkos::create_mirror_view(veg_data.mbbopt);
  veg_host_views["c3psn"] = Kokkos::create_mirror_view(veg_data.c3psn);
  veg_host_views["slatop"] = Kokkos::create_mirror_view(veg_data.slatop);
  veg_host_views["leafcn"] = Kokkos::create_mirror_view(veg_data.leafcn);
  veg_host_views["flnr"] = Kokkos::create_mirror_view(veg_data.flnr);
  veg_host_views["fnitr"] = Kokkos::create_mirror_view(veg_data.fnitr);
  veg_host_views["dleaf"] = Kokkos::create_mirror_view(veg_data.dleaf);
  veg_host_views["smpso"] = Kokkos::create_mirror_view(veg_data.smpso);
  veg_host_views["smpsc"] = Kokkos::create_mirror_view(veg_data.smpsc);
  veg_host_views["tc_stress"] = Kokkos::create_mirror_view(veg_data.tc_stress);
  veg_host_views["z0mr"] = Kokkos::create_mirror_view(veg_data.z0mr);
  veg_host_views["displar"] = Kokkos::create_mirror_view(veg_data.displar);
  veg_host_views["xl"] = Kokkos::create_mirror_view(veg_data.xl);
  veg_host_views["roota_par"] = Kokkos::create_mirror_view(veg_data.roota_par);
  veg_host_views["rootb_par"] = Kokkos::create_mirror_view(veg_data.rootb_par);
  veg_host_views["rholvis"] = Kokkos::create_mirror_view(veg_data.rholvis);
  veg_host_views["rholnir"] = Kokkos::create_mirror_view(veg_data.rholnir);
  veg_host_views["rhosvis"] = Kokkos::create_mirror_view(veg_data.rhosvis);
  veg_host_views["rhosnir"] = Kokkos::create_mirror_view(veg_data.rhosnir);
  veg_host_views["taulvis"] = Kokkos::create_mirror_view(veg_data.taulvis);
  veg_host_views["taulnir"] = Kokkos::create_mirror_view(veg_data.taulnir);
  veg_host_views["tausvis"] = Kokkos::create_mirror_view(veg_data.tausvis);
  veg_host_views["tausnir"] = Kokkos::create_mirror_view(veg_data.tausnir);
  return veg_host_views;
}


void copy_veg_host_views(std::map<std::string, ArrayD1::HostMirror>& veg_host_views, ELM::VegData<ArrayD1, ArrayD2>& veg_data)
{
  Kokkos::deep_copy(veg_data.fnr, veg_host_views["fnr"]);
  Kokkos::deep_copy(veg_data.act25, veg_host_views["act25"]);
  Kokkos::deep_copy(veg_data.kcha, veg_host_views["kcha"]);
  Kokkos::deep_copy(veg_data.koha, veg_host_views["koha"]);
  Kokkos::deep_copy(veg_data.cpha, veg_host_views["cpha"]);
  Kokkos::deep_copy(veg_data.vcmaxha, veg_host_views["vcmaxha"]);
  Kokkos::deep_copy(veg_data.jmaxha, veg_host_views["jmaxha"]);
  Kokkos::deep_copy(veg_data.tpuha, veg_host_views["tpuha"]);
  Kokkos::deep_copy(veg_data.lmrha, veg_host_views["lmrha"]);
  Kokkos::deep_copy(veg_data.vcmaxhd, veg_host_views["vcmaxhd"]);
  Kokkos::deep_copy(veg_data.jmaxhd, veg_host_views["jmaxhd"]);
  Kokkos::deep_copy(veg_data.tpuhd, veg_host_views["tpuhd"]);
  Kokkos::deep_copy(veg_data.lmrhd, veg_host_views["lmrhd"]);
  Kokkos::deep_copy(veg_data.lmrse, veg_host_views["lmrse"]);
  Kokkos::deep_copy(veg_data.qe, veg_host_views["qe"]);
  Kokkos::deep_copy(veg_data.theta_cj, veg_host_views["theta_cj"]);
  Kokkos::deep_copy(veg_data.bbbopt, veg_host_views["bbbopt"]);
  Kokkos::deep_copy(veg_data.mbbopt, veg_host_views["mbbopt"]);
  Kokkos::deep_copy(veg_data.c3psn, veg_host_views["c3psn"]);
  Kokkos::deep_copy(veg_data.slatop, veg_host_views["slatop"]);
  Kokkos::deep_copy(veg_data.leafcn, veg_host_views["leafcn"]);
  Kokkos::deep_copy(veg_data.flnr, veg_host_views["flnr"]);
  Kokkos::deep_copy(veg_data.fnitr, veg_host_views["fnitr"]);
  Kokkos::deep_copy(veg_data.dleaf, veg_host_views["dleaf"]);
  Kokkos::deep_copy(veg_data.smpso, veg_host_views["smpso"]);
  Kokkos::deep_copy(veg_data.smpsc, veg_host_views["smpsc"]);
  Kokkos::deep_copy(veg_data.tc_stress, veg_host_views["tc_stress"]);
  Kokkos::deep_copy(veg_data.z0mr, veg_host_views["z0mr"]);
  Kokkos::deep_copy(veg_data.displar, veg_host_views["displar"]);
  Kokkos::deep_copy(veg_data.xl, veg_host_views["xl"]);
  Kokkos::deep_copy(veg_data.roota_par, veg_host_views["roota_par"]);
  Kokkos::deep_copy(veg_data.rootb_par, veg_host_views["rootb_par"]);
  Kokkos::deep_copy(veg_data.rholvis, veg_host_views["rholvis"]);
  Kokkos::deep_copy(veg_data.rholnir, veg_host_views["rholnir"]);
  Kokkos::deep_copy(veg_data.rhosvis, veg_host_views["rhosvis"]);
  Kokkos::deep_copy(veg_data.rhosnir, veg_host_views["rhosnir"]);
  Kokkos::deep_copy(veg_data.taulvis, veg_host_views["taulvis"]);
  Kokkos::deep_copy(veg_data.taulnir, veg_host_views["taulnir"]);
  Kokkos::deep_copy(veg_data.tausvis, veg_host_views["tausvis"]);
  Kokkos::deep_copy(veg_data.tausnir, veg_host_views["tausnir"]);
}





int main(int argc, char **argv) {

{ // enclosing scope

  Kokkos::initialize(argc, argv);

  { // inner scope

    std::string fname_surfdata(
    "/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc");
    std::string fname_snicar(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
    std::string fname_forc(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc");
    std::string fname_pft(
      "/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata/lnd/clm2/paramdata/clm_params_c180524.nc");
    std::string fname_aerosol(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc");

    int MPI_COMM_WORLD;
    const int n_procs = 1;
    const int ncells = 1;
    const int ntimes = 1008;
    const int myrank = 0;
    const double dtime = 1800.0;
    const double dtime_d = 1800.0 / 86400.0;
    const auto start = ELM::Utils::Date(2014, 7, 15);

    auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
    auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { 1, 1 },
          { 0, 0 });

    // hardwired params
    const double lat = 71.323;
    const double lon = 203.3886;
    const double lat_r = lat * ELM::ELM_PI / 180.0;
    const double lon_r = lon * ELM::ELM_PI / 180.0;
    ELM::LandType Land;
    Land.ltype = 1;
    Land.ctype = 1;
    Land.vtype = 12;
    // hardwired pft type
    // this is the format phenology reader wants
    auto vtype = create<ArrayI1>("vtype", ncells);
    assign(vtype, 12);
    const double dewmx = 0.1;
    const double irrig_rate = 0.0;
    const int n_irrig_steps_left = 0;
    bool do_capsnow = false;
    const bool lakpoi = false;
    auto veg_active = create<ArrayB1>("veg_active", ncells); // need value
    assign(veg_active, true);                                      // hardwired

    // forcing data
    auto forc_tbot = create<ArrayD1>("forc_tbot", ncells);
    auto forc_thbot = create<ArrayD1>("forc_thbot", ncells);
    auto forc_pbot = create<ArrayD1>("forc_pbot", ncells);
    auto forc_qbot = create<ArrayD1>("forc_qbot", ncells);
    auto forc_rh = create<ArrayD1>("forc_rh", ncells);
    auto forc_lwrad = create<ArrayD1>("forc_lwrad", ncells);
    auto forc_solai = create<ArrayD2>("forc_solai", ncells, 2);
    auto forc_solad = create<ArrayD2>("forc_solad", ncells, 2);
    auto forc_rain = create<ArrayD1>("forc_rain", ncells);
    auto forc_snow = create<ArrayD1>("forc_snow", ncells);
    auto forc_u = create<ArrayD1>("forc_u", ncells);
    auto forc_v = create<ArrayD1>("forc_v", ncells);
    auto forc_hgt = create<ArrayD1>("forc_hgt", ncells);
    auto forc_hgt_u = create<ArrayD1>("forc_hgt_u", ncells);
    auto forc_hgt_t = create<ArrayD1>("forc_hgt_t", ncells);
    auto forc_hgt_q = create<ArrayD1>("forc_hgt_q", ncells);
    auto forc_vp = create<ArrayD1>("forc_vp", ncells);
    auto forc_rho = create<ArrayD1>("forc_rho", ncells);
    auto forc_po2 = create<ArrayD1>("forc_po2", ncells);
    auto forc_pco2 = create<ArrayD1>("forc_pco2", ncells);    

    // prescribed sat phenology
    auto tlai = create<ArrayD1>("tlai", ncells);
    auto tsai = create<ArrayD1>("tsai", ncells);
    auto elai = create<ArrayD1>("elai", ncells);
    auto esai = create<ArrayD1>("esai", ncells);
    auto htop = create<ArrayD1>("htop", ncells);
    auto hbot = create<ArrayD1>("hbot", ncells);
    auto frac_veg_nosno_alb = create<ArrayI1>("frac_veg_nosno_alb", ncells);

    // soil hydraulics
    auto watsat = create<ArrayD2>("watsat", ncells, ELM::nlevgrnd);
    auto sucsat = create<ArrayD2>("sucsat", ncells, ELM::nlevgrnd);
    auto bsw = create<ArrayD2>("bsw", ncells, ELM::nlevgrnd);
    auto watdry = create<ArrayD2>("watdry", ncells, ELM::nlevgrnd);
    auto watopt = create<ArrayD2>("watopt", ncells, ELM::nlevgrnd);
    auto watfc = create<ArrayD2>("watfc", ncells, ELM::nlevgrnd);
    
    // topo, microtopography
    auto n_melt = create<ArrayD1>("n_melt", ncells);
    auto micro_sigma = create<ArrayD1>("micro_sigma", ncells);
    auto topo_slope = create<ArrayD1>("topo_slope", ncells);
    assign(topo_slope, 0.070044865858546);
    auto topo_std = create<ArrayD1>("topo_std", ncells);
    assign(topo_std, 3.96141847422387);


    // soil color and texture constants
    auto isoicol = create<ArrayI1>("isoicol", ncells);
    auto albsat = create<ArrayD2>("albsat", ncells, 2);
    auto albdry = create<ArrayD2>("albdry", ncells, 2);
    auto pct_sand = create<ArrayD2>("pct_sand", ncells, ELM::nlevgrnd);
    auto pct_clay = create<ArrayD2>("pct_clay", ncells, ELM::nlevgrnd);
    auto organic = create<ArrayD2>("organic", ncells, ELM::nlevgrnd);

    // snow variables
    auto snl = create<ArrayI1>("snl", ncells);
    assign(snl, 0);
    auto snow_depth = create<ArrayD1>("snow_depth", ncells); // NEED VALUES! - probably always init at 0
    assign(snow_depth, 0.0);
    auto frac_sno = create<ArrayD1>("frac_sno", ncells);     // NEED VALUES!  \ if not glc, icemec, etc, always init these @ 0.0
    assign(frac_sno, 0.0);
    auto int_snow = create<ArrayD1>("int_snow", ncells);     // NEED VALUES!
    assign(int_snow, 0.0);


    // uncategorized
    auto t_grnd = create<ArrayD1>("t_grnd", ncells);
    auto h2ocan = create<ArrayD1>("h2ocan", ncells);
    auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", ncells);
    auto frac_iceold = create<ArrayD2>("frac_iceold", ncells, ELM::nlevsno + ELM::nlevgrnd);
    auto h2osno = create<ArrayD1>("h2osno", ncells);
    auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", ncells, ELM::nlevsno + ELM::nlevgrnd);
    assign(h2osoi_liq, 0.0);
    auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", ncells, ELM::nlevsno + ELM::nlevgrnd);
    assign(h2osoi_ice, 0.0);
    auto snw_rds = create<ArrayD2>("snw_rds", ncells, ELM::nlevsno);


    // for Canopy hydrology
    auto qflx_prec_grnd = create<ArrayD1>("qflx_prec_grnd", ncells);
    auto qflx_snwcp_liq = create<ArrayD1>("qflx_snwcp_liq", ncells);
    auto qflx_snwcp_ice = create<ArrayD1>("qflx_snwcp_ice", ncells);
    assign(qflx_snwcp_ice, 0.0);
    auto qflx_snow_grnd = create<ArrayD1>("qflx_snow_grnd", ncells);
    auto qflx_rain_grnd = create<ArrayD1>("qflx_rain_grnd", ncells);
    auto fwet = create<ArrayD1>("fwet", ncells);
    auto fdry = create<ArrayD1>("fdry", ncells);
    auto qflx_snow_melt = create<ArrayD1>("qflx_snow_melt", ncells);
    auto h2osfc = create<ArrayD1>("h2osfc", ncells);
    auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", ncells);
    auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", ncells);
    auto swe_old = create<ArrayD2>("swe_old", ncells, ELM::nlevsno);
    
    auto t_soisno = create<ArrayD2>("t_soisno", ncells, ELM::nlevsno + ELM::nlevgrnd);
    double tsoi[] = {0.0, 0.0, 0.0, 0.0, 0.0, 278.3081064745931, 276.1568781897738,
      275.55803480737063, 275.2677090940866, 274.7286996980052, 273.15, 272.4187794248787, 270.65049816473027,
      267.8224112387398, 265.7450135695632, 264.49481140089864, 264.14163363048056, 264.3351872934207, 264.1163763444719, 263.88852987294865};


    // for can_sun_shade
    auto nrad = create<ArrayI1>("nrad", ncells);
    auto laisun = create<ArrayD1>("laisun", ncells);
    auto laisha = create<ArrayD1>("laisha", ncells);
    auto tlai_z = create<ArrayD2>("tlai_z", ncells, ELM::nlevcan);
    auto fsun_z = create<ArrayD2>("fsun_z", ncells, ELM::nlevcan);
    auto fabd_sun_z = create<ArrayD2>("fabd_sun_z", ncells, ELM::nlevcan);
    auto fabd_sha_z = create<ArrayD2>("fabd_sha_z", ncells, ELM::nlevcan);
    auto fabi_sun_z = create<ArrayD2>("fabi_sun_z", ncells, ELM::nlevcan);
    auto fabi_sha_z = create<ArrayD2>("fabi_sha_z", ncells, ELM::nlevcan);
    auto parsun_z = create<ArrayD2>("parsun_z", ncells, ELM::nlevcan);
    auto parsha_z = create<ArrayD2>("parsha_z", ncells, ELM::nlevcan);
    auto laisun_z = create<ArrayD2>("laisun_z", ncells, ELM::nlevcan);
    auto laisha_z = create<ArrayD2>("laisha_z", ncells, ELM::nlevcan);

    // for surface rad
    auto sabg_soil = create<ArrayD1>("sabg_soil", ncells);
    auto sabg_snow = create<ArrayD1>("sabg_snow", ncells);
    auto sabg = create<ArrayD1>("sabg", ncells);
    auto sabv = create<ArrayD1>("sabv", ncells);
    auto fsa = create<ArrayD1>("fsa", ncells);
    auto fsr = create<ArrayD1>("fsr", ncells);
    auto sabg_lyr = create<ArrayD2>("sabg_lyr", ncells, ELM::nlevsno + 1);
    auto ftdd = create<ArrayD2>("ftdd", ncells, ELM::numrad);
    auto ftid = create<ArrayD2>("ftid", ncells, ELM::numrad);
    auto ftii = create<ArrayD2>("ftii", ncells, ELM::numrad);
    auto fabd = create<ArrayD2>("fabd", ncells, ELM::numrad);
    auto fabi = create<ArrayD2>("fabi", ncells, ELM::numrad);
    auto albsod = create<ArrayD2>("albsod", ncells, ELM::numrad);
    auto albsoi = create<ArrayD2>("albsoi", ncells, ELM::numrad);
    auto albsnd_hst = create<ArrayD2>("albsnd_hst", ncells, ELM::numrad);
    auto albsni_hst = create<ArrayD2>("albsni_hst", ncells, ELM::numrad);
    auto albgrd = create<ArrayD2>("albgrd", ncells, ELM::numrad);
    auto albgri = create<ArrayD2>("albgri", ncells, ELM::numrad);
    auto flx_absdv = create<ArrayD2>("flx_absdv", ncells, ELM::nlevsno + 1);
    auto flx_absdn = create<ArrayD2>("flx_absdn", ncells, ELM::nlevsno + 1);
    auto flx_absiv = create<ArrayD2>("flx_absiv", ncells, ELM::nlevsno + 1);
    auto flx_absin = create<ArrayD2>("flx_absin", ncells, ELM::nlevsno + 1);
    auto albd = create<ArrayD2>("albd", ncells, ELM::numrad);
    auto albi = create<ArrayD2>("albi", ncells, ELM::numrad);



    // variables for CanopyTemperature
    auto t_h2osfc = create<ArrayD1>("t_h2osfc", ncells);
    assign(t_h2osfc, 274.0);
    auto t_h2osfc_bef = create<ArrayD1>("t_h2osfc_bef", ncells);
    auto z_0_town = create<ArrayD1>("z_0_town", ncells);
    auto z_d_town = create<ArrayD1>("z_d_town", ncells);
    auto soilalpha = create<ArrayD1>("soilalpha", ncells);
    auto soilalpha_u = create<ArrayD1>("soilalpha_u", ncells);
    auto soilbeta = create<ArrayD1>("soilbeta", ncells);
    auto qg_snow = create<ArrayD1>("qg_snow", ncells);
    auto qg_soil = create<ArrayD1>("qg_soil", ncells);
    auto qg = create<ArrayD1>("qg", ncells);
    auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", ncells);
    auto dqgdT = create<ArrayD1>("dqgdT", ncells);
    auto htvp = create<ArrayD1>("htvp", ncells);
    auto emg = create<ArrayD1>("emg", ncells);
    auto emv = create<ArrayD1>("emv", ncells);
    auto z0mg = create<ArrayD1>("z0mg", ncells);
    auto z0hg = create<ArrayD1>("z0hg", ncells);
    auto z0qg = create<ArrayD1>("z0qg", ncells);
    auto z0mv = create<ArrayD1>("z0mv", ncells);
    auto z0hv = create<ArrayD1>("z0hv", ncells);
    auto z0qv = create<ArrayD1>("z0qv", ncells);
    auto thv = create<ArrayD1>("thv", ncells);
    auto z0m = create<ArrayD1>("z0m", ncells);
    auto displa = create<ArrayD1>("displa", ncells);
    auto thm = create<ArrayD1>("thm", ncells);
    auto eflx_sh_tot = create<ArrayD1>("eflx_sh_tot", ncells);
    auto eflx_sh_tot_u = create<ArrayD1>("eflx_sh_tot_u", ncells);
    auto eflx_sh_tot_r = create<ArrayD1>("eflx_sh_tot_r", ncells);
    auto eflx_lh_tot = create<ArrayD1>("eflx_lh_tot", ncells);
    auto eflx_lh_tot_u = create<ArrayD1>("eflx_lh_tot_u", ncells);
    auto eflx_lh_tot_r = create<ArrayD1>("eflx_lh_tot_r", ncells);
    auto eflx_sh_veg = create<ArrayD1>("eflx_sh_veg", ncells);
    auto qflx_evap_tot = create<ArrayD1>("qflx_evap_tot", ncells);
    auto qflx_evap_veg = create<ArrayD1>("qflx_evap_veg", ncells);
    auto qflx_tran_veg = create<ArrayD1>("qflx_tran_veg", ncells);
    auto tssbef = create<ArrayD2>("tssbef", ncells, ELM::nlevgrnd + ELM::nlevsno);
    auto rootfr_road_perv =
        create<ArrayD2>("rootfr_road_perv", ncells, ELM::nlevgrnd); // comes from SoilStateType.F90
    auto rootr_road_perv =
        create<ArrayD2>("rootr_road_perv", ncells, ELM::nlevgrnd); // comes from SoilStateType.F90

    // bareground fluxes
    auto dlrad = create<ArrayD1>("dlrad", ncells);
    auto ulrad = create<ArrayD1>("ulrad", ncells);
    auto eflx_sh_grnd = create<ArrayD1>("eflx_sh_grnd", ncells);
    assign(eflx_sh_grnd, 0.0);
    auto eflx_sh_snow = create<ArrayD1>("eflx_sh_snow", ncells);
    assign(eflx_sh_snow, 0.0);
    auto eflx_sh_soil = create<ArrayD1>("eflx_sh_soil", ncells);
    assign(eflx_sh_soil, 0.0);
    auto eflx_sh_h2osfc = create<ArrayD1>("eflx_sh_h2osfc", ncells);
    assign(eflx_sh_h2osfc, 0.0);
    auto qflx_evap_soi = create<ArrayD1>("qflx_evap_soi", ncells);
    assign(qflx_evap_soi, 0.0);
    auto qflx_ev_snow = create<ArrayD1>("qflx_ev_snow", ncells);
    assign(qflx_ev_snow, 0.0);
    auto qflx_ev_soil = create<ArrayD1>("qflx_ev_soil", ncells);
    assign(qflx_ev_soil, 0.0);
    auto qflx_ev_h2osfc = create<ArrayD1>("qflx_ev_h2osfc", ncells);
    assign(qflx_ev_h2osfc, 0.0);
    auto t_ref2m = create<ArrayD1>("t_ref2m", ncells);
    auto t_ref2m_r = create<ArrayD1>("t_ref2m_r", ncells);
    auto q_ref2m = create<ArrayD1>("q_ref2m", ncells);
    auto rh_ref2m = create<ArrayD1>("rh_ref2m", ncells);
    auto rh_ref2m_r = create<ArrayD1>("rh_ref2m_r", ncells);
    auto cgrnds = create<ArrayD1>("cgrnds", ncells);
    auto cgrndl = create<ArrayD1>("cgrndl", ncells);
    auto cgrnd = create<ArrayD1>("cgrnd", ncells);

    // canopy fluxes
    auto altmax_indx = create<ArrayI1>("altmax_indx", ncells);
    assign(altmax_indx, 5);
    auto altmax_lastyear_indx = create<ArrayI1>("altmax_lastyear_indx", ncells);
    assign(altmax_lastyear_indx, 0);
    auto t10 = create<ArrayD1>("t10", ncells);
    assign(t10, 276.0);
    auto vcmaxcintsha = create<ArrayD1>("vcmaxcintsha", ncells);
    auto vcmaxcintsun = create<ArrayD1>("vcmaxcintsun", ncells);
    auto btran = create<ArrayD1>("btran", ncells);
    auto t_veg = create<ArrayD1>("t_veg", ncells);
    assign(t_veg, 283.0);
    auto rootfr = create<ArrayD2>("rootfr", ncells, ELM::nlevgrnd);
    auto rootr = create<ArrayD2>("rootr", ncells, ELM::nlevgrnd);
    auto eff_porosity = create<ArrayD2>("eff_porosity", ncells, ELM::nlevgrnd);

    // surface albedo and snicar
    // required for SurfaceAlbedo kernels
    // I1
    auto snl_top = create<ArrayI1>("snl_top", ncells);
    auto snl_btm = create<ArrayI1>("snl_btm", ncells);
    auto ncan = create<ArrayI1>("ncan", ncells);
    auto flg_nosnl = create<ArrayI1>("flg_nosnl", ncells);
    // I2
    auto snw_rds_lcl = create<ArrayI2>("snw_rds_lcl", ncells, ELM::nlevsno);
    // D1
    auto mu_not = create<ArrayD1>("mu_not", ncells);
    // D2
    auto fabd_sun = create<ArrayD2>("fabd_sun", ncells, ELM::numrad);
    auto fabd_sha = create<ArrayD2>("fabd_sha", ncells, ELM::numrad);
    auto fabi_sun = create<ArrayD2>("fabi_sun", ncells, ELM::numrad);
    auto fabi_sha = create<ArrayD2>("fabi_sha", ncells, ELM::numrad);
    auto albsnd = create<ArrayD2>("albsnd", ncells, ELM::numrad);
    auto albsni = create<ArrayD2>("albsni", ncells, ELM::numrad);
    auto tsai_z = create<ArrayD2>("tsai_z", ncells, ELM::nlevcan);
    auto h2osoi_vol = create<ArrayD2>("h2osoi_vol", ncells, ELM::nlevgrnd);
    // D3
    auto mss_cnc_aer_in_fdb = create<ArrayD3>("mss_cnc_aer_in_fdb", ncells, ELM::nlevsno, ELM::snow_snicar::sno_nbr_aer);
    auto flx_absd_snw = create<ArrayD3>("flx_absd_snw", ncells, ELM::nlevsno+1, ELM::numrad);
    auto flx_absi_snw = create<ArrayD3>("flx_absi_snw", ncells, ELM::nlevsno+1, ELM::numrad);
    auto flx_abs_lcl = create<ArrayD3>("flx_abs_lcl", ncells, ELM::nlevsno+1, ELM::snow_snicar::numrad_snw);
    // D2
    auto albout_lcl = create<ArrayD2>("albout_lcl", ncells, ELM::snow_snicar::numrad_snw);
    auto flx_slrd_lcl = create<ArrayD2>("flx_slrd_lcl", ncells, ELM::snow_snicar::numrad_snw);
    auto flx_slri_lcl = create<ArrayD2>("flx_slri_lcl", ncells, ELM::snow_snicar::numrad_snw);
    auto h2osoi_ice_lcl = create<ArrayD2>("h2osoi_ice_lcl", ncells, ELM::nlevsno);
    auto h2osoi_liq_lcl = create<ArrayD2>("h2osoi_liq_lcl", ncells, ELM::nlevsno);
    // D3
    auto g_star = create<ArrayD3>("g_star", ncells, ELM::snow_snicar::numrad_snw, ELM::nlevsno);
    auto omega_star = create<ArrayD3>("omega_star", ncells, ELM::snow_snicar::numrad_snw, ELM::nlevsno);
    auto tau_star = create<ArrayD3>("tau_star", ncells, ELM::snow_snicar::numrad_snw, ELM::nlevsno);

    // soil fluxes (outputs)
    auto eflx_soil_grnd = create<ArrayD1>("eflx_soil_grnd", ncells);
    auto qflx_evap_grnd = create<ArrayD1>("qflx_evap_grnd", ncells);
    auto qflx_sub_snow = create<ArrayD1>("qflx_sub_snow", ncells);
    auto qflx_dew_snow = create<ArrayD1>("qflx_dew_snow", ncells);
    auto qflx_dew_grnd = create<ArrayD1>("qflx_dew_grnd", ncells);

    // grid data 
    auto dz = create<ArrayD2>("dz", ncells, ELM::nlevsno + ELM::nlevgrnd);
    auto zsoi = create<ArrayD2>("zsoi", ncells, ELM::nlevsno + ELM::nlevgrnd);
    auto zisoi = create<ArrayD2>("zisoi", ncells, ELM::nlevsno + ELM::nlevgrnd + 1);

    // hardwired grid info
    // this comes from ELM, but is wrong?
    // doesn't matter for now, but definitely wrong
    {
      double dz_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.017512817916255204,
      0.02757896925967625, 0.0454700332424132, 0.07496741098620856,
      0.12360036510228053, 0.20378255101043175, 0.33598062644843263,
      0.5539384053686849, 0.9132900315890611, 1.5057607013992766,
      2.482579696981332, 4.0930819526214, 6.7483512780057175,
      11.12615029420442, 13.851152141963599 };
      auto h_dz = Kokkos::create_mirror_view(dz);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < ELM::nlevsno + ELM::nlevgrnd; ++i) {
          h_dz(n, i) = dz_hardwire[i];
        }
      }
      Kokkos::deep_copy(dz, h_dz);

      double zsoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.007100635417193535,
      0.02792500041531687, 0.06225857393654604, 0.11886506690014327,
      0.21219339590896316, 0.3660657971047043, 0.6197584979298266,
      1.0380270500015696, 1.7276353086671965, 2.8646071131796917,
      4.73915671146575, 7.829766507142356, 12.92532061670855,
      21.32646906315379, 35.17762120511739 };
      auto h_zsoi = Kokkos::create_mirror_view(zsoi);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < ELM::nlevsno + ELM::nlevgrnd; ++i) {
          h_zsoi(n, i) = zsoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(zsoi, h_zsoi);

      double zisoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.017512817916255204, 0.04509178717593146, 0.09056182041834465, 
      0.16552923140455322, 0.28912959650683373, 0.4929121475172655,
      0.8288927739656982, 1.382831179334383, 2.2961212109234443,
      3.8018819123227208, 6.284461609304053, 10.377543561925453,
      17.12589483993117, 28.252045134135592, 42.10319727609919 };
      auto h_zisoi = Kokkos::create_mirror_view(zisoi);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < ELM::nlevsno + ELM::nlevgrnd + 1; ++i) {
          h_zisoi(n, i) = zisoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(zisoi, h_zisoi);
    }
    





    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // init data containers and read time-invariant data from files
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

    // need to modify !!
    int atm_nsteps = 1000;
    auto test_TBOT = create_forc_util<AtmForcType::TBOT>(fname_forc, start, atm_nsteps, ncells);
    auto test_PBOT = create_forc_util<AtmForcType::PBOT>(fname_forc, start, atm_nsteps, ncells);
    auto test_QBOT = create_forc_util<AtmForcType::RH>(fname_forc, start, atm_nsteps, ncells);
    auto test_FLDS = create_forc_util<AtmForcType::FLDS>(fname_forc, start, atm_nsteps, ncells);
    auto test_FSDS = create_forc_util<AtmForcType::FSDS>(fname_forc, start, atm_nsteps, ncells);
    auto test_PREC = create_forc_util<AtmForcType::PREC>(fname_forc, start, atm_nsteps, ncells);
    auto test_WIND = create_forc_util<AtmForcType::WIND>(fname_forc, start, atm_nsteps, ncells);
    auto test_ZBOT = create_forc_util<AtmForcType::ZBOT>(fname_forc, start, atm_nsteps, ncells);
    
    {
      // soil color constants
      auto h_isoicol = Kokkos::create_mirror_view(isoicol);
      auto h_albsat = Kokkos::create_mirror_view(albsat);
      auto h_albdry = Kokkos::create_mirror_view(albdry);
      ELM::read_soil::read_soil_colors(dd, fname_surfdata, h_isoicol, h_albsat, h_albdry);
      if (albsat.extent(0) != h_albsat.extent(0) || albsat.extent(1) != h_albsat.extent(1))
      {
        NS::resize(albsat, h_albsat.extent(0), h_albsat.extent(1));
      }
      if (albdry.extent(0) != h_albdry.extent(0) || albdry.extent(1) != h_albdry.extent(1))
      {
        NS::resize(albdry, h_albdry.extent(0), h_albdry.extent(1));
      }
      Kokkos::deep_copy(isoicol, h_isoicol);
      Kokkos::deep_copy(albsat, h_albsat);
      Kokkos::deep_copy(albdry, h_albdry);
    }

    {
      // soil texture constants
      auto h_pct_sand = Kokkos::create_mirror_view(pct_sand);
      auto h_pct_clay = Kokkos::create_mirror_view(pct_clay);
      auto h_organic = Kokkos::create_mirror_view(organic);
      ELM::read_soil::read_soil_texture(dd, fname_surfdata, h_pct_sand, h_pct_clay, h_organic);
      Kokkos::deep_copy(pct_sand, h_pct_sand);
      Kokkos::deep_copy(pct_clay, h_pct_clay);
      Kokkos::deep_copy(organic, h_organic);
    }

    // snicar radiation parameters
    ELM::SnicarData snicar_data;
    {
      auto host_snicar_d1 = get_snicar_host_views_d1(snicar_data);
      auto host_snicar_d2 = get_snicar_host_views_d2(snicar_data);
      auto host_snicar_d3 = get_snicar_host_views_d3(snicar_data);
      ELM::read_snicar_data(host_snicar_d1, host_snicar_d2,
                            host_snicar_d3, dd.comm, fname_snicar);
      copy_snicar_host_views_d1(host_snicar_d1, snicar_data);
      copy_snicar_host_views_d2(host_snicar_d2, snicar_data);
      copy_snicar_host_views_d3(host_snicar_d3, snicar_data);
    }

    // pft data constants
    ELM::VegData<ArrayD1, ArrayD2> vegdata;
    {
      auto host_veg_views = get_veg_host_views(vegdata);
      vegdata.read_veg_data(host_veg_views, dd.comm, fname_pft);
      copy_veg_host_views(host_veg_views, vegdata);
    }

    // aerosol deposition data manager
    ELM::AerosolDataManager aerosol_data;
    {
      auto host_aero_views = get_aero_host_views(aerosol_data);
      aerosol_data.read_data(host_aero_views, dd.comm, fname_aerosol, lon, lat);
      copy_aero_host_views(host_aero_views, aerosol_data);
    }

    // phenology data manager
    // make host mirrors - need to be persistent
    ELM::PhenologyDataManager phen_data(dd, ncells, 17);
    auto host_phen_views = get_phen_host_views(phen_data);

    // need to modify !!
    //ELM::AerosolMasses aerosol_masses(ncells);
    //ELM::AerosolConcentrations aerosol_concentrations(ncells);







    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // some free functions that calculate initial values
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    
    ELM::init_topo::init_topo_slope(topo_slope[0]);

    ELM::init_topo::init_micro_topo(Land.ltype, topo_slope[0], topo_std[0],
                                    n_melt[0], micro_sigma[0]);

    ELM::init_snow_state::init_snow_layers(snow_depth[0], lakpoi, snl[0],
                                           Kokkos::subview(dz, 0, Kokkos::ALL),
                                           Kokkos::subview(zsoi, 0, Kokkos::ALL),
                                           Kokkos::subview(zisoi, 0, Kokkos::ALL));

    ELM::soil_hydraulics::init_soil_hydraulics(
                                           Kokkos::subview(pct_sand, 0, Kokkos::ALL),
                                           Kokkos::subview(pct_clay, 0, Kokkos::ALL),
                                           Kokkos::subview(organic, 0, Kokkos::ALL),
                                           Kokkos::subview(zsoi, 0, Kokkos::ALL),
                                           Kokkos::subview(watsat, 0, Kokkos::ALL),
                                           Kokkos::subview(bsw, 0, Kokkos::ALL),
                                           Kokkos::subview(sucsat, 0, Kokkos::ALL),
                                           Kokkos::subview(watdry, 0, Kokkos::ALL),
                                           Kokkos::subview(watopt, 0, Kokkos::ALL),
                                           Kokkos::subview(watfc, 0, Kokkos::ALL));


    ELM::init_soil_state::init_vegrootfr(vtype[0], vegdata.roota_par[vtype[0]],
                                         vegdata.rootb_par[vtype[0]],
                                         Kokkos::subview(zisoi, 0, Kokkos::ALL),
                                         Kokkos::subview(rootfr, 0, Kokkos::ALL));

    ELM::init_soil_state::init_soil_temp(Land, snl[0], 
                                         Kokkos::subview(t_soisno, 0, Kokkos::ALL),
                                         t_grnd[0]);

    ELM::init_snow_state::init_snow_state(Land.urbpoi, snl[0], h2osno[0], int_snow[0],
                                          snow_depth[0], h2osfc[0], h2ocan[0], frac_h2osfc[0],
                                          fwet[0], fdry[0], frac_sno[0],
                                          Kokkos::subview(snw_rds, 0, Kokkos::ALL));

    ELM::init_soil_state::init_soilh2o_state(Land, snl[0], 
                                             Kokkos::subview(watsat, 0, Kokkos::ALL),
                                             Kokkos::subview(t_soisno, 0, Kokkos::ALL),
                                             Kokkos::subview(dz, 0, Kokkos::ALL),
                                             Kokkos::subview(h2osoi_vol, 0, Kokkos::ALL),
                                             Kokkos::subview(h2osoi_liq, 0, Kokkos::ALL),
                                             Kokkos::subview(h2osoi_ice, 0, Kokkos::ALL));




    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /*                          TIME LOOP                                                                  */
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // these calls should be inside time loop
    test_TBOT.read_atm_forcing(dd, start, atm_nsteps);
    test_PBOT.read_atm_forcing(dd, start, atm_nsteps);
    test_QBOT.read_atm_forcing(dd, start, atm_nsteps);
    test_FLDS.read_atm_forcing(dd, start, atm_nsteps);
    test_FSDS.read_atm_forcing(dd, start, atm_nsteps);
    test_PREC.read_atm_forcing(dd, start, atm_nsteps);
    test_WIND.read_atm_forcing(dd, start, atm_nsteps);
    test_ZBOT.read_atm_forcing(dd, start, atm_nsteps);

    auto coszen = create<ArrayD1>("coszen", ncells);
    auto cosz_factor = create<ArrayD1>("cosz_factor", ncells);

    ELM::Utils::Date current(start);
    ELM::Utils::Date big(start);
    int idx = 0; // hardwire for ncells = 1
    
    for (int t = 0; t < ntimes; ++t) {
      
      ELM::Utils::Date time_plus_half_dt(current);
      time_plus_half_dt.increment_seconds(900); // assumed 1800s dt

      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // get coszen - should make this a function
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      double max_dayl;
      double dayl;
      {
        ELM::Utils::Date forc_dt_start{test_FSDS.get_data_start_time()};
        forc_dt_start.increment_seconds(round(test_FSDS.forc_t_idx(time_plus_half_dt, test_FSDS.get_data_start_time()) * test_FSDS.get_forc_dt_secs()));
    
        int cosz_doy = current.doy + 1;
        double declin = ELM::incident_shortwave::declination_angle2(cosz_doy); // should be the same for model and forcing - don't cross day barrier in forcing timeseries
        double cosz_decday = ELM::Utils::decimal_doy(current) + 1.0;
        double cosz_forc_decday = ELM::Utils::decimal_doy(forc_dt_start) + 1.0;
        double cosz_forcdt_avg = ELM::incident_shortwave::average_cosz(lat_r, lon_r, declin, test_FSDS.get_forc_dt_secs(), cosz_forc_decday);
        double thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, cosz_decday + dtime_d/2.0);

        cosz_factor(0) = (thiscosz > 0.001) ? 1.0 : 0.0;
        coszen(0) = thiscosz;
        max_dayl = ELM::max_daylength(lat_r);
        dayl = ELM::daylength(lat_r, declin);
      }




      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // timestep init functions
      // read time-variable data
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/


      // read phenology data if required
      // reader will read 3 months of data on first call
      // subsequent calls only read the newest months (when phen_data.need_data() == true)
      // and shift the index of the two remaining older months

      // copy device data to host
      // copying entire views is likely inefficient, but it's currently necessary
      // could be eliminated by shifting older month indices in parallel kernel
      // and reading new data into a mirror of a subview (or a subview of a mirror?)
      // then we would only need one copy from host view into the device view
      // instead of the two we currently have
      // will fix later - too infrequently run to cause concern  
      if (phen_data.need_data()) {
        Kokkos::deep_copy(host_phen_views["MONTHLY_LAI"], phen_data.mlai_);
        Kokkos::deep_copy(host_phen_views["MONTHLY_SAI"], phen_data.msai_);
        Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_TOP"], phen_data.mhtop_);
        Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_BOT"], phen_data.mhbot_);
      }
      // reads three months of data on first call
      // after first call, read new data if phen_data.need_new_data_ == true
      auto phen_updated = phen_data.read_data(host_phen_views, fname_surfdata, current, vtype); // if needed
      // copy host views to device
      // could be made more efficient, see above
      if (phen_updated) {
        Kokkos::deep_copy(phen_data.mlai_, host_phen_views["MONTHLY_LAI"]);
        Kokkos::deep_copy(phen_data.msai_, host_phen_views["MONTHLY_SAI"]);
        Kokkos::deep_copy(phen_data.mhtop_, host_phen_views["MONTHLY_HEIGHT_TOP"]);
        Kokkos::deep_copy(phen_data.mhbot_, host_phen_views["MONTHLY_HEIGHT_BOT"]);
      }
      // run parallel kernel to process phenology data
      phen_data.get_data(current, snow_depth,
                         frac_sno, vtype, elai, esai,
                         htop, hbot, tlai, tsai,
                         frac_veg_nosno_alb);


      //test_TBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_thbot);
      //test_PBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot);
      //test_QBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_pbot, forc_qbot, forc_rh);
      //test_FLDS.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot, forc_qbot, forc_tbot, forc_lwrad);
      //test_FSDS.get_atm_forcing(dtime_d, time_plus_half_dt, cosz_factor, forc_solai, forc_solad);
      //test_PREC.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_rain, forc_snow);
      //test_WIND.get_atm_forcing(dtime_d, time_plus_half_dt, forc_u, forc_v);
      //test_ZBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_hgt, forc_hgt_u, forc_hgt_t,  forc_hgt_q);

      //// calc constitutive air props - call properly later
      //ELM::atm_forcing_physics::ConstitutiveAirProperties air_constitutive(forc_qbot, forc_pbot, 
      //                         forc_tbot, forc_vp,
      //                         forc_rho, forc_po2, forc_pco2);
      //air_constitutive(0); // - call properly

      // get aerosol mss and cnc
      //aerosol_data.invoke_aerosol_source(time_plus_half_dt, dtime, snl, aerosol_masses);
      //ELM::aerosols::invoke_aerosol_concen_and_mass(do_capsnow, dtime, snl, h2osoi_liq,
      //h2osoi_ice, snw_rds, qflx_snwcp_ice, aerosol_masses, aerosol_concentrations);

// only here for testing
auto h2osno_old = create<ArrayD1>("h2osno_old", ncells);   // get rid of this - not needed
auto eflx_bot = create<ArrayD1>("eflx_bot", ncells);       // don't need this, either
auto qflx_glcice = create<ArrayD1>("qflx_glcice", ncells); // don't need this, either
      
      ELM::init_timestep(lakpoi, h2osno[0], veg_active[0], snl[0],
                         Kokkos::subview(h2osoi_ice, 0, Kokkos::ALL),
                         Kokkos::subview(h2osoi_liq, 0, Kokkos::ALL),
                         frac_veg_nosno_alb[0], h2osno_old[0],
                         do_capsnow, eflx_bot[0], qflx_glcice[0],
                         frac_veg_nosno[0],
                         Kokkos::subview(frac_iceold, 0, Kokkos::ALL));

      

      for (int i = 0; i < ncells; ++i) {
      
        std::cout << "elai: " << elai(i) << std::endl;
        std::cout << "esai: " << esai(i) << std::endl;
        std::cout << "tlai: " << tlai(i) << std::endl;
        std::cout << "tsai: " << tsai(i) << std::endl;
        std::cout << "htop: " << htop(i) << std::endl;
        std::cout << "hbot: " << hbot(i) << std::endl;
        std::cout << "frac_veg_nosno_alb: " << frac_veg_nosno_alb(i) << std::endl;
      }

      current.increment_seconds(1800);
    }

  } // inner scope

  Kokkos::finalize();
} // enclosing scope

return 0;
}
