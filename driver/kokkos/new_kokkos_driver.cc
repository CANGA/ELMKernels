
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
#include "land_data.h"
#include "pft_data.h"
#include "atm_data.h"
#include "soil_data.h"
#include "snicar_data.h"
#include "aerosol_data.h"
#include "phenology_data.h"

// initialization routines
#include "init_soil_state.h"
#include "init_snow_state.h"
#include "init_timestep.h"
#include "init_topography.h"

// physics kernels
#include "day_length.h" // serial

#include "incident_shortwave.h"
#include "canopy_hydrology.h"
#include "surface_radiation.h"
#include "canopy_temperature.h"
#include "bareground_fluxes.h"
#include "canopy_fluxes.h"
#include "aerosol_physics.h"
#include "surface_albedo.h"
#include "snow_snicar.h"
#include "surface_fluxes.h"
#include "soil_texture_hydraulic_model.h"
#include "soil_temperature.h"
#include "soil_temp_rhs.h"
#include "soil_temp_lhs.h"
#include "soil_thermal_properties.h"
#include "pentadiagonal_solver.h"
#include "snow_hydrology.h"
#include "transpiration_impl.hh"

#include "helper_functions.hh"

// conditional compilation options
#include "invoke_kernel.hh"
#include "kokkos_includes.hh"

#include "Kokkos_Core.hpp"


using ViewB1 =  Kokkos::View<bool *>;
using ViewI1 =  Kokkos::View<int *>;
using ViewI2 =  Kokkos::View<int **>;
using ViewD1 =  Kokkos::View<double *>;
using ViewD2 =  Kokkos::View<double **>;
using ViewD3 =  Kokkos::View<double ***>;
using h_ViewD1 = ViewD1::HostMirror;
using h_ViewD2 = ViewD2::HostMirror;
using h_ViewD3 = ViewD3::HostMirror;

using ELM::Utils::create;
using ELM::Utils::assign;

using AtmForcType = ELM::AtmForcType;

template<AtmForcType ftype>
using atm_forc_util = ELM::AtmDataManager<ViewD1, ViewD2, ftype>;

template <AtmForcType ftype>
atm_forc_util<ftype> create_forc_util(const std::string& filename,
  const ELM::Utils::Date &file_start_time, const int ntimes,
  const int ncells)
{ return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }


std::map<std::string, h_ViewD2> get_phen_host_views(const ELM::PhenologyDataManager<ViewD2>& phen_data)
{
  std::map<std::string, h_ViewD2> phen_host_views;
  phen_host_views["MONTHLY_LAI"] = Kokkos::create_mirror_view(phen_data.mlai);
  phen_host_views["MONTHLY_SAI"] = Kokkos::create_mirror_view(phen_data.msai);
  phen_host_views["MONTHLY_HEIGHT_TOP"] = Kokkos::create_mirror_view(phen_data.mhtop);
  phen_host_views["MONTHLY_HEIGHT_BOT"] = Kokkos::create_mirror_view(phen_data.mhbot);
  return phen_host_views;
}

std::map<std::string, h_ViewD1> get_aero_host_views(const ELM::AerosolDataManager<ViewD1>& aero_data)
{
  std::map<std::string, h_ViewD1> aero_host_views;
  aero_host_views["BCDEPWET"] = Kokkos::create_mirror_view(aero_data.bcdep);
  aero_host_views["BCPHODRY"] = Kokkos::create_mirror_view(aero_data.bcpho);
  aero_host_views["BCPHIDRY"] = Kokkos::create_mirror_view(aero_data.bcphi);
  aero_host_views["DSTX01DD"] = Kokkos::create_mirror_view(aero_data.dst1_1);
  aero_host_views["DSTX02DD"] = Kokkos::create_mirror_view(aero_data.dst1_2);
  aero_host_views["DSTX03DD"] = Kokkos::create_mirror_view(aero_data.dst2_1);
  aero_host_views["DSTX04DD"] = Kokkos::create_mirror_view(aero_data.dst2_2);
  aero_host_views["DSTX01WD"] = Kokkos::create_mirror_view(aero_data.dst3_1);
  aero_host_views["DSTX02WD"] = Kokkos::create_mirror_view(aero_data.dst3_2);
  aero_host_views["DSTX03WD"] = Kokkos::create_mirror_view(aero_data.dst4_1);
  aero_host_views["DSTX04WD"] = Kokkos::create_mirror_view(aero_data.dst4_2);
  return aero_host_views;
}

void copy_aero_host_views(std::map<std::string, h_ViewD1>& aero_host_views, ELM::AerosolDataManager<ViewD1>& aero_data)
{
  Kokkos::deep_copy(aero_data.bcdep, aero_host_views["BCDEPWET"]);
  Kokkos::deep_copy(aero_data.bcpho, aero_host_views["BCPHODRY"]);
  Kokkos::deep_copy(aero_data.bcphi, aero_host_views["BCPHIDRY"]);
  Kokkos::deep_copy(aero_data.dst1_1, aero_host_views["DSTX01DD"]);
  Kokkos::deep_copy(aero_data.dst1_2, aero_host_views["DSTX02DD"]);
  Kokkos::deep_copy(aero_data.dst2_1, aero_host_views["DSTX03DD"]);
  Kokkos::deep_copy(aero_data.dst2_2, aero_host_views["DSTX04DD"]);
  Kokkos::deep_copy(aero_data.dst3_1, aero_host_views["DSTX01WD"]);
  Kokkos::deep_copy(aero_data.dst3_2, aero_host_views["DSTX02WD"]);
  Kokkos::deep_copy(aero_data.dst4_1, aero_host_views["DSTX03WD"]);
  Kokkos::deep_copy(aero_data.dst4_2, aero_host_views["DSTX04WD"]);
}

std::map<std::string, h_ViewD1> get_snicar_host_views_d1(const ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  std::map<std::string, h_ViewD1> snicar_host_views_d1;
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

std::map<std::string, h_ViewD2> get_snicar_host_views_d2(const ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  std::map<std::string, h_ViewD2> snicar_host_views_d2;
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

std::map<std::string, h_ViewD3> get_snicar_host_views_d3(const ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  std::map<std::string, h_ViewD3> snicar_host_views_d3;
  snicar_host_views_d3["bcint_enh_mam"] = Kokkos::create_mirror_view(snicar_data.bcenh);
  return snicar_host_views_d3;
}

std::map<std::string, h_ViewD3> get_snowage_host_views_d3(const ELM::SnwRdsTable<ViewD3>& snw_table)
{
  std::map<std::string, h_ViewD3> snowage_host_views_d3;
  snowage_host_views_d3["tau"] = Kokkos::create_mirror_view(snw_table.snowage_tau);
  snowage_host_views_d3["kappa"] = Kokkos::create_mirror_view(snw_table.snowage_kappa);
  snowage_host_views_d3["drdsdt0"] = Kokkos::create_mirror_view(snw_table.snowage_drdt0);

  return snowage_host_views_d3;
}




void copy_snicar_host_views_d1(std::map<std::string, h_ViewD1>& snicar_host_views_d1, ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
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


void copy_snicar_host_views_d2(std::map<std::string, h_ViewD2>& snicar_host_views_d2,
  ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
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


void copy_snicar_host_views_d3(std::map<std::string, h_ViewD3>& snicar_host_views_d3,
  ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  Kokkos::deep_copy(snicar_data.bcenh, snicar_host_views_d3["bcint_enh_mam"]);
}

void copy_snowage_host_views_d3(std::map<std::string, h_ViewD3>& snowage_host_views_d3,
  ELM::SnwRdsTable<ViewD3>& snw_table)
{
  Kokkos::deep_copy(snw_table.snowage_tau, snowage_host_views_d3["tau"]);
  Kokkos::deep_copy(snw_table.snowage_kappa, snowage_host_views_d3["kappa"]);
  Kokkos::deep_copy(snw_table.snowage_drdt0, snowage_host_views_d3["drdsdt0"]);
}





std::map<std::string, h_ViewD1> get_pft_host_views(const ELM::PFTData<ViewD1, ViewD2>& pft_data)
{
  std::map<std::string, h_ViewD1> pft_host_views;
  pft_host_views["fnr"] = Kokkos::create_mirror_view(pft_data.fnr);
  pft_host_views["act25"] = Kokkos::create_mirror_view(pft_data.act25);
  pft_host_views["kcha"] = Kokkos::create_mirror_view(pft_data.kcha);
  pft_host_views["koha"] = Kokkos::create_mirror_view(pft_data.koha);
  pft_host_views["cpha"] = Kokkos::create_mirror_view(pft_data.cpha);
  pft_host_views["vcmaxha"] = Kokkos::create_mirror_view(pft_data.vcmaxha);
  pft_host_views["jmaxha"] = Kokkos::create_mirror_view(pft_data.jmaxha);
  pft_host_views["tpuha"] = Kokkos::create_mirror_view(pft_data.tpuha);
  pft_host_views["lmrha"] = Kokkos::create_mirror_view(pft_data.lmrha);
  pft_host_views["vcmaxhd"] = Kokkos::create_mirror_view(pft_data.vcmaxhd);
  pft_host_views["jmaxhd"] = Kokkos::create_mirror_view(pft_data.jmaxhd);
  pft_host_views["tpuhd"] = Kokkos::create_mirror_view(pft_data.tpuhd);
  pft_host_views["lmrhd"] = Kokkos::create_mirror_view(pft_data.lmrhd);
  pft_host_views["lmrse"] = Kokkos::create_mirror_view(pft_data.lmrse);
  pft_host_views["qe"] = Kokkos::create_mirror_view(pft_data.qe);
  pft_host_views["theta_cj"] = Kokkos::create_mirror_view(pft_data.theta_cj);
  pft_host_views["bbbopt"] = Kokkos::create_mirror_view(pft_data.bbbopt);
  pft_host_views["mbbopt"] = Kokkos::create_mirror_view(pft_data.mbbopt);
  pft_host_views["c3psn"] = Kokkos::create_mirror_view(pft_data.c3psn);
  pft_host_views["slatop"] = Kokkos::create_mirror_view(pft_data.slatop);
  pft_host_views["leafcn"] = Kokkos::create_mirror_view(pft_data.leafcn);
  pft_host_views["flnr"] = Kokkos::create_mirror_view(pft_data.flnr);
  pft_host_views["fnitr"] = Kokkos::create_mirror_view(pft_data.fnitr);
  pft_host_views["dleaf"] = Kokkos::create_mirror_view(pft_data.dleaf);
  pft_host_views["smpso"] = Kokkos::create_mirror_view(pft_data.smpso);
  pft_host_views["smpsc"] = Kokkos::create_mirror_view(pft_data.smpsc);
  pft_host_views["tc_stress"] = Kokkos::create_mirror_view(pft_data.tc_stress);
  pft_host_views["z0mr"] = Kokkos::create_mirror_view(pft_data.z0mr);
  pft_host_views["displar"] = Kokkos::create_mirror_view(pft_data.displar);
  pft_host_views["xl"] = Kokkos::create_mirror_view(pft_data.xl);
  pft_host_views["roota_par"] = Kokkos::create_mirror_view(pft_data.roota_par);
  pft_host_views["rootb_par"] = Kokkos::create_mirror_view(pft_data.rootb_par);
  pft_host_views["rholvis"] = Kokkos::create_mirror_view(pft_data.rholvis);
  pft_host_views["rholnir"] = Kokkos::create_mirror_view(pft_data.rholnir);
  pft_host_views["rhosvis"] = Kokkos::create_mirror_view(pft_data.rhosvis);
  pft_host_views["rhosnir"] = Kokkos::create_mirror_view(pft_data.rhosnir);
  pft_host_views["taulvis"] = Kokkos::create_mirror_view(pft_data.taulvis);
  pft_host_views["taulnir"] = Kokkos::create_mirror_view(pft_data.taulnir);
  pft_host_views["tausvis"] = Kokkos::create_mirror_view(pft_data.tausvis);
  pft_host_views["tausnir"] = Kokkos::create_mirror_view(pft_data.tausnir);
  return pft_host_views;
}


void copy_pft_host_views(std::map<std::string, h_ViewD1>& pft_host_views, ELM::PFTData<ViewD1, ViewD2>& pft_data)
{
  Kokkos::deep_copy(pft_data.fnr, pft_host_views["fnr"]);
  Kokkos::deep_copy(pft_data.act25, pft_host_views["act25"]);
  Kokkos::deep_copy(pft_data.kcha, pft_host_views["kcha"]);
  Kokkos::deep_copy(pft_data.koha, pft_host_views["koha"]);
  Kokkos::deep_copy(pft_data.cpha, pft_host_views["cpha"]);
  Kokkos::deep_copy(pft_data.vcmaxha, pft_host_views["vcmaxha"]);
  Kokkos::deep_copy(pft_data.jmaxha, pft_host_views["jmaxha"]);
  Kokkos::deep_copy(pft_data.tpuha, pft_host_views["tpuha"]);
  Kokkos::deep_copy(pft_data.lmrha, pft_host_views["lmrha"]);
  Kokkos::deep_copy(pft_data.vcmaxhd, pft_host_views["vcmaxhd"]);
  Kokkos::deep_copy(pft_data.jmaxhd, pft_host_views["jmaxhd"]);
  Kokkos::deep_copy(pft_data.tpuhd, pft_host_views["tpuhd"]);
  Kokkos::deep_copy(pft_data.lmrhd, pft_host_views["lmrhd"]);
  Kokkos::deep_copy(pft_data.lmrse, pft_host_views["lmrse"]);
  Kokkos::deep_copy(pft_data.qe, pft_host_views["qe"]);
  Kokkos::deep_copy(pft_data.theta_cj, pft_host_views["theta_cj"]);
  Kokkos::deep_copy(pft_data.bbbopt, pft_host_views["bbbopt"]);
  Kokkos::deep_copy(pft_data.mbbopt, pft_host_views["mbbopt"]);
  Kokkos::deep_copy(pft_data.c3psn, pft_host_views["c3psn"]);
  Kokkos::deep_copy(pft_data.slatop, pft_host_views["slatop"]);
  Kokkos::deep_copy(pft_data.leafcn, pft_host_views["leafcn"]);
  Kokkos::deep_copy(pft_data.flnr, pft_host_views["flnr"]);
  Kokkos::deep_copy(pft_data.fnitr, pft_host_views["fnitr"]);
  Kokkos::deep_copy(pft_data.dleaf, pft_host_views["dleaf"]);
  Kokkos::deep_copy(pft_data.smpso, pft_host_views["smpso"]);
  Kokkos::deep_copy(pft_data.smpsc, pft_host_views["smpsc"]);
  Kokkos::deep_copy(pft_data.tc_stress, pft_host_views["tc_stress"]);
  Kokkos::deep_copy(pft_data.z0mr, pft_host_views["z0mr"]);
  Kokkos::deep_copy(pft_data.displar, pft_host_views["displar"]);
  Kokkos::deep_copy(pft_data.xl, pft_host_views["xl"]);
  Kokkos::deep_copy(pft_data.roota_par, pft_host_views["roota_par"]);
  Kokkos::deep_copy(pft_data.rootb_par, pft_host_views["rootb_par"]);
  Kokkos::deep_copy(pft_data.rholvis, pft_host_views["rholvis"]);
  Kokkos::deep_copy(pft_data.rholnir, pft_host_views["rholnir"]);
  Kokkos::deep_copy(pft_data.rhosvis, pft_host_views["rhosvis"]);
  Kokkos::deep_copy(pft_data.rhosnir, pft_host_views["rhosnir"]);
  Kokkos::deep_copy(pft_data.taulvis, pft_host_views["taulvis"]);
  Kokkos::deep_copy(pft_data.taulnir, pft_host_views["taulnir"]);
  Kokkos::deep_copy(pft_data.tausvis, pft_host_views["tausvis"]);
  Kokkos::deep_copy(pft_data.tausnir, pft_host_views["tausnir"]);
}


template <AtmForcType ftype>
void read_atm_data(ELM::AtmDataManager<ViewD1, ViewD2, ftype>& atm_data,
                   const ELM::Utils::DomainDecomposition<2>& dd,
                   const ELM::Utils::Date& model_time,
                   const size_t& ntimes)
{
  auto h_data = Kokkos::create_mirror_view(atm_data.data);
  atm_data.read_atm_forcing(h_data, dd, model_time, ntimes);

  if (atm_data.data.extent(0) != h_data.extent(0) || atm_data.data.extent(1) != h_data.extent(1))
    NS::resize(atm_data.data, h_data.extent(0), h_data.extent(1));

  Kokkos::deep_copy(atm_data.data, h_data);
}


int main(int argc, char **argv) {

{ // enclosing scope

  using namespace ELM::ELMdims;

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
    std::string fname_snowage(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc");
    int MPI_COMM_WORLD;
    const int n_procs = 1;
    const int ncells = 1;
    int idx = 0; // hardwire for ncells = 1
    const int ntimes = 3000;
    const int myrank = 0;
    const double dtime = 1800.0;
    const double dtime_d = 1800.0 / 86400.0;
    const auto start = ELM::Utils::Date(2014, 1, 1);

    auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
    auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { 1, 1 },
          { 0, 0 });

    // hardwired params
    const double lat = 71.323;
    const double lon = 203.3886;
    const double lat_r = lat * ELM::ELMconst::ELM_PI / 180.0;
    const double lon_r = lon * ELM::ELMconst::ELM_PI / 180.0;
    ELM::LandType Land;
    Land.ltype = 1;
    Land.ctype = 1;
    Land.vtype = 12;
    // hardwired pft type
    // this is the format phenology reader wants
    auto vtype = create<ViewI1>("vtype", ncells);
    assign(vtype, 12);
    const double dewmx = 0.1;
    const double irrig_rate = 0.0;
    const int n_irrig_steps_left = 0;
    const int oldfflag = 1;
    const bool lakpoi = false;
    auto veg_active = create<ViewB1>("veg_active", ncells); // need value
    assign(veg_active, true);                               // hardwired
    auto do_capsnow = create<ViewB1>("do_capsnow", ncells); // need value
    assign(do_capsnow, false);                               // hardwired

    // forcing data
    auto forc_tbot = create<ViewD1>("forc_tbot", ncells);
    auto forc_thbot = create<ViewD1>("forc_thbot", ncells);
    auto forc_pbot = create<ViewD1>("forc_pbot", ncells);
    auto forc_qbot = create<ViewD1>("forc_qbot", ncells);
    auto forc_rh = create<ViewD1>("forc_rh", ncells);
    auto forc_lwrad = create<ViewD1>("forc_lwrad", ncells);
    auto forc_solai = create<ViewD2>("forc_solai", ncells, 2);
    auto forc_solad = create<ViewD2>("forc_solad", ncells, 2);
    auto forc_rain = create<ViewD1>("forc_rain", ncells);
    auto forc_snow = create<ViewD1>("forc_snow", ncells);
    auto forc_u = create<ViewD1>("forc_u", ncells);
    auto forc_v = create<ViewD1>("forc_v", ncells);
    auto forc_hgt = create<ViewD1>("forc_hgt", ncells);
    auto forc_hgt_u = create<ViewD1>("forc_hgt_u", ncells);
    auto forc_hgt_t = create<ViewD1>("forc_hgt_t", ncells);
    auto forc_hgt_q = create<ViewD1>("forc_hgt_q", ncells);
    auto forc_vp = create<ViewD1>("forc_vp", ncells);
    auto forc_rho = create<ViewD1>("forc_rho", ncells);
    auto forc_po2 = create<ViewD1>("forc_po2", ncells);
    auto forc_pco2 = create<ViewD1>("forc_pco2", ncells);    

    // prescribed sat phenology
    auto tlai = create<ViewD1>("tlai", ncells);
    auto tsai = create<ViewD1>("tsai", ncells);
    auto elai = create<ViewD1>("elai", ncells);
    auto esai = create<ViewD1>("esai", ncells);
    auto htop = create<ViewD1>("htop", ncells);
    auto hbot = create<ViewD1>("hbot", ncells);
    auto frac_veg_nosno_alb = create<ViewI1>("frac_veg_nosno_alb", ncells);

    // soil hydraulics
    auto watsat = create<ViewD2>("watsat", ncells, nlevgrnd);
    auto sucsat = create<ViewD2>("sucsat", ncells, nlevgrnd);
    auto bsw = create<ViewD2>("bsw", ncells, nlevgrnd);
    auto watdry = create<ViewD2>("watdry", ncells, nlevgrnd);
    auto watopt = create<ViewD2>("watopt", ncells, nlevgrnd);
    auto watfc = create<ViewD2>("watfc", ncells, nlevgrnd);
    
    // topo, microtopography
    auto n_melt = create<ViewD1>("n_melt", ncells);
    auto micro_sigma = create<ViewD1>("micro_sigma", ncells);
    auto topo_slope = create<ViewD1>("topo_slope", ncells);
    assign(topo_slope, 0.070044865858546);
    auto topo_std = create<ViewD1>("topo_std", ncells);
    assign(topo_std, 3.96141847422387);


    // soil color and texture constants
    int max_soil_color = 20; // largest option - can't know until NC file is read
    auto isoicol = create<ViewI1>("isoicol", ncells);
    auto albsat = create<ViewD2>("albsat", max_soil_color, 2);
    auto albdry = create<ViewD2>("albdry", max_soil_color, 2);
    auto pct_sand = create<ViewD2>("pct_sand", ncells, nlevgrnd);
    auto pct_clay = create<ViewD2>("pct_clay", ncells, nlevgrnd);
    auto organic = create<ViewD2>("organic", ncells, nlevgrnd);

    // soil thermal constants
    auto tkmg = create<ViewD2>("tkmg", ncells, nlevgrnd);
    auto tkdry = create<ViewD2>("tkdry", ncells, nlevgrnd);
    auto csol = create<ViewD2>("csol", ncells, nlevgrnd + nlevsno);

    // snow variables
    auto snl = create<ViewI1>("snl", ncells);
    assign(snl, 0);
    auto snow_depth = create<ViewD1>("snow_depth", ncells); // NEED VALUES! - probably always init at 0
    assign(snow_depth, 0.0);
    auto frac_sno = create<ViewD1>("frac_sno", ncells);     // NEED VALUES!  \ if not glc, icemec, etc, always init these @ 0.0
    assign(frac_sno, 0.0);
    auto int_snow = create<ViewD1>("int_snow", ncells);     // NEED VALUES!
    assign(int_snow, 0.0);


    // uncategorized
    auto t_grnd = create<ViewD1>("t_grnd", ncells);
    auto h2ocan = create<ViewD1>("h2ocan", ncells);
    auto frac_veg_nosno = create<ViewI1>("frac_veg_nosno", ncells);
    auto frac_iceold = create<ViewD2>("frac_iceold", ncells, nlevsno + nlevgrnd);
    auto h2osno = create<ViewD1>("h2osno", ncells);
    auto h2osoi_liq = create<ViewD2>("h2osoi_liq", ncells, nlevsno + nlevgrnd);
    assign(h2osoi_liq, 0.0);
    auto h2osoi_ice = create<ViewD2>("h2osoi_ice", ncells, nlevsno + nlevgrnd);
    assign(h2osoi_ice, 0.0);
    auto snw_rds = create<ViewD2>("snw_rds", ncells, nlevsno);


    // for Canopy hydrology
    auto qflx_prec_grnd = create<ViewD1>("qflx_prec_grnd", ncells);
    auto qflx_snwcp_liq = create<ViewD1>("qflx_snwcp_liq", ncells);
    auto qflx_snwcp_ice = create<ViewD1>("qflx_snwcp_ice", ncells);
    assign(qflx_snwcp_ice, 0.0);
    auto qflx_snow_grnd = create<ViewD1>("qflx_snow_grnd", ncells);
    auto qflx_rain_grnd = create<ViewD1>("qflx_rain_grnd", ncells);
    auto fwet = create<ViewD1>("fwet", ncells);
    auto fdry = create<ViewD1>("fdry", ncells);
    auto qflx_snow_melt = create<ViewD1>("qflx_snow_melt", ncells);
    auto h2osfc = create<ViewD1>("h2osfc", ncells);
    auto frac_h2osfc = create<ViewD1>("frac_h2osfc", ncells);
    auto frac_sno_eff = create<ViewD1>("frac_sno_eff", ncells);
    auto swe_old = create<ViewD2>("swe_old", ncells, nlevsno);
    
    auto t_soisno = create<ViewD2>("t_soisno", ncells, nlevsno + nlevgrnd);

    // for can_sun_shade
    auto nrad = create<ViewI1>("nrad", ncells);
    auto laisun = create<ViewD1>("laisun", ncells);
    auto laisha = create<ViewD1>("laisha", ncells);
    auto tlai_z = create<ViewD2>("tlai_z", ncells, nlevcan);
    auto fsun_z = create<ViewD2>("fsun_z", ncells, nlevcan);
    auto fabd_sun_z = create<ViewD2>("fabd_sun_z", ncells, nlevcan);
    auto fabd_sha_z = create<ViewD2>("fabd_sha_z", ncells, nlevcan);
    auto fabi_sun_z = create<ViewD2>("fabi_sun_z", ncells, nlevcan);
    auto fabi_sha_z = create<ViewD2>("fabi_sha_z", ncells, nlevcan);
    auto parsun_z = create<ViewD2>("parsun_z", ncells, nlevcan);
    auto parsha_z = create<ViewD2>("parsha_z", ncells, nlevcan);
    auto laisun_z = create<ViewD2>("laisun_z", ncells, nlevcan);
    auto laisha_z = create<ViewD2>("laisha_z", ncells, nlevcan);

    // for surface rad
    auto sabg_soil = create<ViewD1>("sabg_soil", ncells);
    auto sabg_snow = create<ViewD1>("sabg_snow", ncells);
    auto sabg = create<ViewD1>("sabg", ncells);
    auto sabv = create<ViewD1>("sabv", ncells);
    auto fsa = create<ViewD1>("fsa", ncells);
    auto fsr = create<ViewD1>("fsr", ncells);
    auto sabg_lyr = create<ViewD2>("sabg_lyr", ncells, nlevsno + 1);
    auto ftdd = create<ViewD2>("ftdd", ncells, numrad);
    auto ftid = create<ViewD2>("ftid", ncells, numrad);
    auto ftii = create<ViewD2>("ftii", ncells, numrad);
    auto fabd = create<ViewD2>("fabd", ncells, numrad);
    auto fabi = create<ViewD2>("fabi", ncells, numrad);
    auto albsod = create<ViewD2>("albsod", ncells, numrad);
    auto albsoi = create<ViewD2>("albsoi", ncells, numrad);
    auto albsnd_hst = create<ViewD2>("albsnd_hst", ncells, numrad);
    auto albsni_hst = create<ViewD2>("albsni_hst", ncells, numrad);
    auto albgrd = create<ViewD2>("albgrd", ncells, numrad);
    auto albgri = create<ViewD2>("albgri", ncells, numrad);
    auto flx_absdv = create<ViewD2>("flx_absdv", ncells, nlevsno + 1);
    auto flx_absdn = create<ViewD2>("flx_absdn", ncells, nlevsno + 1);
    auto flx_absiv = create<ViewD2>("flx_absiv", ncells, nlevsno + 1);
    auto flx_absin = create<ViewD2>("flx_absin", ncells, nlevsno + 1);
    auto albd = create<ViewD2>("albd", ncells, numrad);
    auto albi = create<ViewD2>("albi", ncells, numrad);



    // variables for CanopyTemperature
    auto t_h2osfc = create<ViewD1>("t_h2osfc", ncells);
    assign(t_h2osfc, 274.0);
    auto t_h2osfc_bef = create<ViewD1>("t_h2osfc_bef", ncells);
    auto z_0_town = create<ViewD1>("z_0_town", ncells);
    auto z_d_town = create<ViewD1>("z_d_town", ncells);
    auto soilalpha = create<ViewD1>("soilalpha", ncells);
    auto soilalpha_u = create<ViewD1>("soilalpha_u", ncells);
    auto soilbeta = create<ViewD1>("soilbeta", ncells);
    auto qg_snow = create<ViewD1>("qg_snow", ncells);
    auto qg_soil = create<ViewD1>("qg_soil", ncells);
    auto qg = create<ViewD1>("qg", ncells);
    auto qg_h2osfc = create<ViewD1>("qg_h2osfc", ncells);
    auto dqgdT = create<ViewD1>("dqgdT", ncells);
    auto htvp = create<ViewD1>("htvp", ncells);
    auto emg = create<ViewD1>("emg", ncells);
    auto emv = create<ViewD1>("emv", ncells);
    auto z0mg = create<ViewD1>("z0mg", ncells);
    auto z0hg = create<ViewD1>("z0hg", ncells);
    auto z0qg = create<ViewD1>("z0qg", ncells);
    auto z0mv = create<ViewD1>("z0mv", ncells);
    auto z0hv = create<ViewD1>("z0hv", ncells);
    auto z0qv = create<ViewD1>("z0qv", ncells);
    auto thv = create<ViewD1>("thv", ncells);
    auto z0m = create<ViewD1>("z0m", ncells);
    auto displa = create<ViewD1>("displa", ncells);
    auto thm = create<ViewD1>("thm", ncells);
    auto eflx_sh_tot = create<ViewD1>("eflx_sh_tot", ncells);
    auto eflx_sh_tot_u = create<ViewD1>("eflx_sh_tot_u", ncells);
    auto eflx_sh_tot_r = create<ViewD1>("eflx_sh_tot_r", ncells);
    auto eflx_lh_tot = create<ViewD1>("eflx_lh_tot", ncells);
    auto eflx_lh_tot_u = create<ViewD1>("eflx_lh_tot_u", ncells);
    auto eflx_lh_tot_r = create<ViewD1>("eflx_lh_tot_r", ncells);
    auto eflx_sh_veg = create<ViewD1>("eflx_sh_veg", ncells);
    auto qflx_evap_tot = create<ViewD1>("qflx_evap_tot", ncells);
    auto qflx_evap_veg = create<ViewD1>("qflx_evap_veg", ncells);
    auto qflx_tran_veg = create<ViewD1>("qflx_tran_veg", ncells);
    auto tssbef = create<ViewD2>("tssbef", ncells, nlevgrnd + nlevsno);
    auto rootfr_road_perv =
        create<ViewD2>("rootfr_road_perv", ncells, nlevgrnd); // comes from SoilStateType.F90
    auto rootr_road_perv =
        create<ViewD2>("rootr_road_perv", ncells, nlevgrnd); // comes from SoilStateType.F90

    auto forc_hgt_u_patch = create<ViewD1>("forc_hgt_u", ncells);
    auto forc_hgt_t_patch = create<ViewD1>("forc_hgt_t", ncells);
    auto forc_hgt_q_patch = create<ViewD1>("forc_hgt_q", ncells);

    // bareground fluxes
    auto dlrad = create<ViewD1>("dlrad", ncells);
    auto ulrad = create<ViewD1>("ulrad", ncells);
    auto eflx_sh_grnd = create<ViewD1>("eflx_sh_grnd", ncells);
    assign(eflx_sh_grnd, 0.0);
    auto eflx_sh_snow = create<ViewD1>("eflx_sh_snow", ncells);
    assign(eflx_sh_snow, 0.0);
    auto eflx_sh_soil = create<ViewD1>("eflx_sh_soil", ncells);
    assign(eflx_sh_soil, 0.0);
    auto eflx_sh_h2osfc = create<ViewD1>("eflx_sh_h2osfc", ncells);
    assign(eflx_sh_h2osfc, 0.0);
    auto qflx_evap_soi = create<ViewD1>("qflx_evap_soi", ncells);
    assign(qflx_evap_soi, 0.0);
    auto qflx_ev_snow = create<ViewD1>("qflx_ev_snow", ncells);
    assign(qflx_ev_snow, 0.0);
    auto qflx_ev_soil = create<ViewD1>("qflx_ev_soil", ncells);
    assign(qflx_ev_soil, 0.0);
    auto qflx_ev_h2osfc = create<ViewD1>("qflx_ev_h2osfc", ncells);
    assign(qflx_ev_h2osfc, 0.0);
    auto t_ref2m = create<ViewD1>("t_ref2m", ncells);
    auto t_ref2m_r = create<ViewD1>("t_ref2m_r", ncells);
    auto q_ref2m = create<ViewD1>("q_ref2m", ncells);
    auto rh_ref2m = create<ViewD1>("rh_ref2m", ncells);
    auto rh_ref2m_r = create<ViewD1>("rh_ref2m_r", ncells);
    auto cgrnds = create<ViewD1>("cgrnds", ncells);
    auto cgrndl = create<ViewD1>("cgrndl", ncells);
    auto cgrnd = create<ViewD1>("cgrnd", ncells);

    // canopy fluxes
    auto altmax_indx = create<ViewI1>("altmax_indx", ncells);
    assign(altmax_indx, 5);
    auto altmax_lastyear_indx = create<ViewI1>("altmax_lastyear_indx", ncells);
    assign(altmax_lastyear_indx, 0);
    auto t10 = create<ViewD1>("t10", ncells);
    assign(t10, 276.0);
    auto vcmaxcintsha = create<ViewD1>("vcmaxcintsha", ncells);
    auto vcmaxcintsun = create<ViewD1>("vcmaxcintsun", ncells);
    auto btran = create<ViewD1>("btran", ncells);
    auto t_veg = create<ViewD1>("t_veg", ncells);
    assign(t_veg, 283.0);
    auto rootfr = create<ViewD2>("rootfr", ncells, nlevgrnd);
    auto rootr = create<ViewD2>("rootr", ncells, nlevgrnd);
    auto eff_porosity = create<ViewD2>("eff_porosity", ncells, nlevgrnd);

    // surface albedo and snicar
    // required for SurfaceAlbedo kernels
    // I1
    auto snl_top = create<ViewI1>("snl_top", ncells);
    auto snl_btm = create<ViewI1>("snl_btm", ncells);
    auto ncan = create<ViewI1>("ncan", ncells);
    auto flg_nosnl = create<ViewI1>("flg_nosnl", ncells);
    // I2
    auto snw_rds_lcl = create<ViewI2>("snw_rds_lcl", ncells, nlevsno);
    // D1
    auto mu_not = create<ViewD1>("mu_not", ncells);
    // D2
    auto fabd_sun = create<ViewD2>("fabd_sun", ncells, numrad);
    auto fabd_sha = create<ViewD2>("fabd_sha", ncells, numrad);
    auto fabi_sun = create<ViewD2>("fabi_sun", ncells, numrad);
    auto fabi_sha = create<ViewD2>("fabi_sha", ncells, numrad);
    auto albsnd = create<ViewD2>("albsnd", ncells, numrad);
    auto albsni = create<ViewD2>("albsni", ncells, numrad);
    auto tsai_z = create<ViewD2>("tsai_z", ncells, nlevcan);
    auto h2osoi_vol = create<ViewD2>("h2osoi_vol", ncells, nlevgrnd);
    // D3
    auto mss_cnc_aer_in_fdb = create<ViewD3>("mss_cnc_aer_in_fdb", ncells, nlevsno, sno_nbr_aer);
    auto flx_absd_snw = create<ViewD3>("flx_absd_snw", ncells, nlevsno+1, numrad);
    auto flx_absi_snw = create<ViewD3>("flx_absi_snw", ncells, nlevsno+1, numrad);
    auto flx_abs_lcl = create<ViewD3>("flx_abs_lcl", ncells, nlevsno+1, numrad_snw);
    // D2
    auto albout_lcl = create<ViewD2>("albout_lcl", ncells, numrad_snw);
    auto flx_slrd_lcl = create<ViewD2>("flx_slrd_lcl", ncells, numrad_snw);
    auto flx_slri_lcl = create<ViewD2>("flx_slri_lcl", ncells, numrad_snw);
    auto h2osoi_ice_lcl = create<ViewD2>("h2osoi_ice_lcl", ncells, nlevsno);
    auto h2osoi_liq_lcl = create<ViewD2>("h2osoi_liq_lcl", ncells, nlevsno);
    // D3
    auto g_star = create<ViewD3>("g_star", ncells, numrad_snw, nlevsno);
    auto omega_star = create<ViewD3>("omega_star", ncells, numrad_snw, nlevsno);
    auto tau_star = create<ViewD3>("tau_star", ncells, numrad_snw, nlevsno);

    // from soil temp - used in soil_e_balance 
    auto fact = create<ViewD2>("fact", ncells, nlevgrnd + nlevsno); // factors used in computing tridiagonal matrix
    // soil fluxes (outputs)
    auto eflx_soil_grnd = create<ViewD1>("eflx_soil_grnd", ncells);
    auto qflx_evap_grnd = create<ViewD1>("qflx_evap_grnd", ncells);
    auto qflx_sub_snow = create<ViewD1>("qflx_sub_snow", ncells);
    auto qflx_dew_snow = create<ViewD1>("qflx_dew_snow", ncells);
    auto qflx_dew_grnd = create<ViewD1>("qflx_dew_grnd", ncells);
    auto eflx_lwrad_out = create<ViewD1>("eflx_lwrad_out", ncells);
    auto eflx_lwrad_net = create<ViewD1>("eflx_lwrad_net", ncells);
    auto soil_e_balance = create<ViewD1>("soil_e_balance", ncells);

    // outputs from soil temp/snow hydro
    auto snot_top =  create<ViewD1>("snot_top", ncells);
    auto dTdz_top =  create<ViewD1>("dTdz_top", ncells);
    auto snw_rds_top =  create<ViewD1>("snw_rds_top", ncells);
    auto sno_liq_top =  create<ViewD1>("sno_liq_top", ncells);
    auto qflx_sl_top_soil = create<ViewD1>("qflx_sl_top_soil", ncells);
    auto qflx_snow2topsoi = create<ViewD1>("qflx_snow2topsoi", ncells);
    auto mflx_snowlyr_col = create<ViewD1>("mflx_snowlyr_col", ncells);
    auto qflx_top_soil = create<ViewD1>("qflx_top_soil", ncells);
    auto mflx_neg_snow = create<ViewD1>("mflx_neg_snow", ncells);
    auto eflx_snomelt = create<ViewD1>("eflx_snomelt", ncells);
    auto qflx_snomelt = create<ViewD1>("qflx_snomelt", ncells);
    auto xmf_dummy = create<ViewD1>("xmf", ncells);
    auto xmf_h2osfc_dummy = create<ViewD1>("xmf_h2osfc", ncells);
    auto eflx_h2osfc_to_snow_dummy = create<ViewD1>("eflx_h2osfc_to_snow", ncells);
    auto qflx_h2osfc_to_ice_dummy = create<ViewD1>("eflx_h2osfc_to_snow", ncells);
    auto eflx_building_heat_dummy = create<ViewD1>("eflx_building_heat", ncells);
    auto qflx_snofrz = create<ViewD1>("qflx_snofrz", ncells);
    auto qflx_snofrz_lyr = create<ViewD2>("qflx_snofrz_lyr", ncells, nlevsno);
    auto imelt = create<ViewI2>("imelt", ncells, nlevgrnd + nlevsno);
    assign(xmf_dummy, 0.0);
    assign(xmf_h2osfc_dummy, 0.0);
    assign(eflx_h2osfc_to_snow_dummy, 0.0);
    assign(eflx_building_heat_dummy, 0.0);

    // main exchange variables
    // transpiration
    // vegetation/soil water exchange (m H2O/s) (+ = to atm)
    auto qflx_rootsoi = create<ViewD2>("qflx_rootsoi", ncells, nlevgrnd);

    // grid data 
    auto dz = create<ViewD2>("dz", ncells, nlevsno + nlevgrnd);
    auto zsoi = create<ViewD2>("zsoi", ncells, nlevsno + nlevgrnd);
    auto zisoi = create<ViewD2>("zisoi", ncells, nlevsno + nlevgrnd + 1);

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
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
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
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
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
        for (int i = 0; i < nlevsno + nlevgrnd + 1; ++i) {
          h_zisoi(n, i) = zisoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(zisoi, h_zisoi);
    }


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // init data containers and read time-invariant data from files
    // these call only need to occur once @ beginning of simulation
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

    // need to modify !!
    int atm_nsteps = ntimes + 1;
    const auto fstart = ELM::Utils::Date(1985, 1, 1);
    auto forc_TBOT = create_forc_util<AtmForcType::TBOT>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_PBOT = create_forc_util<AtmForcType::PBOT>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_QBOT = create_forc_util<AtmForcType::RH>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_FLDS = create_forc_util<AtmForcType::FLDS>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_FSDS = create_forc_util<AtmForcType::FSDS>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_PREC = create_forc_util<AtmForcType::PREC>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_WIND = create_forc_util<AtmForcType::WIND>(fname_forc, fstart, atm_nsteps, ncells);
    auto forc_ZBOT = create_forc_util<AtmForcType::ZBOT>(fname_forc, fstart, atm_nsteps, ncells);
    
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
    ELM::SnicarData<ViewD1, ViewD2, ViewD3> snicar_data;
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

    // snow aging parameters
    // this is a large lookup table
    // need to figure out appropriate method to handle 
    // this table and Mie aerosol model lookup table
    ELM::SnwRdsTable<ViewD3> snw_rds_table;
    {
      auto host_snowage_d3 = get_snowage_host_views_d3(snw_rds_table);
      ELM::read_snowrds_data(host_snowage_d3, dd.comm, fname_snowage);
      copy_snowage_host_views_d3(host_snowage_d3, snw_rds_table);
    }



    // pft data constants
    ELM::PFTData<ViewD1, ViewD2> pft_data;
    {
      auto host_pft_views = get_pft_host_views(pft_data);
      ELM::read_pft_data(host_pft_views, dd.comm, fname_pft);
      copy_pft_host_views(host_pft_views, pft_data);
    }
    // Kokkos view of struct PSNVegData
    auto psn_pft = create<Kokkos::View<ELM::PFTDataPSN *>>("psn_pft", ncells);
    Kokkos::parallel_for("pft_psn_init", ncells, KOKKOS_LAMBDA (const int i) {
      psn_pft(i) = pft_data.get_pft_psn(vtype(i));
    });


    // aerosol deposition data manager
    ELM::AerosolDataManager<ViewD1> aerosol_data;
    {
      auto host_aero_views = get_aero_host_views(aerosol_data);
      ELM::read_aerosol_data(host_aero_views, dd.comm, fname_aerosol, lon, lat);
      copy_aero_host_views(host_aero_views, aerosol_data);
    }

    // phenology data manager
    // make host mirrors - need to be persistent
    ELM::PhenologyDataManager<ViewD2> phen_data(dd, ncells, 17);
    auto host_phen_views = get_phen_host_views(phen_data);

    // containers for aerosol deposition and concentration within snowpack layers
    ELM::AerosolMasses<ViewD2> aerosol_masses(ncells);
    ELM::AerosolConcentrations<ViewD2> aerosol_concentrations(ncells);




    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // some free functions that calculate initial values
    // need to be run once for each cell
    // can be parallel
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

    Kokkos::parallel_for("init_functions", ncells, KOKKOS_LAMBDA (const int idx) {
    
      ELM::init_topo_slope(topo_slope(idx));

      ELM::init_micro_topo(
                          Land.ltype, topo_slope(idx),
                          topo_std(idx), n_melt(idx),
                          micro_sigma(idx));

      ELM::init_snow_layers(
                            snow_depth(idx), lakpoi, snl(idx),
                            Kokkos::subview(dz, idx, Kokkos::ALL),
                            Kokkos::subview(zsoi, idx, Kokkos::ALL),
                            Kokkos::subview(zisoi, idx, Kokkos::ALL));

      ELM::init_soil_hydraulics(
                              Kokkos::subview(pct_sand, idx, Kokkos::ALL),
                              Kokkos::subview(pct_clay, idx, Kokkos::ALL),
                              Kokkos::subview(organic, idx, Kokkos::ALL),
                              Kokkos::subview(zsoi, idx, Kokkos::ALL),
                              Kokkos::subview(watsat, idx, Kokkos::ALL),
                              Kokkos::subview(bsw, idx, Kokkos::ALL),
                              Kokkos::subview(sucsat, idx, Kokkos::ALL),
                              Kokkos::subview(watdry, idx, Kokkos::ALL),
                              Kokkos::subview(watopt, idx, Kokkos::ALL),
                              Kokkos::subview(watfc, idx, Kokkos::ALL),
                              Kokkos::subview(tkmg, idx, Kokkos::ALL),
                              Kokkos::subview(tkdry, idx, Kokkos::ALL),
                              Kokkos::subview(csol, idx, Kokkos::ALL));


      ELM::init_vegrootfr(
                          vtype(idx), pft_data.roota_par(vtype(idx)),
                          pft_data.rootb_par(vtype(idx)),
                          Kokkos::subview(zisoi, idx, Kokkos::ALL),
                          Kokkos::subview(rootfr, idx, Kokkos::ALL));

      ELM::init_soil_temp(
                          Land, snl(idx),
                          Kokkos::subview(t_soisno, idx, Kokkos::ALL),
                          t_grnd(idx));

      ELM::init_snow_state(
                          Land.urbpoi, snl(idx), h2osno(idx), int_snow(idx),
                          snow_depth(idx), h2osfc(idx), h2ocan(idx), frac_h2osfc(idx),
                          fwet(idx), fdry(idx), frac_sno(idx),
                          Kokkos::subview(snw_rds, idx, Kokkos::ALL));

      ELM::init_soilh2o_state(
                              Land, snl(idx),
                              Kokkos::subview(watsat, idx, Kokkos::ALL),
                              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
                              Kokkos::subview(dz, idx, Kokkos::ALL),
                              Kokkos::subview(h2osoi_vol, idx, Kokkos::ALL),
                              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
                              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL));
    });


    // hardwired soil/snow water state
    {
      double h2osoi_ice_hardwire[] = {
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 51.095355179469955, 131.99213225849098,
        17.829256395227745, 95.72899575304584, 155.31526899797177,
        0.01, 0.01, 0.01,
        0.01, 0.01 };

      double h2osoi_liq_hardwire[] = {
        0.0, 0.0, 0.0,
        0.0, 0.0, 7.045411435071487,
        14.353496179256807, 36.308518784697064, 62.46145027256513,
        97.14000248023912, 97.47148319510016, 78.52160092062527,
        65.63904088905001, 41.25305599181871, 70.8566046019581,
        0.01, 0.01, 0.01,
        0.01, 0.01 };

      auto h_soi_ice = Kokkos::create_mirror_view(h2osoi_ice);
      auto h_soi_liq = Kokkos::create_mirror_view(h2osoi_liq);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
          h_soi_ice(n, i) = h2osoi_ice_hardwire[i];
          h_soi_liq(n, i) = h2osoi_liq_hardwire[i];
        }
      }
      Kokkos::deep_copy(h2osoi_ice, h_soi_ice);
      Kokkos::deep_copy(h2osoi_liq, h_soi_liq);


      double h2osoi_vol_hardwire[] = {
        0.4016484663460637, 0.5196481455614503, 0.7967166638201649,
        0.8331813710901114, 0.7859200286330449, 0.7517405589446893,
        0.6621235242027332, 0.1535948180493002, 0.15947477948341815,
        0.15954052527228618, 8.420726808634413e-06, 5.107428986500891e-06,
        3.0978122726178113e-06, 1.8789181213767733e-06, 1.5092697845407248e-06 };
      auto h_soi_vol = Kokkos::create_mirror_view(h2osoi_vol);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevgrnd; ++i) {
          h_soi_vol(n, i) = h2osoi_vol_hardwire[i];
        }
      }
      Kokkos::deep_copy(h2osoi_vol, h_soi_vol);

    }


    // hardwired soil/snow temp info
    {
     double tsoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 278.3081064745931,
      276.1568781897738, 275.55803480737063, 275.2677090940866,
      274.7286996980052, 273.15, 272.4187794248787,
      270.65049816473027, 267.8224112387398, 265.7450135695632,
      264.49481140089864, 264.14163363048056, 264.3351872934207,
      264.1163763444719, 263.88852987294865 };
      auto h_tsoi = Kokkos::create_mirror_view(t_soisno);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
          h_tsoi(n, i) = tsoi_hardwire[i];
        }
      }
      auto h_tgrnd = Kokkos::create_mirror_view(t_grnd);
      h_tgrnd(idx) = h_tsoi(idx, nlevsno - snl(idx));
      Kokkos::deep_copy(t_soisno, h_tsoi);
      Kokkos::deep_copy(t_grnd, h_tgrnd);
    }




    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /*                          TIME LOOP                                                                  */
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // these calls should be inside time loop
    // leave here for now - only run kernels (atm_nsteps - 1) times
    read_atm_data(forc_TBOT, dd, start, atm_nsteps);
    read_atm_data(forc_PBOT, dd, start, atm_nsteps);
    read_atm_data(forc_QBOT, dd, start, atm_nsteps);
    read_atm_data(forc_FLDS, dd, start, atm_nsteps);
    read_atm_data(forc_FSDS, dd, start, atm_nsteps);
    read_atm_data(forc_PREC, dd, start, atm_nsteps);
    read_atm_data(forc_WIND, dd, start, atm_nsteps);
    read_atm_data(forc_ZBOT, dd, start, atm_nsteps);

    auto coszen = create<ViewD1>("coszen", ncells);

    ELM::Utils::Date current(start);

    for (int t = 0; t < ntimes; ++t) {

      ELM::Utils::Date time_plus_half_dt(current);
      time_plus_half_dt.increment_seconds(dtime/2);

      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // get coszen
      // only one value currently
      // will change when slope aspect modifier is completed
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      double max_dayl;
      double dayl;
      {
        // there are three methods to calculate zenith angle
        // they all produce similar results for the lat/lon tested here.

        // for now a single value for coszen is appropriate
        // but a slope based factor would necessitate per-cell values

        // first method - average cosz for dt_start to dt_end
        auto decday = ELM::Utils::decimal_doy(current) + 1.0;
        assign(coszen, ELM::incident_shortwave::average_cosz(lat_r, lon_r, dtime, decday));

        // second method - point cosz at dt_start + dt/2
        //auto thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, decday + dtime / 86400.0 /2.0);

        // third method - calc avg cosz over forcing dt (larger than model dt)
        // then calc point dt at start + dt/2
        // and use to calculate cosz_factor
        //ELM::Utils::Date forc_dt_start{forc_FSDS.get_data_start_time()};
        //forc_dt_start.increment_seconds(round(forc_FSDS.forc_t_idx(time_plus_half_dt, forc_FSDS.get_data_start_time()) * forc_FSDS.get_forc_dt_secs()));
        //double cosz_forc_decday = ELM::Utils::decimal_doy(forc_dt_start) + 1.0;
        //auto cosz_forcdt_avg = ELM::incident_shortwave::average_cosz(lat_r, lon_r, forc_FSDS.get_forc_dt_secs(), cosz_forc_decday);
        //auto thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, decday + dtime / 86400.0 /2.0);
        //cosz_factor = (thiscosz > 0.001) ? std::min(thiscosz/cosz_forcdt_avg, 10.0) : 0.0;



        max_dayl = ELM::max_daylength(lat_r);
        dayl = ELM::daylength(lat_r, ELM::incident_shortwave::declination_angle2(current.doy + 1));
      }


      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // timestep init functions
      // read time-variable data
      // these are all self-invoking parallel kernels
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
      // will fix later - too infrequently run (once per month) to cause concern
      {
        if (phen_data.need_data()) {
          Kokkos::deep_copy(host_phen_views["MONTHLY_LAI"], phen_data.mlai);
          Kokkos::deep_copy(host_phen_views["MONTHLY_SAI"], phen_data.msai);
          Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_TOP"], phen_data.mhtop);
          Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_BOT"], phen_data.mhbot);
        }
        // reads three months of data on first call
        // after first call, read new data if phen_data.need_new_data_ == true
        auto phen_updated = phen_data.read_data(host_phen_views, fname_surfdata, current, vtype); // if needed
        // copy host views to device
        // could be made more efficient, see above
        if (phen_updated) {
          Kokkos::deep_copy(phen_data.mlai, host_phen_views["MONTHLY_LAI"]);
          Kokkos::deep_copy(phen_data.msai, host_phen_views["MONTHLY_SAI"]);
          Kokkos::deep_copy(phen_data.mhtop, host_phen_views["MONTHLY_HEIGHT_TOP"]);
          Kokkos::deep_copy(phen_data.mhbot, host_phen_views["MONTHLY_HEIGHT_BOT"]);
        }
        // run parallel kernel to process phenology data
        phen_data.get_data(current, snow_depth,
                           frac_sno, vtype, elai, esai,
                           htop, hbot, tlai, tsai,
                           frac_veg_nosno_alb);
      }


      // get forcing data and process in parallel
      {
        forc_TBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_thbot);
        forc_PBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot);
        forc_QBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_pbot, forc_qbot, forc_rh);
        forc_FLDS.get_atm_forcing(dtime_d, time_plus_half_dt, forc_pbot, forc_qbot, forc_tbot, forc_lwrad);
        forc_FSDS.get_atm_forcing(dtime_d, time_plus_half_dt, coszen, forc_solai, forc_solad);
        forc_PREC.get_atm_forcing(dtime_d, time_plus_half_dt, forc_tbot, forc_rain, forc_snow);
        forc_WIND.get_atm_forcing(dtime_d, time_plus_half_dt, forc_u, forc_v);
        forc_ZBOT.get_atm_forcing(dtime_d, time_plus_half_dt, forc_hgt, forc_hgt_u, forc_hgt_t,  forc_hgt_q);
      }

      // calculate constitutive air properties
      {
        ELM::atm_forcing_physics::ConstitutiveAirProperties
          compute_air(forc_qbot, forc_pbot,
                      forc_tbot, forc_vp,
                      forc_rho, forc_po2,
                      forc_pco2);

        invoke_kernel(compute_air, std::make_tuple(forc_pbot.extent(0)), "ConstitutiveAirProperties");
      }


      // get aerosol mss and cnc
      {
        ELM::aerosols::invoke_aerosol_source(time_plus_half_dt, dtime, snl, aerosol_data, aerosol_masses);
        ELM::aerosols::invoke_aerosol_concen_and_mass(dtime, do_capsnow, snl, h2osoi_liq,
        h2osoi_ice, snw_rds, qflx_snwcp_ice, aerosol_masses, aerosol_concentrations);
      }



      {
        Kokkos::parallel_for("init_spatial_loop", ncells, KOKKOS_LAMBDA (const int idx) {

          ELM::init_timestep(lakpoi, veg_active(idx),
                             frac_veg_nosno_alb(idx),
                             snl(idx), h2osno(idx),
                             Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
                             Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
                             do_capsnow(idx),
                             frac_veg_nosno(idx),
                             Kokkos::subview(frac_iceold, idx, Kokkos::ALL));
        });
      }


      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // Main physics calls
      // single large loop to call all physics kernels
      // will move to less naive approach soon
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

      Kokkos::parallel_for("first_spatial_loop", ncells, KOKKOS_LAMBDA (const int idx) {


        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call surface albedo and SNICAR kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          // parse pft data for Land.vtype
          ELM::PFTDataAlb alb_pft = pft_data.get_pft_alb(vtype(idx));

          ELM::surface_albedo::init_timestep(
              Land.urbpoi,
              elai(idx),
              Kokkos::subview(aerosol_concentrations.mss_cnc_bcphi, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_concentrations.mss_cnc_bcpho, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_concentrations.mss_cnc_dst1, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_concentrations.mss_cnc_dst2, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_concentrations.mss_cnc_dst3, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_concentrations.mss_cnc_dst4, idx, Kokkos::ALL),
              vcmaxcintsun(idx),
              vcmaxcintsha(idx),
              Kokkos::subview(albsod, idx, Kokkos::ALL),
              Kokkos::subview(albsoi, idx, Kokkos::ALL),
              Kokkos::subview(albgrd, idx, Kokkos::ALL),
              Kokkos::subview(albgri, idx, Kokkos::ALL),
              Kokkos::subview(albd, idx, Kokkos::ALL),
              Kokkos::subview(albi, idx, Kokkos::ALL),
              Kokkos::subview(fabd, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sun, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sha, idx, Kokkos::ALL),
              Kokkos::subview(fabi, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sun, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sha, idx, Kokkos::ALL),
              Kokkos::subview(ftdd, idx, Kokkos::ALL),
              Kokkos::subview(ftid, idx, Kokkos::ALL),
              Kokkos::subview(ftii, idx, Kokkos::ALL),
              Kokkos::subview(flx_absdv, idx, Kokkos::ALL),
              Kokkos::subview(flx_absdn, idx, Kokkos::ALL),
              Kokkos::subview(flx_absiv, idx, Kokkos::ALL),
              Kokkos::subview(flx_absin, idx, Kokkos::ALL),
              Kokkos::subview(mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL));

          ELM::surface_albedo::soil_albedo(
              Land,
              snl(idx),
              t_grnd(idx),
              coszen(idx),
              Kokkos::subview(h2osoi_vol, idx, Kokkos::ALL),
              Kokkos::subview(albsat, isoicol(idx), Kokkos::ALL),
              Kokkos::subview(albdry, isoicol(idx), Kokkos::ALL),
              Kokkos::subview(albsod, idx, Kokkos::ALL),
              Kokkos::subview(albsoi, idx, Kokkos::ALL));

          {
            int flg_slr_in = 1; // direct-beam

            ELM::snow_snicar::init_timestep (
                Land.urbpoi,
                flg_slr_in,
                coszen(idx),
                h2osno(idx),
                snl(idx),
                Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
                Kokkos::subview(snw_rds, idx, Kokkos::ALL),
                snl_top(idx),
                snl_btm(idx),
                Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL),
                flg_nosnl(idx),
                Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
                Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
                mu_not(idx),
                Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL));

            ELM::snow_snicar::snow_aerosol_mie_params(
                Land.urbpoi,
                flg_slr_in,
                snl_top(idx),
                snl_btm(idx),
                coszen(idx),
                h2osno(idx),
                Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
                snicar_data.ss_alb_oc1,
                snicar_data.asm_prm_oc1,
                snicar_data.ext_cff_mss_oc1,
                snicar_data.ss_alb_oc2,
                snicar_data.asm_prm_oc2,
                snicar_data.ext_cff_mss_oc2,
                snicar_data.ss_alb_dst1,
                snicar_data.asm_prm_dst1,
                snicar_data.ext_cff_mss_dst1,
                snicar_data.ss_alb_dst2,
                snicar_data.asm_prm_dst2,
                snicar_data.ext_cff_mss_dst2,
                snicar_data.ss_alb_dst3,
                snicar_data.asm_prm_dst3,
                snicar_data.ext_cff_mss_dst3,
                snicar_data.ss_alb_dst4,
                snicar_data.asm_prm_dst4,
                snicar_data.ext_cff_mss_dst4,
                snicar_data.ss_alb_snw_drc,
                snicar_data.asm_prm_snw_drc,
                snicar_data.ext_cff_mss_snw_drc,
                snicar_data.ss_alb_snw_dfs,
                snicar_data.asm_prm_snw_dfs,
                snicar_data.ext_cff_mss_snw_dfs,
                snicar_data.ss_alb_bc1,
                snicar_data.asm_prm_bc1,
                snicar_data.ext_cff_mss_bc1,
                snicar_data.ss_alb_bc2,
                snicar_data.asm_prm_bc2,
                snicar_data.ext_cff_mss_bc2,
                snicar_data.bcenh,
                Kokkos::subview(mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL));

            ELM::snow_snicar::snow_radiative_transfer_solver(
                Land.urbpoi,
                flg_slr_in,
                flg_nosnl(idx),
                snl_top(idx),
                snl_btm(idx),
                coszen(idx),
                h2osno(idx),
                mu_not(idx),
                Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL),
                Kokkos::subview(albsoi, idx, Kokkos::ALL),
                Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL));    

            ELM::snow_snicar::snow_albedo_radiation_factor(
                Land.urbpoi,
                flg_slr_in,
                snl_top(idx),
                coszen(idx),
                mu_not(idx),
                h2osno(idx),
                Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
                Kokkos::subview(albsoi, idx, Kokkos::ALL),
                Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(albsnd, idx, Kokkos::ALL),
                Kokkos::subview(flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL));
          }


          {
            int flg_slr_in = 2; // diffuse

            ELM::snow_snicar::init_timestep (
                Land.urbpoi,
                flg_slr_in,
                coszen(idx),
                h2osno(idx),
                snl(idx),
                Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
                Kokkos::subview(snw_rds, idx, Kokkos::ALL),
                snl_top(idx),
                snl_btm(idx),
                Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL),
                flg_nosnl(idx),
                Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
                Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
                mu_not(idx),
                Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL));

            ELM::snow_snicar::snow_aerosol_mie_params(
                Land.urbpoi,
                flg_slr_in,
                snl_top(idx),
                snl_btm(idx),
                coszen(idx),
                h2osno(idx),
                Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
                Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
                snicar_data.ss_alb_oc1,
                snicar_data.asm_prm_oc1,
                snicar_data.ext_cff_mss_oc1,
                snicar_data.ss_alb_oc2,
                snicar_data.asm_prm_oc2,
                snicar_data.ext_cff_mss_oc2,
                snicar_data.ss_alb_dst1,
                snicar_data.asm_prm_dst1,
                snicar_data.ext_cff_mss_dst1,
                snicar_data.ss_alb_dst2,
                snicar_data.asm_prm_dst2,
                snicar_data.ext_cff_mss_dst2,
                snicar_data.ss_alb_dst3,
                snicar_data.asm_prm_dst3,
                snicar_data.ext_cff_mss_dst3,
                snicar_data.ss_alb_dst4,
                snicar_data.asm_prm_dst4,
                snicar_data.ext_cff_mss_dst4,
                snicar_data.ss_alb_snw_drc,
                snicar_data.asm_prm_snw_drc,
                snicar_data.ext_cff_mss_snw_drc,
                snicar_data.ss_alb_snw_dfs,
                snicar_data.asm_prm_snw_dfs,
                snicar_data.ext_cff_mss_snw_dfs,
                snicar_data.ss_alb_bc1,
                snicar_data.asm_prm_bc1,
                snicar_data.ext_cff_mss_bc1,
                snicar_data.ss_alb_bc2,
                snicar_data.asm_prm_bc2,
                snicar_data.ext_cff_mss_bc2,
                snicar_data.bcenh,
                Kokkos::subview(mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL));

            ELM::snow_snicar::snow_radiative_transfer_solver(
                Land.urbpoi,
                flg_slr_in,
                flg_nosnl(idx),
                snl_top(idx),
                snl_btm(idx),
                coszen(idx),
                h2osno(idx),
                mu_not(idx),
                Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL),
                Kokkos::subview(albsoi, idx, Kokkos::ALL),
                Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL));

            ELM::snow_snicar::snow_albedo_radiation_factor(
                Land.urbpoi,
                flg_slr_in,
                snl_top(idx),
                coszen(idx),
                mu_not(idx),
                h2osno(idx),
                Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
                Kokkos::subview(albsoi, idx, Kokkos::ALL),
                Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
                Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(albsni, idx, Kokkos::ALL),
                Kokkos::subview(flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL));
          }

          ELM::surface_albedo::ground_albedo(
              Land.urbpoi,
              coszen(idx),
              frac_sno(idx),
              Kokkos::subview(albsod, idx, Kokkos::ALL),
              Kokkos::subview(albsoi, idx, Kokkos::ALL),
              Kokkos::subview(albsnd, idx, Kokkos::ALL),
              Kokkos::subview(albsni, idx, Kokkos::ALL),
              Kokkos::subview(albgrd, idx, Kokkos::ALL),
              Kokkos::subview(albgri, idx, Kokkos::ALL));

          ELM::surface_albedo::flux_absorption_factor(
              Land,
              coszen(idx),
              frac_sno(idx),
              Kokkos::subview(albsod, idx, Kokkos::ALL),
              Kokkos::subview(albsoi, idx, Kokkos::ALL),
              Kokkos::subview(albsnd, idx, Kokkos::ALL),
              Kokkos::subview(albsni, idx, Kokkos::ALL),
              Kokkos::subview(flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL),
              Kokkos::subview(flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL),
              Kokkos::subview(flx_absdv, idx, Kokkos::ALL),
              Kokkos::subview(flx_absdn, idx, Kokkos::ALL),
              Kokkos::subview(flx_absiv, idx, Kokkos::ALL),
              Kokkos::subview(flx_absin, idx, Kokkos::ALL));

          ELM::surface_albedo::canopy_layer_lai(
              Land.urbpoi,
              elai(idx),
              esai(idx),
              tlai(idx),
              tsai(idx),
              nrad(idx),
              ncan(idx),
              Kokkos::subview(tlai_z, idx, Kokkos::ALL),
              Kokkos::subview(tsai_z, idx, Kokkos::ALL),
              Kokkos::subview(fsun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sha_z, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sha_z, idx, Kokkos::ALL));

          ELM::surface_albedo::two_stream_solver(
              Land,
              nrad(idx),
              coszen(idx),
              t_veg(idx),
              fwet(idx),
              elai(idx),
              esai(idx),
              Kokkos::subview(tlai_z, idx, Kokkos::ALL),
              Kokkos::subview(tsai_z, idx, Kokkos::ALL),
              Kokkos::subview(albgrd, idx, Kokkos::ALL),
              Kokkos::subview(albgri, idx, Kokkos::ALL),
              alb_pft,
              vcmaxcintsun(idx),
              vcmaxcintsha(idx),
              Kokkos::subview(albd, idx, Kokkos::ALL),
              Kokkos::subview(ftid, idx, Kokkos::ALL),
              Kokkos::subview(ftdd, idx, Kokkos::ALL),
              Kokkos::subview(fabd, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sun, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sha, idx, Kokkos::ALL),
              Kokkos::subview(albi, idx, Kokkos::ALL),
              Kokkos::subview(ftii, idx, Kokkos::ALL),
              Kokkos::subview(fabi, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sun, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sha, idx, Kokkos::ALL),
              Kokkos::subview(fsun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sha_z, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sha_z, idx, Kokkos::ALL));
        }

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call canopy_hydrology kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          // local vars - these need to be thread local in parallel runs
          double qflx_candrip;
          double qflx_through_snow;
          double qflx_through_rain;
          double fracsnow;
          double fracrain;

          double qflx_irrig = 0.0; // hardwired here

          ELM::canopy_hydrology::interception(
              Land,
              frac_veg_nosno(idx),
              forc_rain(idx),
              forc_snow(idx),
              dewmx,
              elai(idx),
              esai(idx),
              dtime,
              h2ocan(idx),
              qflx_candrip,
              qflx_through_snow,
              qflx_through_rain,
              fracsnow,
              fracrain);

          ELM::canopy_hydrology::ground_flux(
              Land,
              do_capsnow(idx),
              frac_veg_nosno(idx),
              forc_rain(idx),
              forc_snow(idx),
              qflx_irrig,
              qflx_candrip,
              qflx_through_snow,
              qflx_through_rain,
              fracsnow,
              fracrain,
              qflx_prec_grnd(idx),
              qflx_snwcp_liq(idx),
              qflx_snwcp_ice(idx),
              qflx_snow_grnd(idx),
              qflx_rain_grnd(idx));

          ELM::canopy_hydrology::fraction_wet(
              Land,
              frac_veg_nosno(idx),
              dewmx,
              elai(idx),
              esai(idx),
              h2ocan(idx),
              fwet(idx),
              fdry(idx));

          ELM::canopy_hydrology::snow_init(
              Land,
              dtime,
              do_capsnow(idx),
              oldfflag,
              forc_tbot(idx),
              t_grnd(idx),
              qflx_snow_grnd(idx),
              qflx_snow_melt(idx),
              n_melt(idx),
              snow_depth(idx),
              h2osno(idx),
              int_snow(idx),
              Kokkos::subview(swe_old, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(frac_iceold, idx, Kokkos::ALL),
              snl(idx),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(zsoi, idx, Kokkos::ALL),
              Kokkos::subview(zisoi, idx, Kokkos::ALL),
              Kokkos::subview(snw_rds, idx, Kokkos::ALL),
              frac_sno_eff(idx),
              frac_sno(idx));

          ELM::canopy_hydrology::fraction_h2osfc(
              Land,
              micro_sigma(idx),
              h2osno(idx),
              h2osfc(idx),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              frac_sno(idx),
              frac_sno_eff(idx),
              frac_h2osfc(idx));
        }

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call surface_radiation kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          // local to these kernel calls
            double trd[numrad] = {0.0,0.0};
            double tri[numrad] = {0.0,0.0};

          // call canopy_sunshade_fractions kernel
          ELM::surface_radiation::canopy_sunshade_fractions(
              Land,
              nrad(idx),
              elai(idx),
              Kokkos::subview(tlai_z, idx, Kokkos::ALL),
              Kokkos::subview(fsun_z, idx, Kokkos::ALL),
              Kokkos::subview(forc_solad, idx, Kokkos::ALL),
              Kokkos::subview(forc_solai, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabd_sha_z, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sun_z, idx, Kokkos::ALL),
              Kokkos::subview(fabi_sha_z, idx, Kokkos::ALL),
              Kokkos::subview(parsun_z, idx, Kokkos::ALL),
              Kokkos::subview(parsha_z, idx, Kokkos::ALL),
              Kokkos::subview(laisun_z, idx, Kokkos::ALL),
              Kokkos::subview(laisha_z, idx, Kokkos::ALL),
              laisun(idx),
              laisha(idx));

          ELM::surface_radiation::initialize_flux(
              Land,
              sabg_soil(idx),
              sabg_snow(idx),
              sabg(idx),
              sabv(idx),
              fsa(idx),
              Kokkos::subview(sabg_lyr, idx, Kokkos::ALL));

          ELM::surface_radiation::total_absorbed_radiation(
              Land,
              snl(idx),
              Kokkos::subview(ftdd, idx, Kokkos::ALL),
              Kokkos::subview(ftid, idx, Kokkos::ALL),
              Kokkos::subview(ftii, idx, Kokkos::ALL),
              Kokkos::subview(forc_solad, idx, Kokkos::ALL),
              Kokkos::subview(forc_solai, idx, Kokkos::ALL),
              Kokkos::subview(fabd, idx, Kokkos::ALL),
              Kokkos::subview(fabi, idx, Kokkos::ALL),
              Kokkos::subview(albsod, idx, Kokkos::ALL),
              Kokkos::subview(albsoi, idx, Kokkos::ALL),
              Kokkos::subview(albsnd_hst, idx, Kokkos::ALL),
              Kokkos::subview(albsni_hst, idx, Kokkos::ALL),
              Kokkos::subview(albgrd, idx, Kokkos::ALL),
              Kokkos::subview(albgri, idx, Kokkos::ALL),
              sabv(idx),
              fsa(idx),
              sabg(idx),
              sabg_soil(idx),
              sabg_snow(idx),
              trd,
              tri);

          ELM::surface_radiation::layer_absorbed_radiation(
              Land,
              snl(idx),
              sabg(idx),
              sabg_snow(idx),
              snow_depth(idx),
              Kokkos::subview(flx_absdv, idx, Kokkos::ALL),
              Kokkos::subview(flx_absdn, idx, Kokkos::ALL),
              Kokkos::subview(flx_absiv, idx, Kokkos::ALL),
              Kokkos::subview(flx_absin, idx, Kokkos::ALL),
              trd,
              tri,
              Kokkos::subview(sabg_lyr, idx, Kokkos::ALL));

          ELM::surface_radiation::reflected_radiation(
              Land,
              Kokkos::subview(albd, idx, Kokkos::ALL),
              Kokkos::subview(albi, idx, Kokkos::ALL),
              Kokkos::subview(forc_solad, idx, Kokkos::ALL),
              Kokkos::subview(forc_solai, idx, Kokkos::ALL),
              fsr(idx));
        }

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call canopy_temperature kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          double qred; // soil surface relative humidity
          double hr;   // relative humidity
          ELM::canopy_temperature::old_ground_temp(
              Land,
              t_h2osfc(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              t_h2osfc_bef(idx),
              Kokkos::subview(tssbef, idx, Kokkos::ALL));

          ELM::canopy_temperature::ground_temp(
              Land,
              snl(idx),
              frac_sno_eff(idx),
              frac_h2osfc(idx),
              t_h2osfc(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              t_grnd(idx));

          ELM::canopy_temperature::calc_soilalpha(
              Land,
              frac_sno(idx),
              frac_h2osfc(idx),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(watsat, idx, Kokkos::ALL),
              Kokkos::subview(sucsat, idx, Kokkos::ALL),
              Kokkos::subview(bsw, idx, Kokkos::ALL),
              Kokkos::subview(watdry, idx, Kokkos::ALL),
              Kokkos::subview(watopt, idx, Kokkos::ALL),
              Kokkos::subview(rootfr_road_perv, idx, Kokkos::ALL),
              Kokkos::subview(rootr_road_perv, idx, Kokkos::ALL),
              qred, hr,
              soilalpha(idx),
              soilalpha_u(idx));

          ELM::canopy_temperature::calc_soilbeta(
              Land,
              frac_sno(idx),
              frac_h2osfc(idx),
              Kokkos::subview(watsat, idx, Kokkos::ALL),
              Kokkos::subview(watfc, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              soilbeta(idx));

          ELM::canopy_temperature::humidities(
              Land,
              snl(idx),
              forc_qbot(idx),
              forc_pbot(idx),
              t_h2osfc(idx),
              t_grnd(idx),
              frac_sno(idx),
              frac_sno_eff(idx),
              frac_h2osfc(idx),
              qred,
              hr,
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              qg_snow(idx),
              qg_soil(idx),
              qg(idx),
              qg_h2osfc(idx),
              dqgdT(idx));

          ELM::canopy_temperature::ground_properties(
              Land,
              snl(idx),
              frac_sno(idx),
              forc_thbot(idx),
              forc_qbot(idx),
              elai(idx),
              esai(idx),
              htop(idx),
              pft_data.displar,
              pft_data.z0mr,
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              emg(idx),
              emv(idx),
              htvp(idx),
              z0mg(idx),
              z0hg(idx),
              z0qg(idx),
              z0mv(idx),
              z0hv(idx),
              z0qv(idx),
              thv(idx),
              z0m(idx),
              displa(idx));

          ELM::canopy_temperature::forcing_height(
              Land,
              veg_active(idx),
              frac_veg_nosno(idx),
              forc_hgt_u(idx),
              forc_hgt_t(idx),
              forc_hgt_q(idx),
              z0m(idx),
              z0mg(idx),
              z_0_town(idx),
              z_d_town(idx),
              forc_tbot(idx),
              displa(idx),
              forc_hgt_u_patch(idx),
              forc_hgt_t_patch(idx),
              forc_hgt_q_patch(idx),
              thm(idx));

          ELM::canopy_temperature::init_energy_fluxes(
              Land,
              eflx_sh_tot(idx),
              eflx_sh_tot_u(idx),
              eflx_sh_tot_r(idx),
              eflx_lh_tot(idx),
              eflx_lh_tot_u(idx),
              eflx_lh_tot_r(idx),
              eflx_sh_veg(idx),
              qflx_evap_tot(idx),
              qflx_evap_veg(idx),
              qflx_tran_veg(idx));
        }

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call bareground_fluxes kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          // temporary data to pass between functions
          double zldis;   // reference height "minus" zero displacement height [m]
          double displa;  // displacement height [m]
          double dth;     // diff of virtual temp. between ref. height and surface
          double dqh;     // diff of humidity between ref. height and surface
          double obu;     // Monin-Obukhov length (m)
          double ur;      // wind speed at reference height [m/s]
          double um;      // wind speed including the stablity effect [m/s]
          double temp1;   // relation for potential temperature profile
          double temp2;   // relation for specific humidity profile
          double temp12m; // relation for potential temperature profile applied at 2-m
          double temp22m; // relation for specific humidity profile applied at 2-m
          double ustar;   // friction velocity [m/s]

          ELM::bareground_fluxes::initialize_flux(
              Land,
              frac_veg_nosno(idx),
              forc_u(idx),
              forc_v(idx),
              forc_qbot(idx),
              forc_thbot(idx),
              forc_hgt_u_patch(idx),
              thm(idx),
              thv(idx),
              t_grnd(idx),
              qg(idx),
              z0mg(idx),
              dlrad(idx),
              ulrad(idx),
              zldis,
              displa,
              dth,
              dqh,
              obu,
              ur,
              um);

          ELM::bareground_fluxes::stability_iteration(
              Land,
              frac_veg_nosno(idx),
              forc_hgt_t_patch(idx),
              forc_hgt_u_patch(idx),
              forc_hgt_q_patch(idx),
              z0mg(idx),
              zldis,
              displa,
              dth,
              dqh,
              ur,
              forc_qbot(idx),
              forc_thbot(idx),
              thv(idx),
              z0hg(idx),
              z0qg(idx),
              obu,
              um,
              temp1,
              temp2,
              temp12m,
              temp22m,
              ustar);

          ELM::bareground_fluxes::compute_flux(
              Land,
              frac_veg_nosno(idx),
              snl(idx),
              forc_rho(idx),
              soilbeta(idx),
              dqgdT(idx),
              htvp(idx),
              t_h2osfc(idx),
              qg_snow(idx),
              qg_soil(idx),
              qg_h2osfc(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              forc_pbot(idx),
              dth,
              dqh,
              temp1,
              temp2,
              temp12m,
              temp22m,
              ustar,
              forc_qbot(idx),
              thm(idx),
              cgrnds(idx),
              cgrndl(idx),
              cgrnd(idx),
              eflx_sh_grnd(idx),
              eflx_sh_tot(idx),
              eflx_sh_snow(idx),
              eflx_sh_soil(idx),
              eflx_sh_h2osfc(idx),
              qflx_evap_soi(idx),
              qflx_evap_tot(idx),
              qflx_ev_snow(idx),
              qflx_ev_soil(idx),
              qflx_ev_h2osfc(idx),
              t_ref2m(idx),
              t_ref2m_r(idx),
              q_ref2m(idx),
              rh_ref2m(idx),
              rh_ref2m_r(idx));
        }

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call canopy_fluxes kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          // temporary data to pass between functions
          double wtg = 0.0;         // heat conductance for ground [m/s]
          double wtgq = 0.0;        // latent heat conductance for ground [m/s]
          double wtalq = 0.0;       // normalized latent heat cond. for air and leaf [-]
          double wtlq0 = 0.0;       // normalized latent heat conductance for leaf [-]
          double wtaq0 = 0.0;       // normalized latent heat conductance for air [-]
          double wtl0 = 0.0;        // normalized heat conductance for leaf [-]
          double wta0 = 0.0;        // normalized heat conductance for air [-]
          double wtal = 0.0;        // normalized heat conductance for air and leaf [-]
          double dayl_factor = 0.0; // scalar (0-1) for daylength effect on Vcmax
          double air = 0.0;         // atmos. radiation temporay set
          double bir = 0.0;         // atmos. radiation temporay set
          double cir = 0.0;         // atmos. radiation temporay set
          double el = 0.0;          // vapor pressure on leaf surface [pa]
          double qsatl = 0.0;       // leaf specific humidity [kg/kg]
          double qsatldT = 0.0;     // derivative of "qsatl" on "t_veg"
          double taf = 0.0;         // air temperature within canopy space [K]
          double qaf = 0.0;         // humidity of canopy air [kg/kg]
          double um = 0.0;          // wind speed including the stablity effect [m/s]
          double ur = 0.0;          // wind speed at reference height [m/s]
          double dth = 0.0;         // diff of virtual temp. between ref. height and surface
          double dqh = 0.0;         // diff of humidity between ref. height and surface
          double obu = 0.0;         // Monin-Obukhov length (m)
          double zldis = 0.0;       // reference height "minus" zero displacement height [m]
          double temp1 = 0.0;       // relation for potential temperature profile
          double temp2 = 0.0;       // relation for specific humidity profile
          double temp12m = 0.0;     // relation for potential temperature profile applied at 2-m
          double temp22m = 0.0;     // relation for specific humidity profile applied at 2-m
          double tlbef = 0.0;       // leaf temperature from previous iteration [K]
          double delq = 0.0;        // temporary
          double dt_veg = 0.0;      // change in t_veg, last iteration (Kelvin)

          ELM::canopy_fluxes::initialize_flux(
              Land,
              snl(idx),
              frac_veg_nosno(idx),
              frac_sno(idx),
              forc_hgt_u_patch(idx),
              thm(idx),
              thv(idx),
              max_dayl,
              dayl,
              altmax_indx(idx),
              altmax_lastyear_indx(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(rootfr, idx, Kokkos::ALL),
              psn_pft(idx).tc_stress,
              Kokkos::subview(sucsat, idx, Kokkos::ALL),
              Kokkos::subview(watsat, idx, Kokkos::ALL),
              Kokkos::subview(bsw, idx, Kokkos::ALL),
              psn_pft(idx).smpso,
              psn_pft(idx).smpsc,
              elai(idx),
              esai(idx),
              emv(idx),
              emg(idx),
              qg(idx),
              t_grnd(idx),
              forc_tbot(idx),
              forc_pbot(idx),
              forc_lwrad(idx),
              forc_u(idx),
              forc_v(idx),
              forc_qbot(idx),
              forc_thbot(idx),
              z0mg(idx),
              btran(idx),
              displa(idx),
              z0mv(idx),
              z0hv(idx),
              z0qv(idx),
              Kokkos::subview(rootr, idx, Kokkos::ALL),
              Kokkos::subview(eff_porosity, idx, Kokkos::ALL),
              dayl_factor,
              air,
              bir,
              cir,
              el,
              qsatl,
              qsatldT,
              taf,
              qaf,
              um,
              ur,
              obu,
              zldis,
              delq,
              t_veg(idx));

          ELM::canopy_fluxes::stability_iteration(
              Land,
              dtime,
              snl(idx),
              frac_veg_nosno(idx),
              frac_sno(idx),
              forc_hgt_u_patch(idx),
              forc_hgt_t_patch(idx),
              forc_hgt_q_patch(idx),
              fwet(idx),
              fdry(idx),
              laisun(idx),
              laisha(idx),
              forc_rho(idx),
              snow_depth(idx),
              soilbeta(idx),
              frac_h2osfc(idx),
              t_h2osfc(idx),
              sabv(idx),
              h2ocan(idx),
              htop(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              air,
              bir,
              cir,
              ur,
              zldis,
              displa(idx),
              elai(idx),
              esai(idx),
              t_grnd(idx),
              forc_pbot(idx),
              forc_qbot(idx),
              forc_thbot(idx),
              z0mg(idx),
              z0mv(idx),
              z0hv(idx),
              z0qv(idx),
              thm(idx),
              thv(idx),
              qg(idx),
              psn_pft(idx),
              nrad(idx),
              t10(idx),
              Kokkos::subview(tlai_z, idx, Kokkos::ALL),
              vcmaxcintsha(idx),
              vcmaxcintsun(idx),
              Kokkos::subview(parsha_z, idx, Kokkos::ALL),
              Kokkos::subview(parsun_z, idx, Kokkos::ALL),
              Kokkos::subview(laisha_z, idx, Kokkos::ALL),
              Kokkos::subview(laisun_z, idx, Kokkos::ALL),
              forc_pco2(idx),
              forc_po2(idx),
              dayl_factor,
              btran(idx),
              qflx_tran_veg(idx),
              qflx_evap_veg(idx),
              eflx_sh_veg(idx),
              wtg,
              wtl0,
              wta0,
              wtal,
              el,
              qsatl,
              qsatldT,
              taf,
              qaf,
              um,
              dth,
              dqh,
              obu,
              temp1,
              temp2,
              temp12m,
              temp22m,
              tlbef,
              delq,
              dt_veg,
              t_veg(idx),
              wtgq,
              wtalq,
              wtlq0,
              wtaq0);

          ELM::canopy_fluxes::compute_flux(
              Land,
              dtime,
              snl(idx),
              frac_veg_nosno(idx),
              frac_sno(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              frac_h2osfc(idx),
              t_h2osfc(idx),
              sabv(idx),
              qg_snow(idx),
              qg_soil(idx),
              qg_h2osfc(idx),
              dqgdT(idx),
              htvp(idx),
              wtg,
              wtl0,
              wta0,
              wtal,
              air,
              bir,
              cir,
              qsatl,
              qsatldT,
              dth,
              dqh,
              temp1,
              temp2,
              temp12m,
              temp22m,
              tlbef,
              delq,
              dt_veg,
              t_veg(idx),
              t_grnd(idx),
              forc_pbot(idx),
              qflx_tran_veg(idx),
              qflx_evap_veg(idx),
              eflx_sh_veg(idx),
              forc_qbot(idx),
              forc_rho(idx),
              thm(idx),
              emv(idx),
              emg(idx),
              forc_lwrad(idx),
              wtgq,
              wtalq,
              wtlq0,
              wtaq0,
              h2ocan(idx),
              eflx_sh_grnd(idx),
              eflx_sh_snow(idx),
              eflx_sh_soil(idx),
              eflx_sh_h2osfc(idx),
              qflx_evap_soi(idx),
              qflx_ev_snow(idx),
              qflx_ev_soil(idx),
              qflx_ev_h2osfc(idx),
              dlrad(idx),
              ulrad(idx),
              cgrnds(idx),
              cgrndl(idx),
              cgrnd(idx),
              t_ref2m(idx),
              t_ref2m_r(idx),
              q_ref2m(idx),
              rh_ref2m(idx),
              rh_ref2m_r(idx));
        }

      }); // parallel for over cells
      


      ELM::soil_temp::solve_temperature<ViewD3>(
          dtime,
          snl,
          frac_veg_nosno,
          dlrad,
          emg,
          forc_lwrad,
          htvp,
          cgrnd,
          eflx_sh_soil,
          qflx_ev_soil,
          eflx_sh_h2osfc,
          qflx_ev_h2osfc,
          eflx_sh_grnd,
          qflx_evap_soi,
          eflx_sh_snow,
          qflx_ev_snow,
          frac_sno_eff,
          frac_sno,
          frac_h2osfc,
          sabg_snow,
          sabg_soil,
          sabg_lyr,
          watsat,
          sucsat,
          bsw,
          tkmg,
          tkdry,
          csol,
          dz,
          zsoi,
          zisoi,
          h2osfc,
          h2osno,
          snow_depth,
          int_snow,
          t_h2osfc,
          t_grnd,
          xmf_h2osfc_dummy,
          xmf_dummy,
          qflx_h2osfc_to_ice_dummy,
          eflx_h2osfc_to_snow_dummy,
          qflx_snofrz,
          qflx_snow_melt,
          qflx_snomelt,
          eflx_snomelt,
          imelt,
          h2osoi_liq,
          h2osoi_ice,
          qflx_snofrz_lyr,
          t_soisno,
          fact);

      Kokkos::parallel_for("second_spatial_loop", ncells, KOKKOS_LAMBDA (const int idx) {

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call snow hydrology kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {

          ELM::snow::snow_water(
              do_capsnow(idx),
              snl(idx),
              dtime,
              frac_sno_eff(idx),
              h2osno(idx),
              qflx_sub_snow(idx),
              qflx_evap_grnd(idx),
              qflx_dew_snow(idx),
              qflx_dew_grnd(idx),
              qflx_rain_grnd(idx),
              qflx_snomelt(idx),
              qflx_snow_melt(idx),
              qflx_top_soil(idx),
              int_snow(idx),
              frac_sno(idx),
              mflx_neg_snow(idx),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst1, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst2, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst3, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst4, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL));
        }

      });

      // aerosol deposition must be called between these two snow hydrology functions
      // look at combining aerosol_phase_change and aerosol deposition
      ELM::aerosols::invoke_aerosol_source(time_plus_half_dt, dtime, snl, aerosol_data, aerosol_masses);

      Kokkos::parallel_for("third_spatial_loop", ncells, KOKKOS_LAMBDA (const int idx) {

        {

          ELM::snow::aerosol_phase_change(
              snl(idx),
              dtime,
              qflx_sub_snow(idx),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL));


          ELM::trans::transpiration(veg_active(idx), qflx_tran_veg(idx),
              Kokkos::subview(rootr, idx, Kokkos::ALL),
              Kokkos::subview(qflx_rootsoi, idx, Kokkos::ALL));

          ELM::snow::snow_compaction(snl(idx),
              ELM::subgridflag,
              Land.ltype,
              dtime,
              int_snow(idx),
              n_melt(idx),
              frac_sno(idx),
              Kokkos::subview(imelt, idx, Kokkos::ALL),
              Kokkos::subview(swe_old, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(frac_iceold, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL));


          ELM::snow::combine_layers(
              Land.urbpoi,
              Land.ltype,
              dtime,
              snl(idx),
              h2osno(idx),
              snow_depth(idx),
              frac_sno_eff(idx),
              frac_sno(idx),
              int_snow(idx),
              qflx_sl_top_soil(idx),
              qflx_snow2topsoi(idx),
              mflx_snowlyr_col(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(snw_rds, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst1, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst2, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst3, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst4, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(zsoi, idx, Kokkos::ALL),
              Kokkos::subview(zisoi, idx, Kokkos::ALL));


          ELM::snow::divide_layers(
              frac_sno(idx),
              snl(idx),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(snw_rds, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcphi, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_bcpho, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst1, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst2, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst3, idx, Kokkos::ALL),
              Kokkos::subview(aerosol_masses.mss_dst4, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(zsoi, idx, Kokkos::ALL),
              Kokkos::subview(zisoi, idx, Kokkos::ALL));


          ELM::snow::prune_snow_layers(
              snl(idx),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(zsoi, idx, Kokkos::ALL),
              Kokkos::subview(zisoi, idx, Kokkos::ALL));


          ELM::snow::snow_aging(
              do_capsnow(idx),
              snl(idx),
              frac_sno(idx),
              dtime,
              qflx_snwcp_ice(idx),
              qflx_snow_grnd(idx),
              h2osno(idx),
              Kokkos::subview(dz, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, idx, Kokkos::ALL),
              Kokkos::subview(h2osoi_ice, idx, Kokkos::ALL),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(qflx_snofrz_lyr, idx, Kokkos::ALL),
              snw_rds_table,
              snot_top(idx),
              dTdz_top(idx),
              snw_rds_top(idx),
              sno_liq_top(idx),
              Kokkos::subview(snw_rds, idx, Kokkos::ALL));

        }


        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        // call surface_fluxes kernels
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
        {
          const auto snotop = nlevsno-snl(idx);
          const auto& soitop = nlevsno;
          ELM::surface_fluxes::initial_flux_calc(
              Land.urbpoi,
              snl(idx),
              frac_sno_eff(idx),
              frac_h2osfc(idx),
              t_h2osfc_bef(idx),
              tssbef(idx, snotop),
              tssbef(idx, soitop),
              t_grnd(idx),
              cgrnds(idx),
              cgrndl(idx),
              eflx_sh_grnd(idx),
              qflx_evap_soi(idx),
              qflx_ev_snow(idx),
              qflx_ev_soil(idx),
              qflx_ev_h2osfc(idx));

          ELM::surface_fluxes::update_surface_fluxes(
              Land.urbpoi,
              do_capsnow(idx),
              snl(idx),
              dtime,
              t_grnd(idx),
              htvp(idx),
              frac_sno_eff(idx),
              frac_h2osfc(idx),
              t_h2osfc_bef(idx),
              sabg_soil(idx),
              sabg_snow(idx),
              dlrad(idx),
              frac_veg_nosno(idx),
              emg(idx),
              forc_lwrad(idx),
              tssbef(idx, snotop),
              tssbef(idx, soitop),
              h2osoi_ice(idx, snotop),
              h2osoi_liq(idx, soitop),
              eflx_sh_veg(idx),
              qflx_evap_veg(idx),
              qflx_evap_soi(idx),
              eflx_sh_grnd(idx),
              qflx_ev_snow(idx),
              qflx_ev_soil(idx),
              qflx_ev_h2osfc(idx),
              eflx_soil_grnd(idx),
              eflx_sh_tot(idx),
              qflx_evap_tot(idx),
              eflx_lh_tot(idx),
              qflx_evap_grnd(idx),
              qflx_sub_snow(idx),
              qflx_dew_snow(idx),
              qflx_dew_grnd(idx),
              qflx_snwcp_liq(idx),
              qflx_snwcp_ice(idx));

          ELM::surface_fluxes::lwrad_outgoing(
              Land.urbpoi,
              snl(idx),
              frac_veg_nosno(idx),
              forc_lwrad(idx),
              frac_sno_eff(idx),
              tssbef(idx, snotop),
              tssbef(idx, soitop),
              frac_h2osfc(idx),
              t_h2osfc_bef(idx),
              t_grnd(idx),
              ulrad(idx),
              emg(idx),
              eflx_lwrad_out(idx),
              eflx_lwrad_net(idx));

          soil_e_balance(idx) = ELM::surface_fluxes::soil_energy_balance(
              Land.ctype,
              snl(idx),
              eflx_soil_grnd(idx),
              xmf_dummy(idx),
              xmf_h2osfc_dummy(idx),
              frac_h2osfc(idx),
              t_h2osfc(idx),
              t_h2osfc_bef(idx),
              dtime,
              eflx_h2osfc_to_snow_dummy(idx),
              eflx_building_heat_dummy(idx),
              frac_sno_eff(idx),
              Kokkos::subview(t_soisno, idx, Kokkos::ALL),
              Kokkos::subview(tssbef, idx, Kokkos::ALL),
              Kokkos::subview(fact, idx, Kokkos::ALL));
        }

      }); // parallel for over cells

  //    std::cout << "soil temp fluxes:  "
  //        << hs_soil(0) << "  "
  //        << hs_h2osfc(0) << "  "
  //        << hs_top(0) << "  "
  //        << hs_top_snow(0) << "  "
  //        << dhsdT(0) << "  "
  //        << soil_e_balance(0) << "  "
  //        << sabg_chk(0) << std::endl;

      std::cout << "lwrad_out:  "
          << eflx_lwrad_out(0) << "  "
          << eflx_lwrad_net(0) << std::endl;

      for (int i = 0; i < nlevsno + nlevgrnd; ++i)
        std::cout << "column vars:  " << i <<
         "  t_soisno:  " << t_soisno(0, i) <<
         "  h2osoi_ice:  " << h2osoi_ice(0, i) <<
         "  h2osoi_liq:  " << h2osoi_liq(0, i) <<
         "  dz:  " << dz(0, i) <<
         "  zsoi:  " << zsoi(0, i) <<
         "  zisoi:  " << zisoi(0, i) << std::endl;


      for (int i = 0; i < ncells; ++i) {
        std::cout << "h2osno: " << h2osno(i) << std::endl;
        std::cout << "t_grnd: " << t_grnd(i) << std::endl;
        std::cout << "snow_depth: " << snow_depth(i) << std::endl;
        std::cout << "frac_sno: " << frac_sno(i) << std::endl;
        std::cout << "frac_sno_eff: " << frac_sno_eff(i) << std::endl;
        std::cout << "frac_veg_nosno_alb: " << frac_veg_nosno_alb(i) << std::endl;
      }

      current.increment_seconds(dtime);

    } // time loop

  } // inner scope

  Kokkos::finalize();
} // enclosing scope
return 0;
}
