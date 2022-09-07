#include "helper_functions.hh"

// conditional compilation options
#include "invoke_kernel.hh"
#include "compile_options.hh"

#include "soil_data.h"
#include "elm_constants.h"

// initialization routines
#include "init_soil_state.h"
#include "init_snow_state.h"
#include "init_topography.h"
#include "soil_texture_hydraulic_model.h"

#include "initialize_elm_kokkos.hh"


namespace ELM::kokkos_init {
std::unordered_map<std::string, h_ViewD1> get_snicar_host_views_d1(const ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  std::unordered_map<std::string, h_ViewD1> snicar_host_views_d1;
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

void copy_snicar_host_views_d1(std::unordered_map<std::string, h_ViewD1>& snicar_host_views_d1, ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
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



std::unordered_map<std::string, h_ViewD2> get_snicar_host_views_d2(const ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  std::unordered_map<std::string, h_ViewD2> snicar_host_views_d2;
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

void copy_snicar_host_views_d2(std::unordered_map<std::string, h_ViewD2>& snicar_host_views_d2,
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



std::unordered_map<std::string, h_ViewD3> get_snicar_host_views_d3(const ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  std::unordered_map<std::string, h_ViewD3> snicar_host_views_d3;
  snicar_host_views_d3["bcint_enh_mam"] = Kokkos::create_mirror_view(snicar_data.bcenh);
  return snicar_host_views_d3;
}

void copy_snicar_host_views_d3(std::unordered_map<std::string, h_ViewD3>& snicar_host_views_d3,
  ELM::SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data)
{
  Kokkos::deep_copy(snicar_data.bcenh, snicar_host_views_d3["bcint_enh_mam"]);
}



std::unordered_map<std::string, h_ViewD3> get_snowage_host_views_d3(const ELM::SnwRdsTable<ViewD3>& snw_table)
{
  std::unordered_map<std::string, h_ViewD3> snowage_host_views_d3;
  snowage_host_views_d3["tau"] = Kokkos::create_mirror_view(snw_table.snowage_tau);
  snowage_host_views_d3["kappa"] = Kokkos::create_mirror_view(snw_table.snowage_kappa);
  snowage_host_views_d3["drdsdt0"] = Kokkos::create_mirror_view(snw_table.snowage_drdt0);

  return snowage_host_views_d3;
}

void copy_snowage_host_views_d3(std::unordered_map<std::string, h_ViewD3>& snowage_host_views_d3,
  ELM::SnwRdsTable<ViewD3>& snw_table)
{
  Kokkos::deep_copy(snw_table.snowage_tau, snowage_host_views_d3["tau"]);
  Kokkos::deep_copy(snw_table.snowage_kappa, snowage_host_views_d3["kappa"]);
  Kokkos::deep_copy(snw_table.snowage_drdt0, snowage_host_views_d3["drdsdt0"]);
}



std::unordered_map<std::string, h_ViewD1> get_pft_host_views(const ELM::PFTData<ViewD1>& pft_data)
{
  std::unordered_map<std::string, h_ViewD1> pft_host_views;
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

void copy_pft_host_views(std::unordered_map<std::string, h_ViewD1>& pft_host_views, ELM::PFTData<ViewD1>& pft_data)
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



std::unordered_map<std::string, h_ViewD1> get_aero_host_views(const ELM::AerosolDataManager<ViewD1>& aero_data)
{
  std::unordered_map<std::string, h_ViewD1> aero_host_views;
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

void copy_aero_host_views(std::unordered_map<std::string, h_ViewD1>& aero_host_views, ELM::AerosolDataManager<ViewD1>& aero_data)
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

} // namespace ELM::kokkos_init


void ELM::initialize_kokkos_elm (
  ELMStateType& S,
  SnicarData<ViewD1, ViewD2, ViewD3>& snicar_data,
  SnwRdsTable<ViewD3>& snw_rds_table,
  PFTData<ViewD1>& pft_data,
  AerosolDataManager<ViewD1>& aerosol_data,
  const Utils::DomainDecomposition<2>& dd,
  const std::string& fname_surfdata,
  const std::string& fname_param,
  const std::string& fname_snicar,
  const std::string& fname_snowage,
  const std::string& fname_aerosol)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevgrnd;

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // read data from files into host mirrors
  // and then copy to kokkos device
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

  {
    // soil color constants
    auto h_isoicol = Kokkos::create_mirror_view(S.isoicol);
    auto h_albsat = Kokkos::create_mirror_view(S.albsat);
    auto h_albdry = Kokkos::create_mirror_view(S.albdry);
    ELM::read_soil::read_soil_colors(dd, fname_surfdata, h_isoicol, h_albsat, h_albdry);
    if (S.albsat.extent(0) != h_albsat.extent(0) || S.albsat.extent(1) != h_albsat.extent(1))
    {
      NS::resize(S.albsat, h_albsat.extent(0), h_albsat.extent(1));
    }
    if (S.albdry.extent(0) != h_albdry.extent(0) || S.albdry.extent(1) != h_albdry.extent(1))
    {
      NS::resize(S.albdry, h_albdry.extent(0), h_albdry.extent(1));
    }
    Kokkos::deep_copy(S.isoicol, h_isoicol);
    Kokkos::deep_copy(S.albsat, h_albsat);
    Kokkos::deep_copy(S.albdry, h_albdry);
  }

  const int ncells = S.snl.extent(0);
  auto pct_sand = ELM::Utils::create<ViewD2>("pct_sand", ncells, nlevgrnd); // only used in init_soil_hydraulics()
  auto pct_clay = ELM::Utils::create<ViewD2>("pct_clay", ncells, nlevgrnd); // only used in init_soil_hydraulics()
  auto organic = ELM::Utils::create<ViewD2>("organic", ncells, nlevgrnd); // only used in init_soil_hydraulics()
  auto organic_max = ELM::Utils::create<ViewD1>("organic_max", 1); // only used in init_soil_hydraulics()

  {
    // soil texture constants
    auto h_pct_sand = Kokkos::create_mirror_view(pct_sand);
    auto h_pct_clay = Kokkos::create_mirror_view(pct_clay);
    auto h_organic = Kokkos::create_mirror_view(organic);
    auto h_organic_max = Kokkos::create_mirror_view(organic_max);
    read_soil::read_soil_texture(dd, fname_surfdata, fname_param,
                    h_organic_max, h_pct_sand, h_pct_clay, h_organic);
    Kokkos::deep_copy(pct_sand, h_pct_sand);
    Kokkos::deep_copy(pct_clay, h_pct_clay);
    Kokkos::deep_copy(organic, h_organic);
    Kokkos::deep_copy(organic_max, h_organic_max);
  }

  // snicar radiation parameters
  {
    auto host_snicar_d1 = kokkos_init::get_snicar_host_views_d1(snicar_data);
    auto host_snicar_d2 = kokkos_init::get_snicar_host_views_d2(snicar_data);
    auto host_snicar_d3 = kokkos_init::get_snicar_host_views_d3(snicar_data);
    read_snicar_data(host_snicar_d1, host_snicar_d2,
                          host_snicar_d3, dd.comm, fname_snicar);
    kokkos_init::copy_snicar_host_views_d1(host_snicar_d1, snicar_data);
    kokkos_init::copy_snicar_host_views_d2(host_snicar_d2, snicar_data);
    kokkos_init::copy_snicar_host_views_d3(host_snicar_d3, snicar_data);
  }

  // snow aging parameters
  // this is a large lookup table
  // need to figure out appropriate method to handle 
  // this table and Mie aerosol model lookup table
  {
    auto host_snowage_d3 = kokkos_init::get_snowage_host_views_d3(snw_rds_table);
    read_snowrds_data(host_snowage_d3, dd.comm, fname_snowage);
    kokkos_init::copy_snowage_host_views_d3(host_snowage_d3, snw_rds_table);
  }

  // pft data constants
  {
    auto host_pft_views = kokkos_init::get_pft_host_views(pft_data);
    read_pft_data(host_pft_views, dd.comm, fname_param);
    kokkos_init::copy_pft_host_views(host_pft_views, pft_data);
  }

  // aerosol deposition data manager
  {
    auto host_aero_views = kokkos_init::get_aero_host_views(aerosol_data);
    read_aerosol_data(host_aero_views, dd.comm, fname_aerosol, S.lon, S.lat);
    kokkos_init::copy_aero_host_views(host_aero_views, aerosol_data);
  }

  // Kokkos view of struct PSNVegData
  auto psn_pft = ELM::Utils::create<Kokkos::View<PFTDataPSN *>>("psn_pft", ncells);

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // free functions that calculate initial values
  // need to be run once for each cell
  // wrap in parallel lambda
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto init_functions = ELM_LAMBDA (const int& idx) {

    S.psn_pft(idx) = pft_data.get_pft_psn(S.vtype(idx));
  
    init_topo_slope(S.topo_slope(idx));

    init_micro_topo(
                    S.Land.ltype, S.topo_slope(idx),
                    S.topo_std(idx), S.n_melt(idx),
                    S.micro_sigma(idx));

    init_snow_layers(
                    S.snow_depth(idx), S.Land.lakpoi, S.snl(idx),
                    Kokkos::subview(S.dz, idx, Kokkos::ALL),
                    Kokkos::subview(S.zsoi, idx, Kokkos::ALL),
                    Kokkos::subview(S.zisoi, idx, Kokkos::ALL));

    init_soil_hydraulics(
                        organic_max(0),
                        Kokkos::subview(pct_sand, idx, Kokkos::ALL),
                        Kokkos::subview(pct_clay, idx, Kokkos::ALL),
                        Kokkos::subview(organic, idx, Kokkos::ALL),
                        Kokkos::subview(S.zsoi, idx, Kokkos::ALL),
                        Kokkos::subview(S.watsat, idx, Kokkos::ALL),
                        Kokkos::subview(S.bsw, idx, Kokkos::ALL),
                        Kokkos::subview(S.sucsat, idx, Kokkos::ALL),
                        Kokkos::subview(S.watdry, idx, Kokkos::ALL),
                        Kokkos::subview(S.watopt, idx, Kokkos::ALL),
                        Kokkos::subview(S.watfc, idx, Kokkos::ALL),
                        Kokkos::subview(S.tkmg, idx, Kokkos::ALL),
                        Kokkos::subview(S.tkdry, idx, Kokkos::ALL),
                        Kokkos::subview(S.csol, idx, Kokkos::ALL));


    init_vegrootfr(
                  S.vtype(idx), pft_data.roota_par(S.vtype(idx)),
                  pft_data.rootb_par(S.vtype(idx)),
                  Kokkos::subview(S.zisoi, idx, Kokkos::ALL),
                  Kokkos::subview(S.rootfr, idx, Kokkos::ALL));

    init_soil_temp(
                  S.Land, S.snl(idx),
                  Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
                  S.t_grnd(idx));

    init_snow_state(
                   S.Land.urbpoi, S.snl(idx), S.h2osno(idx), S.int_snow(idx),
                   S.snow_depth(idx), S.h2osfc(idx), S.h2ocan(idx), S.frac_h2osfc(idx),
                   S.fwet(idx), S.fdry(idx), S.frac_sno(idx),
                   Kokkos::subview(S.snw_rds, idx, Kokkos::ALL));

    init_soilh2o_state(
                      S.Land, S.snl(idx),
                      Kokkos::subview(S.watsat, idx, Kokkos::ALL),
                      Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
                      Kokkos::subview(S.dz, idx, Kokkos::ALL),
                      Kokkos::subview(S.h2osoi_vol, idx, Kokkos::ALL),
                      Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
                      Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL));
  };
  invoke_kernel(init_functions, std::make_tuple(ncells), "init functions");
}
