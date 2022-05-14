
#include <map>
#include <vector>
#include <string>

#include "snicar_data.h"
#include "snow_snicar.h"

ELM::SnicarData::SnicarData()
    : ss_alb_oc1("ss_alb_oc1", snow_snicar::numrad_snw),
      asm_prm_oc1("asm_prm_oc1", snow_snicar::numrad_snw),
      ext_cff_mss_oc1("ext_cff_mss_oc1", snow_snicar::numrad_snw),
      ss_alb_oc2("ss_alb_oc2", snow_snicar::numrad_snw),
      asm_prm_oc2("asm_prm_oc2", snow_snicar::numrad_snw),
      ext_cff_mss_oc2("ext_cff_mss_oc2", snow_snicar::numrad_snw),
      ss_alb_dst1("ss_alb_dst1", snow_snicar::numrad_snw),
      asm_prm_dst1("asm_prm_dst1", snow_snicar::numrad_snw),
      ext_cff_mss_dst1("ext_cff_mss_dst1", snow_snicar::numrad_snw),
      ss_alb_dst2("ss_alb_dst2", snow_snicar::numrad_snw),
      asm_prm_dst2("asm_prm_dst2", snow_snicar::numrad_snw),
      ext_cff_mss_dst2("ext_cff_mss_dst2", snow_snicar::numrad_snw),
      ss_alb_dst3("ss_alb_dst3", snow_snicar::numrad_snw),
      asm_prm_dst3("asm_prm_dst3", snow_snicar::numrad_snw),
      ext_cff_mss_dst3("ext_cff_mss_dst3", snow_snicar::numrad_snw),
      ss_alb_dst4("ss_alb_dst4", snow_snicar::numrad_snw),
      asm_prm_dst4("asm_prm_dst4", snow_snicar::numrad_snw),
      ext_cff_mss_dst4("ext_cff_mss_dst4", snow_snicar::numrad_snw),
      ss_alb_snw_drc("ss_alb_snw_drc", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
      asm_prm_snw_drc("asm_prm_snw_drc", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
      ext_cff_mss_snw_drc("ext_cff_mss_snw_drc", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
      ss_alb_snw_dfs("ss_alb_snw_dfs", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
      asm_prm_snw_dfs("asm_prm_snw_dfs", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
      ext_cff_mss_snw_dfs("ext_cff_mss_snw_dfs", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
      ss_alb_bc1("ss_alb_bc1", ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw),
      asm_prm_bc1("asm_prm_bc1", ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw),
      ext_cff_mss_bc1("ext_cff_mss_bc1", ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw),
      ss_alb_bc2("ss_alb_bc2", ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw),
      asm_prm_bc2("asm_prm_bc2", ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw),
      ext_cff_mss_bc2("ext_cff_mss_bc2", ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw),
      bcenh("bcenh", snow_snicar::idx_bcint_icerds_max + 1, snow_snicar::idx_bc_nclrds_max + 1, snow_snicar::numrad_snw)
    {}

void ELM::read_snicar_data(
  std::map<std::string, h_ArrayD1>& snicar_views_d1,
  std::map<std::string, h_ArrayD2>& snicar_views_d2,
  std::map<std::string, h_ArrayD3>& snicar_views_d3,
  const Comm_type& comm, const std::string& filename) {

  using namespace ELM::snicar_utils;

  // read 1D arrays of size snow_snicar::numrad_snw
  {
    std::array<size_t, 1> start = {0};
    std::array<size_t, 1> count = {snow_snicar::numrad_snw};
    Array<double, 1> arr_for_read(snow_snicar::numrad_snw);

    for (auto& [varname, arr] : snicar_views_d1) {
      read_and_fill_array(comm, filename, varname, arr_for_read, arr);
    }
  }

  // read 2D arrays of size [numrad_snw, idx_Mie_snw_mx]
  {
    std::array<size_t, 2> start = {0};
    std::array<size_t, 2> count = {snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx};
    Array<double, 2> arr_for_read(snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx);

    std::vector<std::string> names =
    {
      "ss_alb_ice_drc",
      "asm_prm_ice_drc",
      "ext_cff_mss_ice_drc",
      "ss_alb_ice_dfs",
      "asm_prm_ice_dfs",
      "ext_cff_mss_ice_dfs"
    };

    for (const auto& varname : names) {
      auto& arr = snicar_views_d2[varname];
      read_and_fill_array(comm, filename, varname, arr_for_read, arr);
    }
  }

  // read 2D arrays of size [idx_bc_nclrds_max+1, numrad_snw]
  {
    std::array<size_t, 2> start = {0};
    std::array<size_t, 2> count = {ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw};
    Array<double, 2> arr_for_read(ELM::snow_snicar::idx_bc_nclrds_max + 1, ELM::snow_snicar::numrad_snw);

    std::vector<std::string> names =
    {
      "ss_alb_bc_mam",
      "asm_prm_bc_mam",
      "ext_cff_mss_bc_mam",
      "ss_alb_bc_mam",
      "asm_prm_bc_mam",
      "ext_cff_mss_bc_mam"
    };

    for (const auto& varname : names) {
      auto& arr = snicar_views_d2[varname];
      read_and_fill_array(comm, filename, varname, arr_for_read, arr);
    }
  }

  // read 3D array of size [idx_bcint_icerds_max+1, idx_bc_nclrds_max+1, numrad_snw]
  {
    std::array<size_t, 3> start = {0};
    std::array<size_t, 3> count = {snow_snicar::idx_bcint_icerds_max + 1, snow_snicar::idx_bc_nclrds_max + 1,
                                   snow_snicar::numrad_snw};
    Array<double, 3> arr_for_read(snow_snicar::idx_bcint_icerds_max + 1, snow_snicar::idx_bc_nclrds_max + 1,
                                  snow_snicar::numrad_snw);

    std::string varname("bcint_enh_mam");
    auto& arr = snicar_views_d3[varname];
    read_and_fill_array(comm, filename, varname, arr_for_read, arr);
  }
}
