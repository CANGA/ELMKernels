
#pragma once

#include "elm_constants.h"
#include "array.hh"
#include "read_input.hh"
#include "utils.hh"

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "snow_snicar.h"
#include "compile_options.hh"

namespace ELM::snicar_utils {

// read and manipulate data
template <typename ArrayT, typename T, size_t D>
void read_and_fill_array(const Comm_type& comm, const std::string& filename, const std::string& varname,
                         Array<T, D>& arr_for_read, ArrayT& arrout);

} // namespace ELM::snicar_utils

namespace ELM {

template <typename ArrayD1, typename ArrayD2, typename ArrayD3>
struct SnicarData {

  const int numrad_snw_{ELMdims::numrad_snw};
  const int idx_Mie_snw_mx_{snow_snicar::detail::idx_Mie_snw_mx};
  const int idx_bc_nclrds_max_{snow_snicar::detail::idx_bc_nclrds_max};
  const int idx_bcint_icerds_max_{snow_snicar::detail::idx_bcint_icerds_max};

  SnicarData();
  ~SnicarData() = default;

  // arrays to store snicar parameters
  ArrayD1 ss_alb_oc1;
  ArrayD1 asm_prm_oc1;
  ArrayD1 ext_cff_mss_oc1;
  ArrayD1 ss_alb_oc2;
  ArrayD1 asm_prm_oc2;
  ArrayD1 ext_cff_mss_oc2;
  ArrayD1 ss_alb_dst1;
  ArrayD1 asm_prm_dst1;
  ArrayD1 ext_cff_mss_dst1;
  ArrayD1 ss_alb_dst2;
  ArrayD1 asm_prm_dst2;
  ArrayD1 ext_cff_mss_dst2;
  ArrayD1 ss_alb_dst3;
  ArrayD1 asm_prm_dst3;
  ArrayD1 ext_cff_mss_dst3;
  ArrayD1 ss_alb_dst4;
  ArrayD1 asm_prm_dst4;
  ArrayD1 ext_cff_mss_dst4;
  ArrayD2 ss_alb_snw_drc;
  ArrayD2 asm_prm_snw_drc;
  ArrayD2 ext_cff_mss_snw_drc;
  ArrayD2 ss_alb_snw_dfs;
  ArrayD2 asm_prm_snw_dfs;
  ArrayD2 ext_cff_mss_snw_dfs;
  ArrayD2 ss_alb_bc1;
  ArrayD2 asm_prm_bc1;
  ArrayD2 ext_cff_mss_bc1;
  ArrayD2 ss_alb_bc2;
  ArrayD2 asm_prm_bc2;
  ArrayD2 ext_cff_mss_bc2;
  ArrayD3 bcenh;
};


template <typename ArrayD3>
struct SnwRdsTable {
  const int idx_T_max_{snow_snicar::detail::idx_T_max};
  const int idx_Tgrd_max_{snow_snicar::detail::idx_Tgrd_max};
  const int idx_rhos_max_{snow_snicar::detail::idx_rhos_max};
  SnwRdsTable();
  ~SnwRdsTable() = default;
  ArrayD3 snowage_tau;
  ArrayD3 snowage_kappa;
  ArrayD3 snowage_drdt0;
};

// read all fields in SnicarData
template <typename h_ArrayD1, typename h_ArrayD2, typename h_ArrayD3>
void read_snicar_data(
  std::unordered_map<std::string, h_ArrayD1>& snicar_views_d1,
  std::unordered_map<std::string, h_ArrayD2>& snicar_views_d2,
  std::unordered_map<std::string, h_ArrayD3>& snicar_views_d3,
  const Comm_type& comm, const std::string& filename);

template <typename h_ArrayD3>
void read_snowrds_data(std::unordered_map<std::string, h_ArrayD3>& snowage_views_d3,
  const Comm_type& comm, const std::string& filename);

} // namespace ELM

#include "snicar_data_impl.hh"
