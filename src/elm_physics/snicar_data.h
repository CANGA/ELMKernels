
#pragma once

#include "array.hh"
#include "read_input.hh"
#include "utils.hh"

#include <array>
#include <string>
#include <map>
#include <vector>

#include "snow_snicar.h"
#include "kokkos_includes.hh"

namespace ELM::snicar_utils {

// read and manipulate data
template <typename ArrayT, typename T, size_t D>
void read_and_fill_array(const Comm_type& comm, const std::string& filename, const std::string& varname,
                         Array<T, D>& arr_for_read, ArrayT& arrout);

} // namespace ELM::snicar_utils

namespace ELM {

template <typename ArrayD1, typename ArrayD2, typename ArrayD3>
struct SnicarData {

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

  SnicarData();
};

// read all fields in SnicarData
template <typename h_ArrayD1, typename h_ArrayD2, typename h_ArrayD3>
void read_snicar_data(
  std::map<std::string, h_ArrayD1>& snicar_views_d1,
  std::map<std::string, h_ArrayD2>& snicar_views_d2,
  std::map<std::string, h_ArrayD3>& snicar_views_d3,
  const Comm_type& comm, const std::string& filename);

} // namespace ELM

#include "snicar_data_impl.hh"
