
#pragma once

template <typename ArrayD1, typename ArrayD2, typename ArrayD3>
ELM::SnicarData<ArrayD1, ArrayD2, ArrayD3>::
SnicarData()
    : ss_alb_oc1("ss_alb_oc1", numrad_snw_),
      asm_prm_oc1("asm_prm_oc1", numrad_snw_),
      ext_cff_mss_oc1("ext_cff_mss_oc1", numrad_snw_),
      ss_alb_oc2("ss_alb_oc2", numrad_snw_),
      asm_prm_oc2("asm_prm_oc2", numrad_snw_),
      ext_cff_mss_oc2("ext_cff_mss_oc2", numrad_snw_),
      ss_alb_dst1("ss_alb_dst1", numrad_snw_),
      asm_prm_dst1("asm_prm_dst1", numrad_snw_),
      ext_cff_mss_dst1("ext_cff_mss_dst1", numrad_snw_),
      ss_alb_dst2("ss_alb_dst2", numrad_snw_),
      asm_prm_dst2("asm_prm_dst2", numrad_snw_),
      ext_cff_mss_dst2("ext_cff_mss_dst2", numrad_snw_),
      ss_alb_dst3("ss_alb_dst3", numrad_snw_),
      asm_prm_dst3("asm_prm_dst3", numrad_snw_),
      ext_cff_mss_dst3("ext_cff_mss_dst3", numrad_snw_),
      ss_alb_dst4("ss_alb_dst4", numrad_snw_),
      asm_prm_dst4("asm_prm_dst4", numrad_snw_),
      ext_cff_mss_dst4("ext_cff_mss_dst4", numrad_snw_),
      ss_alb_snw_drc("ss_alb_snw_drc", numrad_snw_, idx_Mie_snw_mx_),
      asm_prm_snw_drc("asm_prm_snw_drc", numrad_snw_, idx_Mie_snw_mx_),
      ext_cff_mss_snw_drc("ext_cff_mss_snw_drc", numrad_snw_, idx_Mie_snw_mx_),
      ss_alb_snw_dfs("ss_alb_snw_dfs", numrad_snw_, idx_Mie_snw_mx_),
      asm_prm_snw_dfs("asm_prm_snw_dfs", numrad_snw_, idx_Mie_snw_mx_),
      ext_cff_mss_snw_dfs("ext_cff_mss_snw_dfs", numrad_snw_, idx_Mie_snw_mx_),
      ss_alb_bc1("ss_alb_bc1", idx_bc_nclrds_max_ + 1, numrad_snw_),
      asm_prm_bc1("asm_prm_bc1", idx_bc_nclrds_max_ + 1, numrad_snw_),
      ext_cff_mss_bc1("ext_cff_mss_bc1", idx_bc_nclrds_max_ + 1, numrad_snw_),
      ss_alb_bc2("ss_alb_bc2", idx_bc_nclrds_max_ + 1, numrad_snw_),
      asm_prm_bc2("asm_prm_bc2", idx_bc_nclrds_max_ + 1, numrad_snw_),
      ext_cff_mss_bc2("ext_cff_mss_bc2", idx_bc_nclrds_max_ + 1, numrad_snw_),
      bcenh("bcenh", idx_bcint_icerds_max_ + 1, idx_bc_nclrds_max_ + 1, numrad_snw_)
    {}


template <typename ArrayD3>
ELM::SnwRdsTable<ArrayD3>::
SnwRdsTable()
    : snowage_tau("snowage_tau", idx_T_max_+1,idx_Tgrd_max_+1,idx_rhos_max_+1),
      snowage_kappa("snowage_kappa", idx_T_max_+1,idx_Tgrd_max_+1,idx_rhos_max_+1),
      snowage_drdt0("snowage_drdt0", idx_T_max_+1,idx_Tgrd_max_+1,idx_rhos_max_+1)
    {}


template <typename h_ArrayD1, typename h_ArrayD2, typename h_ArrayD3>
void ELM::read_snicar_data(
  std::unordered_map<std::string, h_ArrayD1>& snicar_views_d1,
  std::unordered_map<std::string, h_ArrayD2>& snicar_views_d2,
  std::unordered_map<std::string, h_ArrayD3>& snicar_views_d3,
  const Comm_type& comm, const std::string& filename)
{
  using namespace ELM::snicar_utils;
  using ELMdims::numrad_snw;
  using snow_snicar::detail::idx_Mie_snw_mx;
  using snow_snicar::detail::idx_bc_nclrds_max;
  using snow_snicar::detail::idx_bcint_icerds_max;

  // read 1D arrays of size numrad_snw()
  {
    std::array<size_t, 1> start = {0};
    std::array<size_t, 1> count = {numrad_snw()};
    Array<double, 1> arr_for_read(numrad_snw());

    for (auto& [varname, arr] : snicar_views_d1) {
      read_and_fill_array(comm, filename, varname, arr_for_read, arr);
    }
  }

  // read 2D arrays of size [numrad_snw(), idx_Mie_snw_mx]
  {
    std::array<size_t, 2> start = {0};
    std::array<size_t, 2> count = {numrad_snw(), idx_Mie_snw_mx};
    Array<double, 2> arr_for_read(numrad_snw(), idx_Mie_snw_mx);

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

  // read 2D arrays of size [idx_bc_nclrds_max+1, numrad_snw()]
  {
    std::array<size_t, 2> start = {0};
    std::array<size_t, 2> count = {idx_bc_nclrds_max + 1, numrad_snw()};
    Array<double, 2> arr_for_read(idx_bc_nclrds_max + 1, numrad_snw());

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

  // read 3D array of size [idx_bcint_icerds_max+1, idx_bc_nclrds_max+1, numrad_snw()]
  {
    std::array<size_t, 3> start = {0};
    std::array<size_t, 3> count = {idx_bcint_icerds_max + 1, idx_bc_nclrds_max + 1,
                                  numrad_snw()};
    Array<double, 3> arr_for_read(idx_bcint_icerds_max + 1, idx_bc_nclrds_max + 1,
                                 numrad_snw());

    std::string varname("bcint_enh_mam");
    auto& arr = snicar_views_d3[varname];
    read_and_fill_array(comm, filename, varname, arr_for_read, arr);
  }
}


template <typename h_ArrayD3>
void ELM::read_snowrds_data(std::unordered_map<std::string, h_ArrayD3>& snowage_views_d3,
  const Comm_type& comm, const std::string& filename)
{
  using namespace ELM::snicar_utils;
  using snow_snicar::detail::idx_T_max;
  using snow_snicar::detail::idx_Tgrd_max;
  using snow_snicar::detail::idx_rhos_max;

  // read 3D arrays of size [idx_T_max+1, idx_Tgrd_max+1, idx_rhos_max+1]
  {
    std::array<size_t, 3> start = {0};
    std::array<size_t, 3> count = {idx_T_max + 1, idx_Tgrd_max + 1, idx_rhos_max + 1};
    Array<double, 3> arr_for_read(idx_T_max + 1, idx_Tgrd_max + 1, idx_rhos_max + 1);

    std::vector<std::string> names =
    {
      "tau",
      "kappa",
      "drdsdt0"
    };

    for (const auto& varname : names) {
      auto& arr = snowage_views_d3[varname];
      read_and_fill_array(comm, filename, varname, arr_for_read, arr);
    }
  }

}

template <typename ArrayT, typename T, size_t D>
inline void ELM::snicar_utils::
read_and_fill_array(const Comm_type& comm, const std::string& filename,
                    const std::string& varname, Array<T, D>& arr_for_read,
                    ArrayT& arrout)
{
  std::array<size_t, D> start{0};
  std::array<size_t, D> count;
  for (int i = 0; i != D; ++i)
    count[i] = static_cast<size_t>(arr_for_read.extent(i));
  IO::read_netcdf(comm, filename, varname, start, count, arr_for_read.data());
  if constexpr (D == 1) {
    for (int i = 0; i < arr_for_read.extent(0); ++i)
      arrout(i) = arr_for_read(i);
  } else if constexpr (D == 2) {
    for (int i = 0; i < arr_for_read.extent(0); ++i)
      for (int j = 0; j < arr_for_read.extent(1); ++j)
        arrout(i, j) = arr_for_read(i, j);
  } else if constexpr (D == 3) {
    for (int i = 0; i < arr_for_read.extent(0); ++i)
      for (int j = 0; j < arr_for_read.extent(1); ++j)
        for (int k = 0; k < arr_for_read.extent(2); ++k)
          arrout(i, j, k) = arr_for_read(i, j, k);
  }
  return;
}


