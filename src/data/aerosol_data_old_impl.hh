
#pragma once

namespace ELM::aerosol_utils {

template <typename h_ArrayD1>
void read_variable_slice(const Comm_type& comm, const std::string& filename,
                         const std::string& varname, const size_t& lon_idx,
                         const size_t& lat_idx, h_ArrayD1 arr)
{
  Array<double, 3> file_data(12, 1, 1); // one year of monthly data
  std::array<size_t, 3> start{0, lat_idx, lon_idx};
  std::array<size_t, 3> count{12, 1, 1};
  IO::read_netcdf(comm, filename, varname, start, count, file_data.data());

  for (int i = 0; i < 12; ++i) {
    arr(i) = file_data(i, 0, 0);
  }
}

} // namespace ELM::aerosol_utils

template <typename ArrayD1>
ELM::AerosolDataManager<ArrayD1>::AerosolDataManager()
    : bcdep("bcdep", 12), bcpho("bcpho", 12), bcphi("bcphi", 12),
      dst1_1("dst1_1", 12), dst1_2("dst1_2", 12), dst2_1("dst2_1", 12),
      dst2_2("dst2_2", 12), dst3_1("dst3_1", 12), dst3_2("dst3_2", 12),
      dst4_1("dst4_1", 12), dst4_2("dst4_2", 12)
    {}

template <typename ArrayD1>
auto ELM::AerosolDataManager<ArrayD1>::
get_aerosol_source(const Utils::Date& model_time, const double& dtime) const
{
  auto [wt1, wt2] = monthly_data::monthly_data_weights(model_time);
  auto [m1, m2] = monthly_data::month_indices(model_time);
  double forc_bcphi = (wt1 * bcphi(m1) + wt2 * bcphi(m2)) * dtime;
  double forc_bcpho = (wt1 * (bcdep(m1) + bcpho(m1)) + wt2 * (bcdep(m2) + bcpho(m2))) * dtime;
  double forc_dst1 = (wt1 * (dst1_1(m1) + dst1_2(m1)) + wt2 * (dst1_1(m2) + dst1_2(m2))) * dtime;
  double forc_dst2 = (wt1 * (dst2_1(m1) + dst2_2(m1)) + wt2 * (dst2_1(m2) + dst2_2(m2))) * dtime;
  double forc_dst3 = (wt1 * (dst3_1(m1) + dst3_2(m1)) + wt2 * (dst3_1(m2) + dst3_2(m2))) * dtime;
  double forc_dst4 = (wt1 * (dst4_1(m1) + dst4_2(m1)) + wt2 * (dst4_1(m2) + dst4_2(m2))) * dtime;
  return std::make_tuple(forc_bcphi, forc_bcpho, forc_dst1, forc_dst2, forc_dst3, forc_dst4);
}

template <typename h_ArrayD1>
void ELM::read_aerosol_data(std::unordered_map<std::string, h_ArrayD1>& aerosol_views, 
  const Comm_type& comm, const std::string& filename, const double& lon_d, const double& lat_d)
{
  const auto [lon_idx, lat_idx] = aerosol_utils::get_nearest_indices(comm, filename, lon_d, lat_d);

  for (auto& [varname, arr] : aerosol_views) {
    aerosol_utils::read_variable_slice(comm, filename, varname, lon_idx, lat_idx, arr);
  }
}
