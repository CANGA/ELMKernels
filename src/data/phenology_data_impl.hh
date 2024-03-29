
#pragma once


namespace ELM {

// this is derived from SatellitePhenologyMod.F90
template <typename ArrayD2>
PhenologyDataManager<ArrayD2>::
PhenologyDataManager(const Utils::DomainDecomposition<2>& dd,
                     const size_t& ncells, const size_t& npfts)
    : mlai("mlai", 3, ncells), msai("msai", 3, ncells),
      mhtop("mhtop", 3, ncells), mhbot("mhbot", 3, ncells),
      dd_{dd}, ncells_{ncells}, npfts_{npfts}, initialized_{false},
      need_new_data_{false}
    {}

// read data from file
// either all three months of data, or new single month
template <typename ArrayD2>
template <typename ArrayI1, typename h_ArrayD2>
bool PhenologyDataManager<ArrayD2>::
read_data(std::unordered_map<std::string, h_ArrayD2>& phenology_views, const std::string& filename,
          const Utils::Date& model_time, const ArrayI1 vtype)
{
  if (!initialized_) {
    read_initial(phenology_views, filename, model_time, vtype);
    initialized_ = true;
    need_new_data_ = false;
    auto m1 = monthly_data::first_month_idx(model_time);
    data_m1_ = m1;
    return true;
  } else if (need_new_data_) {
    read_new_month(phenology_views, filename, model_time, vtype);
    need_new_data_ = false;
    auto m1 = monthly_data::first_month_idx(model_time);
    data_m1_ = m1;
    return true;
  }
  return false;
}

template <typename ArrayD2>
template <typename ArrayI1, typename ArrayD1>
void PhenologyDataManager<ArrayD2>::
get_data(const Utils::Date& model_time, const ArrayD1 snow_depth,
         const ArrayD1 frac_sno, const ArrayI1 vtype, ArrayD1 elai, ArrayD1 esai,
         ArrayD1 htop, ArrayD1 hbot, ArrayD1 tlai, ArrayD1 tsai,
         ArrayI1 frac_veg_nosno_alb)
{
  auto [wt1, wt2] = monthly_data::monthly_data_weights(model_time);
  auto m1 = monthly_data::first_month_idx(model_time);
  int start_idx = 0;
  if (data_m1_ != m1) {
    start_idx += 1;
    need_new_data_ = true;
  }

  phenology::ComputePhenology compute_phen(mlai, msai, mhtop, mhbot, snow_depth, frac_sno, vtype, wt1, wt2,
                                           start_idx, elai, esai, htop, hbot, tlai, tsai, frac_veg_nosno_alb);

  apply_parallel_for(compute_phen, "ComputePhenology", elai.extent(0));
}

// read 1 month of data from file (1, npfts, nlat, nlon) for input param month
// and place into 2D array arr(ntimes, ncells) where ntimes = arr_idx
// by combining nlat & nlon and filtering by pft type (vtype) - only one pft per grid cell currently
template <typename ArrayD2>
template <typename ArrayI1, typename h_ArrayD2>
void PhenologyDataManager<ArrayD2>::
read_month(const std::string& filename, const std::string& varname,
           const size_t& month, const ArrayI1 vtype,
           const int& arr_idx, h_ArrayD2 arr)
{
  // allocate one month of data
  Array<double, 4> arr_for_read(1, npfts_, dd_.n_local[0], dd_.n_local[1]);
  std::array<size_t, 4> start = {month, 0, dd_.start[0], dd_.start[1]};
  std::array<size_t, 4> count = {1, npfts_, dd_.n_local[0], dd_.n_local[1]};
  IO::read_netcdf(dd_.comm, filename, varname, start, count, arr_for_read.data());
  for (size_t i = 0; i != dd_.n_local[0]; ++i) {
    for (size_t j = 0; j != dd_.n_local[1]; ++j) {
      size_t ncell_idx = i * dd_.n_local[1] + j;
      int pft = vtype(ncell_idx);
      arr(arr_idx, ncell_idx) = arr_for_read(0, pft, i, j);
    }
  }
}

template <typename ArrayD2>
template <typename ArrayI1, typename h_ArrayD2>
void PhenologyDataManager<ArrayD2>::
read_initial(std::unordered_map<std::string, h_ArrayD2>& phenology_views,
             const std::string& filename, const Utils::Date& model_time,
             const ArrayI1 vtype)
{
  const auto [m1, m2, m3] = monthly_data::triple_month_indices(model_time);
  std::unordered_map<int, int> months = {{0, m1}, {1, m2}, {2, m3}};
  for (auto [idx, mon] : months) {
    for (auto& [varname, arr] : phenology_views) {
      read_month(filename, varname, mon, vtype, idx, arr);
    }
  }
}

template <typename ArrayD2>
template <typename ArrayI1, typename h_ArrayD2>
void PhenologyDataManager<ArrayD2>::read_new_month(std::unordered_map<std::string, h_ArrayD2>& phenology_views, const std::string& filename,
                                               const Utils::Date& model_time, const ArrayI1 vtype) {

  auto advance_month_idx = [] (h_ArrayD2& arr) {
    for (size_t mon = 0; mon < arr.extent(0) - 1; ++mon) {
      for (size_t cell = 0; cell < arr.extent(1); ++cell) {
        arr(mon, cell) = arr(mon + 1, cell);
      }
    }
  };
  const auto m3 = monthly_data::third_month_idx(model_time);
  for (auto& [varname, arr] : phenology_views) {
    advance_month_idx(arr);
    read_month(filename, varname, m3, vtype, 2, arr);
  }
}

} // namespace ELM
