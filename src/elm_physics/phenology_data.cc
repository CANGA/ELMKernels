
#include "phenology_data.h"

// this is derived from SatellitePhenologyMod.F90
namespace ELM::phenology_data {


  PhenologyDataManager::PhenologyDataManager(const size_t ncells, const size_t npfts) 
    : mlai_("mlai", 3, ncells),
      msai_("msai", 3, ncells),
      mhtop_("mhtop", 3, ncells),
      mhbot_("mhbot", 3, ncells),
      ncells_(ncells),
      npfts_(npfts),
      initialized_(false),
      need_new_data_(false),
      phenology_names_{ {{"MONTHLY_LAI", mlai_}, {"MONTHLY_SAI", msai_}, 
      {"MONTHLY_HEIGHT_TOP", mhtop_}, {"MONTHLY_HEIGHT_BOT", mhbot_}} } {}
  

  // read data from file
  // either all three months of data, or new single month
  void PhenologyDataManager::read_data(const Utils::DomainDecomposition<2> &dd, const std::string& filename, const Utils::Date& model_time,
    const ArrayI1& vtype) {
    if (!initialized_) {
      read_initial(dd, filename, model_time, vtype);
      initialized_ = true;
      auto m1 = monthly_data::first_month_idx(model_time);
      data_m1_ = m1;
    } else if (need_new_data_) {
      read_new_month(dd, filename, model_time, vtype);
      need_new_data_ = false;
      auto m1 = monthly_data::first_month_idx(model_time);
      data_m1_ = m1;
    }
  }

  void PhenologyDataManager::get_data(const Utils::Date& model_time, const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype,
    ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot, ArrayD1& tlai, ArrayD1&tsai, ArrayI1& frac_veg_nosno_alb) {
    auto [wt1, wt2] = monthly_data::monthly_data_weights(model_time);
    auto m1 = monthly_data::first_month_idx(model_time);
    int start_idx = 0;
    if (data_m1_ != m1) { 
      start_idx += 1;
      need_new_data_= true;
    }
    ComputePhenology compute_phen(mlai_, msai_, mhtop_, mhbot_, snow_depth, frac_sno,
      vtype, wt1, wt2, start_idx, elai, esai, htop, hbot, tlai, tsai, frac_veg_nosno_alb);
    for (int i = 0; i < elai.extent(0); ++i) {
      std::invoke(compute_phen, i);
    }
  }

  // read 1 month of data from file (1, npfts, nlat, nlon) for input param month 
  // and place into 2D array arr(ntimes, ncells) where ntimes = arr_idx
  // by combining nlat & nlon and filtering by pft type (vtype) - only one pft per grid cell currently
  void PhenologyDataManager::read_month(const Utils::DomainDecomposition<2> &dd, const std::string& filename, const std::string& varname,
    const size_t month, const ArrayI1& vtype, const int arr_idx, ArrayD2& arr) {
    // allocate one month of data
    Array<double, 4> arr_for_read(1, npfts_, dd.n_local[0], dd.n_local[1]);
    std::array<size_t, 4> start = {month, 0, dd.start[0], dd.start[1]};
    std::array<size_t, 4> count = {1, npfts_, dd.n_local[0], dd.n_local[1]};
    IO::read_netcdf(dd.comm, filename, varname, start, count, arr_for_read.data());
    for (size_t i = 0; i != dd.n_local[0]; ++i) {
      for (size_t j = 0; j != dd.n_local[1]; ++j) {
        size_t ncell_idx = i * dd.n_local[1] + j;
        int pft = vtype(ncell_idx);
        arr(arr_idx, ncell_idx) = arr_for_read(0, pft, i, j);
      }
    }
  }

  void PhenologyDataManager::read_initial(const Utils::DomainDecomposition<2> &dd, const std::string& filename,
    const Utils::Date& model_time, const ArrayI1& vtype) {
    const auto [m1, m2, m3] = monthly_data::triple_month_indices(model_time);
    std::map<int, int> months = {{0, m1}, {1, m2}, {2, m3}};
    for (auto [idx, mon] : months) {
      for (auto& [varname, arr] : phenology_names_) {
        read_month(dd, filename, varname, mon, vtype, idx, arr);
      }
    }
  }

  void PhenologyDataManager::read_new_month(const Utils::DomainDecomposition<2> &dd, const std::string& filename,
    const Utils::Date& model_time, const ArrayI1& vtype) {

        auto advance_month_idx = [] (ArrayD2& arr) {
          for (int mon = 0; mon < arr.extent(0)-1; ++mon) {
            for (int cell = 0; cell < arr.extent(1); ++cell) {
              arr(mon,cell) = arr(mon+1,cell);
            }
          }
        };
    const auto m3 = monthly_data::third_month_idx(model_time);
    for (auto& [varname, arr] : phenology_names_) {
      advance_month_idx(arr);
      read_month(dd, filename, varname, m3, vtype, 2, arr);
    }
  }

} // namespace ELM::phenology_data

