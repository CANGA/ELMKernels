
#pragma once

#include "monthly_data.h"
#include "elm_constants.h"

#include "array.hh"
#include "utils.hh"
#include "read_input.hh"

#include <iostream>
#include <map>
#include <tuple>
#include <cmath>
#include <string>
#include <array>

// this is derived from SatellitePhenologyMod.F90
namespace ELM::phenology_data {

// functor to calculate phenology parameters for time = model_time
template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
struct ComputePhenology {
  ComputePhenology(const ArrayD2& mlai, const ArrayD2& msai, const ArrayD2& mhtop, const ArrayD2& mhbot,
    const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype, const double& wt1,
    const double wt2, const int start_idx, ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot,
    ArrayD1& tlai, ArrayD1&tsai, ArrayI1& frac_veg_nosno_alb)
      : mlai_(mlai),
        msai_(msai),
        mhtop_(mhtop),
        mhbot_(mhbot),
        snow_depth_(snow_depth),
        frac_sno_(frac_sno),
        vtype_(vtype),
        wt1_(wt1),
        wt2_(wt2),
        start_idx_(start_idx),
        elai_(elai),
        esai_(esai),
        htop_(htop),
        hbot_(hbot),
        tlai_(tlai),
        tsai_(tsai),
        frac_veg_nosno_alb_(frac_veg_nosno_alb) {}

  void operator()(const int i) const {
  // leaf phenology
  // Set leaf and stem areas based on day of year
  // Interpolate leaf area index, stem area index, and vegetation heights
  // between two monthly values using weights, (wt1, wt2)
  if (vtype_(i) != ELM::noveg) {
    tlai_(i) = wt1_ * mlai_(start_idx_, i) + wt2_ * mlai_(start_idx_+1, i);
    tsai_(i) = wt1_ * msai_(start_idx_, i) + wt2_ * msai_(start_idx_+1, i);
    htop_(i) = wt1_ * mhtop_(start_idx_, i) + wt2_ * mhtop_(start_idx_+1, i);
    hbot_(i) = wt1_ * mhbot_(start_idx_, i) + wt2_ * mhbot_(start_idx_+1, i);
  } else {
    tlai_(i) =  0.0;
    tsai_(i) =  0.0;
    htop_(i) =  0.0;
    hbot_(i) =  0.0;
  }
  // adjust lai and sai for burying by snow. if exposed lai and sai
  // are less than 0.05, set equal to zero to prevent numerical
  // problems associated with very small lai and sai.
  // snow burial fraction for short vegetation (e.g. grasses) as in
  // Wang and Zeng, 2007.
  double fb;
  if (vtype_(i) > ELM::noveg && vtype_(i) <= ELM::nbrdlf_dcd_brl_shrub) {
    double ol = std::min(std::max(snow_depth_(i) - hbot_(i), 0.0), htop_(i) - hbot_(i));
    fb = 1.0 - ol / std::max(1.e-06, htop_(i) - hbot_(i));
  } else {
    // 0.2m is assumed depth of snow required for complete burial of grasses
    fb = 1.0 - std::max(std::min(snow_depth_(i), 0.2), 0.0) / 0.2; 
  }
  // area weight by snow covered fraction
  elai_(i) = std::max(tlai_(i) * (1.0 - frac_sno_(i)) + tlai_(i) * fb * frac_sno_(i), 0.0);
  esai_(i) = std::max(tsai_(i) * (1.0 - frac_sno_(i)) + tsai_(i) * fb * frac_sno_(i), 0.0);
  if (elai_(i) < 0.05) {
    elai_(i) = 0.0;
  }
  if (esai_(i) < 0.05) {
    esai_(i) = 0.0;
  }
  // Fraction of vegetation free of snow
  if ((elai_(i) + esai_(i)) >= 0.05) {
    frac_veg_nosno_alb_(i) = 1;
  } else {
    frac_veg_nosno_alb_(i) = 0;
  }
}

private:
  ArrayD2 mlai_, msai_, mhtop_, mhbot_;
  ArrayD1 snow_depth_, frac_sno_;
  ArrayI1 vtype_;
  double wt1_, wt2_;
  int start_idx_;
  ArrayD1 elai_, esai_, htop_, hbot_, tlai_, tsai_;
  ArrayI1 frac_veg_nosno_alb_;
};


// class to manage phenology data
class PhenologyDataManager : public MonthlyDataManager {

using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

public:

  PhenologyDataManager(const size_t ncells, const size_t npfts) 
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

  // read 1 month of data from file (1, npfts, nlat, nlon) for input param month 
  // and place into 2D array arr(ntimes, ncells) where ntimes = arr_idx
  // by combining nlat & nlon and filtering by pft type (vtype) - only one pft per grid cell currently
  void read_month(const Utils::DomainDecomposition<2> &dd, const std::string& filename, const std::string& varname,
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


  auto triple_month_indices(const Utils::Date& model_time) {
    auto [m1, m2] = month_indices(model_time);
    int m3 = m2 + 1;
    if (m3 > 11) { m3 = 0; }
    return std::make_tuple(m1, m2, m3);
  }

  void read_initial(const Utils::DomainDecomposition<2> &dd, const std::string& filename,
    const Utils::Date& model_time, const ArrayI1& vtype) {
    auto [m1, m2, m3] = triple_month_indices(model_time);
    std::cout << std::setprecision(15) << "m1,2,3 " << m1 << "  " << m2 << "  " << m3 << std::endl;
    std::map<int, int> months = {{0, m1}, {1, m2}, {2, m3}};
    for (const auto& [idx, mon] : months) {
      for (auto& [varname, arr] : phenology_names_) {
        read_month(dd, filename, varname, mon, vtype, idx, arr);
      }
    }
    data_m1_ = m1;
  }


  void read_new_month(const Utils::DomainDecomposition<2> &dd, const std::string& filename,
    const Utils::Date& model_time, const ArrayI1& vtype) {

        auto advance_month_idx = [] (ArrayD2& arr) {
          for (int mon = 0; mon < arr.extent(0)-1; ++mon) {
            for (int cell = 0; cell < arr.extent(1); ++cell) {
              arr(mon,cell) = arr(mon+1,cell);
            }
          }
        };
    auto [m1, m2, m3] = triple_month_indices(model_time);
    for (auto& [varname, arr] : phenology_names_) {
      advance_month_idx(arr);
      read_month(dd, filename, varname, m3, vtype, 2, arr);
    }
  }
  

  // read data from file
  // either all three months of data, or new single month
  void read_data(const Utils::DomainDecomposition<2> &dd, const std::string& filename, const Utils::Date& model_time,
    const ArrayI1& vtype) {
    if (!initialized_) {
      read_initial(dd, filename, model_time, vtype);
      initialized_ = true;
    } else if (need_new_data_) {
      auto [m1, m2] = month_indices(model_time);
      data_m1_ = m1;
      read_new_month(dd, filename, model_time, vtype);
    }
  }


  void get_data(const Utils::Date& model_time, const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype,
    ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot, ArrayD1& tlai, ArrayD1&tsai, ArrayI1& frac_veg_nosno_alb) {
    auto [wt1, wt2] = monthly_data_weights(model_time);
    auto [m1, m2] = month_indices(model_time);
    int start_idx = 0;
    if (data_m1_ != m1) { 
      start_idx = 1;
      need_new_data_= true;
    }
    ComputePhenology compute_phen(mlai_, msai_, mhtop_, mhbot_, snow_depth, frac_sno,
      vtype, wt1, wt2, start_idx, elai, esai, htop, hbot, tlai, tsai, frac_veg_nosno_alb);
    for (int i = 0; i < elai.extent(0); ++i) {
      std::invoke(compute_phen, i);
    }
  }

private:
  ArrayD2 mlai_, msai_, mhtop_, mhbot_;
  size_t ncells_, npfts_;
  bool initialized_, need_new_data_;
  std::map<std::string, ArrayD2&> phenology_names_;
  int data_m1_;
};

} // namespace ELM::phenology_data

//#include "phenology_data_impl.hh"
