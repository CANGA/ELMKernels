
#pragma once

#include "monthly_data.h"
#include "elm_constants.h"

#include "array.hh"
#include "date_time.hh"
#include "utils.hh"
#include "read_input.hh"

#include <map>
#include <tuple>
#include <cmath>
#include <string>
#include <array>
#include <functional>

// this is derived from SatellitePhenologyMod.F90
namespace ELM::phenology_data {

// functor to calculate phenology parameters for time = model_time
template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
struct ComputePhenology {
  ComputePhenology(const ArrayD2& mlai, const ArrayD2& msai, const ArrayD2& mhtop, const ArrayD2& mhbot,
    const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype, const double& wt1,
    const double wt2, const int start_idx, ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot,
    ArrayD1& tlai, ArrayD1& tsai, ArrayI1& frac_veg_nosno_alb);

  void operator()(const int i) const;

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
class PhenologyDataManager {

using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

public:
  PhenologyDataManager(const size_t ncells, const size_t npfts);

  // read data from file
  // either all three months of data, or new single month
  void read_data(const Utils::DomainDecomposition<2> &dd, const std::string& filename, const Utils::Date& model_time,
    const ArrayI1& vtype);

  // get phenology data - call parallel physics kernel - return phenology data for this timestep
  void get_data(const Utils::Date& model_time, const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype,
    ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot, ArrayD1& tlai, ArrayD1&tsai, ArrayI1& frac_veg_nosno_alb);

private:
  // read 1 month of data from file (1, npfts, nlat, nlon) for input param month 
  // and place into 2D array arr(ntimes, ncells) where ntimes = arr_idx
  // by combining nlat & nlon and filtering by pft type (vtype) - only one pft per grid cell currently
  void read_month(const Utils::DomainDecomposition<2> &dd, const std::string& filename, const std::string& varname,
    const size_t month, const ArrayI1& vtype, const int arr_idx, ArrayD2& arr);

  // read 3 months of data into member arrays
  void read_initial(const Utils::DomainDecomposition<2> &dd, const std::string& filename,
    const Utils::Date& model_time, const ArrayI1& vtype);

  // advance index - move member data at month idx  1 & 2 to month 0 & 1 and then read new data into month idx 2 
  void read_new_month(const Utils::DomainDecomposition<2> &dd, const std::string& filename,
    const Utils::Date& model_time, const ArrayI1& vtype);

  ArrayD2 mlai_, msai_, mhtop_, mhbot_;
  size_t ncells_, npfts_;
  bool initialized_, need_new_data_;
  std::map<std::string, ArrayD2> phenology_names_;
  int data_m1_;
};

} // namespace ELM::phenology_data

#include "phenology_data_impl.hh"
