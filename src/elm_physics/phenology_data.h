
#pragma once

#include "array.hh"
#include "date_time.hh"
#include "elm_constants.h"
#include "utils.hh"

#include "monthly_data.h"
#include "phenology_physics.h"
#include "read_input.hh"

#include <array>
#include <map>
#include <string>

#include "kokkos_includes.hh"
#include "invoke_kernel.hh"

// this is derived from SatellitePhenologyMod.F90
namespace ELM {

// class to manage phenology data
template <typename ArrayD2>
class PhenologyDataManager {

public:

  // these are public to provide access from driver
  // for eg hostview creation
  ArrayD2 mlai, msai, mhtop, mhbot;

  PhenologyDataManager(const Utils::DomainDecomposition<2>& dd, const size_t ncells, const size_t npfts);

  // read data from file
  // either all three months of data, or new single month
  template <typename ArrayI1, typename h_ArrayD2>
  bool read_data(std::map<std::string, h_ArrayD2>& phenology_views, const std::string& filename, const Utils::Date& model_time,
                 const ArrayI1& vtype);

  // get phenology data - call parallel physics kernel - return phenology data for this timestep
  template <typename ArrayI1, typename ArrayD1>
  void get_data(const Utils::Date& model_time, const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype,
                ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot, ArrayD1& tlai, ArrayD1& tsai,
                ArrayI1& frac_veg_nosno_alb);

  // will data be read if read_data is called?
  inline bool need_data() { return need_new_data_; }

private:
  // read 1 month of data from file (1, npfts, nlat, nlon) for input param month
  // and place into 2D array arr(ntimes, ncells) where ntimes = arr_idx
  // by combining nlat & nlon and filtering by pft type (vtype) - only one pft per grid cell currently
  template <typename ArrayI1, typename h_ArrayD2>
  void read_month(const std::string& filename, const std::string& varname,
                  const size_t month, const ArrayI1& vtype, const int arr_idx, h_ArrayD2& arr);

  // read 3 months of data into member arrays
  template <typename ArrayI1, typename h_ArrayD2>
  void read_initial(std::map<std::string, h_ArrayD2>& phenology_views, const std::string& filename, const Utils::Date& model_time,
                    const ArrayI1& vtype);

  // advance index - move member data at month idx  1 & 2 to month 0 & 1 and then read new data into month idx 2
  template <typename ArrayI1, typename h_ArrayD2>
  void read_new_month(std::map<std::string, h_ArrayD2>& phenology_views, const std::string& filename,
                      const Utils::Date& model_time, const ArrayI1& vtype);

  const Utils::DomainDecomposition<2> dd_;
  size_t ncells_, npfts_;
  bool initialized_;
  int data_m1_;
  bool need_new_data_;

};

} // namespace ELM

#include "phenology_data_impl.hh"
