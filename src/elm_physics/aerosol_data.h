
#pragma once

#include "monthly_data.h"
#include "read_input.hh"

#include "array.hh"
#include "date_time.hh"
#include "utils.hh"

#include <string>
#include <utility>

#include <array>
#include <cmath>
#include <functional>
#include <tuple>

#include "kokkos_includes.hh"
#include "invoke_kernel.hh"




//#include "array.hh"
//#include "date_time.hh"
//#include "utils.hh"
//
//#include <string>
//#include <utility>
//
//#include "kokkos_includes.hh"

/*
setup to read 12 months of aerosol forcing
it is currently the drivers responsibility 
to keep track of time 
ELM MAPPING:
IDX   VARNAME
1    "BCDEPWET"
2    "BCPHODRY"
3    "BCPHIDRY"
4    "OCDEPWET"  -- DON'T NEED
5    "OCPHODRY"  -- DON'T NEED
6    "OCPHIDRY"  -- DON'T NEED
7    "DSTX01DD"
8    "DSTX02DD"
9    "DSTX03DD"
10   "DSTX04DD"
11   "DSTX01WD"
12   "DSTX02WD"
13   "DSTX03WD"
14   "DSTX04WD"

VARS USED IN THIS MODEL      input dependencies(elm index)
mss_cnc_bcphi                3
mss_cnc_bcpho                1 & 2
mss_cnc_dst1                 7 & 8
mss_cnc_dst2                 9 & 10
mss_cnc_dst3                 11 & 12
mss_cnc_dst4                 13 & 14

ALIAS TO DATA FROM FILE
"BCDEPWET"  bcdep
"BCPHODRY"  bcpho
"BCPHIDRY"  bcphi
"DSTX01DD"  dst1_1
"DSTX02DD"  dst1_2
"DSTX03DD"  dst2_1
"DSTX04DD"  dst2_2
"DSTX01WD"  dst3_1
"DSTX02WD"  dst3_2
"DSTX03WD"  dst4_1
"DSTX04WD"  dst4_2

DEPENDENCY ALIASES
mss_cnc_bcphi - bcphi
mss_cnc_bcpho - bcdep & bcpho
mss_cnc_dst1 - dst1_1 & dst1_2
mss_cnc_dst2 - dst2_1 & dst2_2
mss_cnc_dst3 - dst3_1 & dst3_2
mss_cnc_dst4 - dst4_1 & dst4_2

*/
namespace ELM {


// class to manage aerosol data
template <typename ArrayD1>
class AerosolDataManager {

public:

  // public to allow access from driver
  ArrayD1 bcdep, bcpho, bcphi, dst1_1, dst1_2, dst2_1;
  ArrayD1 dst2_2, dst3_1, dst3_2, dst4_1, dst4_2;

  AerosolDataManager();
  ~AerosolDataManager() = default;

  // interpolate and accumulate aerosol forcing data to get aerosol sources for this timestep
  auto get_aerosol_source(const Utils::Date& model_time, const double& dtime) const; 
};


// read 12 months of values at closest point to lon_d and lat_d
template <typename h_ArrayD1>
void read_aerosol_data(std::map<std::string, h_ArrayD1>& aerosol_views, const Comm_type& comm,
  const std::string& filename, const double& lon_d, const double& lat_d);


namespace aerosol_utils {

  // get closest indices to [lon, lat]
  std::pair<size_t, size_t> get_nearest_indices(const Comm_type& comm, const std::string& filename, const double& lon_d,
                                                const double& lat_d);

  // read a slice from file and reshape into 1D array
  template <typename h_ArrayD1>
  void read_variable_slice(const Comm_type& comm, const std::string& filename, const std::string& varname,
                           const size_t& lon_idx, const size_t& lat_idx, h_ArrayD1 arr);
}

} // namespace ELM

#include "aerosol_data_impl.hh"
