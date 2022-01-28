
#pragma once

#include "monthly_data.h"
#include "elm_constants.h"

#include "array.hh"
#include "utils.hh"
#include "read_input.hh"

#include <utility>
#include <tuple>
#include <cmath>
#include <string>
#include <array>

/*
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
namespace ELM::aerosol_data {

struct AerosolMasses {
using ArrayD2 = ELM::Array<double, 2>;
  ArrayD2 mss_bcphi, mss_bcpho, mss_dst1, mss_dst2, mss_dst3, mss_dst4;
  AerosolMasses(const int ncells);
};


struct AerosolConcentrations {
using ArrayD2 = ELM::Array<double, 2>;
  ArrayD2 mss_cnc_bcphi, mss_cnc_bcpho, mss_cnc_dst1, mss_cnc_dst2, mss_cnc_dst3, mss_cnc_dst4;
  AerosolConcentrations(const int ncells);
};


// functor to calculate aerosol deposition
// this is from AerosolMod.F90
// should be run after inter-layer aerosol fluxes are accounted for in 
// SnowHydrologyMod.F90::SnowWater()
template<typename T, typename ArrayI1>
struct ComputeAerosolDeposition {
  ComputeAerosolDeposition(const T& aerosol_forc, const ArrayI1& snl, AerosolMasses& aerosol_masses);

  void operator()(const int i) const;

private:
  ArrayI1 snl_;
  AerosolMasses aerosol_masses_;
  double forc_bcphi_, forc_bcpho_, forc_dst1_, forc_dst2_, forc_dst3_, forc_dst4_;
};


// functor to calculate aerosol mass and concentration in snow layers
// this is from AerosolMod.F90
// gets run in ELM directly after hydrology is called
template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
struct ComputeAerosolConcenAndMass {
  ComputeAerosolConcenAndMass(const bool& do_capsnow, const double& dtime, const ArrayI1& snl, const ArrayD2& h2osoi_liq,
    const ArrayD2& h2osoi_ice, const ArrayD2& snw_rds, const ArrayD1& qflx_snwcp_ice, AerosolMasses& aerosol_masses,
    AerosolConcentrations& aerosol_concentrations);

  void operator()(const int i) const;

private:
  bool do_capsnow_;
  double dtime_;
  ArrayI1 snl_;
  ArrayD2 h2osoi_liq_;
  ArrayD2 h2osoi_ice_;
  ArrayD2 snw_rds_;
  ArrayD1 qflx_snwcp_ice_;
  AerosolMasses aerosol_masses_;
  AerosolConcentrations aerosol_concentrations_;
};

// convenience function to invoke aerosol mass and concen functor
template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
void invoke_aerosol_concen_and_mass(const bool& do_capsnow, const double& dtime, const ArrayI1& snl, const ArrayD2& h2osoi_liq,
  const ArrayD2& h2osoi_ice, const ArrayD2& snw_rds, const ArrayD1& qflx_snwcp_ice, AerosolMasses& aerosol_masses,
  AerosolConcentrations& aerosol_concentrations);

// class to manage aerosol data
class AerosolDataManager : public MonthlyDataManager {

using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;

public:

  AerosolDataManager();

  // read a slice from file and reshape into 1D array
  void read_variable_slice(const Comm_type &comm, const std::string& filename, const std::string& varname,
    const size_t& lon_idx, const size_t& lat_idx, ArrayD1& arr);

  // read 12 months of values at closest point to lon_d and lat_d
  void read_data(const Comm_type &comm, const std::string& filename, const double& lon_d, const double& lat_d);

  // get closest indices to [lon, lat] 
  std::pair<size_t, size_t> get_nearest_indices(const Comm_type &comm, const std::string& filename,
    const double& lon_d, const double& lat_d);

  // interpolate and accumulate aerosol forcing data to get aerosol sources for this timestep
  auto get_aerosol_source(const Utils::Date& model_time, const double& dtime);

  // convenience function to invoke aerosol deposition source functor
  void invoke_aerosol_source(const Utils::Date& model_time, const double& dtime, const ArrayI1& snl,
    AerosolMasses& aerosol_masses);

protected:
  ArrayD1 bcdep_, bcpho_, bcphi_, dst1_1_, dst1_2_, dst2_1_; 
  ArrayD1 dst2_2_, dst3_1_, dst3_2_, dst4_1_, dst4_2_;

};

} // namespace ELM::aerosol_data

#include "aerosol_data_impl.hh"
