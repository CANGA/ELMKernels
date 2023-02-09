
#pragma once

namespace ELM {

  namespace aero_data {

    template<typename ArrayD1>
    struct AerosolFileInput;

    template<typename ArrayD2>
    struct AerosolMasses;

    template<typename ArrayD2>
    struct AerosolConcentrations;

  } // namespace ELM::aero_data



  // computes aerosol forcing rate for current timestep
  // input data in AerosolFileInput has already been interpolated in time and space
  template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
  void compute_aerosol_deposition(const double& dtime,
                                  const ArrayI1& snl,
                                  const aero_data::AerosolFileInput<ArrayD1>& aero_input,
                                  aero_data::AerosolMasses<ArrayD2>& aero_mass);
  
  // updates aerosol per-layer mass and concentration
  // due to atmospheric deposition and changes in snowpack
  template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
  void update_aerosol_mass_and_concen(const double& dtime, const ArrayI1& snl,
                                      const ArrayI1& do_capsnow, const ArrayD1& qflx_snwcp_ice,
                                      const ArrayD2& h2osoi_ice, const ArrayD2& h2osoi_liq,
                                      ArrayD2& snw_rds,
                                      aero_data::AerosolMasses<ArrayD2>& aero_mass,
                                      aero_data::AerosolConcentrations<ArrayD2>& aero_conc);

} // namespace ELM

namespace ELM::aero_impl {
  
  // variadic summation for aerosols
  template<typename... Args>
  ACCELERATE
  auto aerosol_sum(const Args&&... args);
  
  // aerosol mass deposited during timestep
  template<typename T, typename... Args>
  ACCELERATE
  T aero_mass_deposited_during_dt(const T& dtime, const Args&&... args);
  
  // add aerosol deposition to top layer of snowpack
  template<typename T>
  ACCELERATE
  void add_aero_deposition_to_snow(const T& dtime, const T& bcphi_in,
                                   const T& bcpho_in, const T& bcdep_in,
                                   const T& dst1_1_in, const T& dst1_2_in,
                                   const T& dst2_1_in, const T& dst2_2_in,
                                   const T& dst3_1_in, const T& dst3_2_in,
                                   const T& dst4_1_in, const T& dst4_2_in,
                                   T& mss_bcphi, T& mss_bcpho, T& mss_dst1,
                                   T& mss_dst2, T& mss_dst3, T& mss_dst4);
  
  // get mass of water in snow layer
  template<typename T>
  ACCELERATE
  T get_snow_mass(const int& snow_idx, const int& snotop, const T& h2osoi_ice, const T& h2osoi_liq);
  
  // get snowcap scaling factor
  // this approach conserves the aerosol mass concentration
  // (but not the aerosol mass) when snow-capping is invoked
  template<typename T>
  ACCELERATE
  T get_snowcap_scl_fct(const int& snow_idx, const int& snotop,
                        const int& do_capsnow, const T& snowmass,
                        const T& qflx_snwcp_ice, const T& dtime);
  
  // update aerosol layer mass
  template<typename T>
  ACCELERATE
  void update_aerosol_mass(const T& snowcap_scl_fct, T& mss_bcphi, T& mss_bcpho,
                           T& mss_dst1, T& mss_dst2, T& mss_dst3, T& mss_dst4);
  
  // update aerosol layer concentration
  template<typename T>
  ACCELERATE
  void update_aerosol_concen(const T& snowmass, const T& mss_bcphi, const T& mss_bcpho,
                             const T& mss_dst1, const T& mss_dst2, const T& mss_dst3,
                             const T& mss_dst4, T& cnc_bcphi, T& cnc_bcpho, T& cnc_dst1,
                             T& cnc_dst2, T& cnc_dst3, T& cnc_dst4);

  // this shouldn't be here
  // move during refactor
  // reset snw_rds(j) to 0.0 if no snow in layer j
  template<typename T>
  ACCELERATE
  void reset_snow_rds(const int& snow_idx, const int& snotop, T& snw_rds);

} // namespace ELM::aero_impl

#include "aerosol_physics_new_impl.hh"
