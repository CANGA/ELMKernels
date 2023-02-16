
#pragma once

// forward declarations
namespace ELM::aero_data {

  template<typename ArrayD1>
  struct AerosolFileInput;

  template<typename ArrayD2>
  struct AerosolMasses;

  template<typename ArrayD2>
  struct AerosolConcentrations;

} // namespace ELM::aero_data


namespace ELM::aero_impl {

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

} // namespace ELM::aero_impl

namespace ELM {
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
                                      aero_data::AerosolMasses<ArrayD2>& aero_mass,
                                      aero_data::AerosolConcentrations<ArrayD2>& aero_conc);

} // namespace ELM

#include "aerosol_physics_new_impl.hh"
