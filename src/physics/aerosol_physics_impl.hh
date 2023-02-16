
#pragma once

#include <cassert>
#include <utility>

#include "invoke_kernel.hh"
#include "aerosol_physics_new.h"


template<typename T>
ACCELERATE
T ELM::aero_impl::
get_snow_mass(const int& snow_idx, const int& snotop, const T& h2osoi_ice, const T& h2osoi_liq)
{
  return (snow_idx < snotop) ? 1.e-12 : h2osoi_ice + h2osoi_liq;
}


template<typename T>
ACCELERATE
T ELM::aero_impl::
get_snowcap_scl_fct(const int& snow_idx, const int& snotop,
                    const int& do_capsnow, const T& snowmass,
                    const T& qflx_snwcp_ice, const T& dtime)
{
  return (snow_idx == snotop && do_capsnow) ?
         (snowmass / (snowmass + qflx_snwcp_ice * dtime)) : // snowcap_scl_fct
         (snow_idx < snotop) ?
         0.0 : 1.0; // 0.0 for layers above snow, 1.0 (no change) for layers within snowpack
}


template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
void
ELM::compute_aerosol_deposition(const double& dtime,
                                const ArrayI1& snl,
                                const aero_data::AerosolFileInput<ArrayD1>& aero_input,
                                aero_data::AerosolMasses<ArrayD2>& aero_mass)
{
  const int aero_mss_size = aero_mass.mss_bcphi.extent(0);
  assert( (aero_mss_size == aero_input.bcphi.extent(0)) &&
        "Error: different size views in AerosolFileInput and AerosolMasses");
  assert( (aero_mss_size == snl.extent(0)) &&
    "Error: different size ncols dimension (0) in AerosolMasses and snl (snow layers)");

  auto aero_deposition = ELM_LAMBDA (const int i) {

    if (snl(i) > 0) {
      const int j = ELMdims::nlevsno() - snl(i);
      aero_mass.mss_bcphi(i,j) += (aero_input.bcphi(i) * dtime);
      aero_mass.mss_bcpho(i,j) += ((aero_input.bcpho(i) + aero_input.bcdep(i)) * dtime);
      aero_mass.mss_dst1(i,j)  += ((aero_input.dst1_1(i) + aero_input.dst1_2(i)) * dtime);
      aero_mass.mss_dst2(i,j)  += ((aero_input.dst2_1(i) + aero_input.dst2_2(i)) * dtime);
      aero_mass.mss_dst3(i,j)  += ((aero_input.dst3_1(i) + aero_input.dst3_2(i)) * dtime);
      aero_mass.mss_dst4(i,j)  += ((aero_input.dst4_1(i) + aero_input.dst4_2(i)) * dtime);
    }
  };
  apply_parallel_for(aero_deposition, "aero_deposition", aero_mss_size);
}


template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
void
ELM::update_aerosol_mass_and_concen(const double& dtime, const ArrayI1& snl,
                                    const ArrayI1& do_capsnow, const ArrayD1& qflx_snwcp_ice,
                                    const ArrayD2& h2osoi_ice, const ArrayD2& h2osoi_liq,
                                    aero_data::AerosolMasses<ArrayD2>& aero_mass,
                                    aero_data::AerosolConcentrations<ArrayD2>& aero_conc)
{
  const int aero_mass_size = aero_mass.mss_bcphi.extent(0);
  const int aero_conc_size = aero_conc.cnc_bcphi.extent(0);
  assert( (aero_mass_size == aero_conc_size) &&
        "Error: different size views in AerosolConcentrations and AerosolMasses");
  assert( (aero_mass_size == snl.extent(0)) &&
    "Error: different size ncols dimension (0) in AerosolMasses and snl (snow layers)");

  auto aero_mass_concen = ELM_LAMBDA (const int i) {
    // top snow layer index
    const int snotop = ELMdims::nlevsno() - snl(i);
    // loop through snow layers
    for (int sl = 0; sl < ELMdims::nlevsno(); sl++) {

      // get mass of water in snow layer and snowcap scaling factor
      const double snowmass = aero_impl::get_snow_mass(sl, snotop, h2osoi_ice(i,sl), h2osoi_liq(i,sl));
      const double snowcap_scl_fct = aero_impl::get_snowcap_scl_fct(sl, snotop, do_capsnow(i),
                                                         snowmass, qflx_snwcp_ice(i), dtime);
      // update aerosol masses
      aero_mass.mss_bcphi(i,sl) *= snowcap_scl_fct;
      aero_mass.mss_bcpho(i,sl) *= snowcap_scl_fct;
      aero_mass.mss_dst1(i,sl)  *= snowcap_scl_fct;
      aero_mass.mss_dst2(i,sl)  *= snowcap_scl_fct;
      aero_mass.mss_dst3(i,sl)  *= snowcap_scl_fct;
      aero_mass.mss_dst4(i,sl)  *= snowcap_scl_fct;

      // update aerosol concentration
      const double snwmss_inv = 1.0 / snowmass;
      aero_conc.cnc_bcphi(i,sl) = aero_mass.mss_bcphi(i,sl) * snwmss_inv;
      aero_conc.cnc_bcpho(i,sl) = aero_mass.mss_bcpho(i,sl) * snwmss_inv;
      aero_conc.cnc_dst1(i,sl)  = aero_mass.mss_dst1(i,sl)  * snwmss_inv;
      aero_conc.cnc_dst2(i,sl)  = aero_mass.mss_dst2(i,sl)  * snwmss_inv;
      aero_conc.cnc_dst3(i,sl)  = aero_mass.mss_dst3(i,sl)  * snwmss_inv;
      aero_conc.cnc_dst4(i,sl)  = aero_mass.mss_dst4(i,sl)  * snwmss_inv;
    }
  };
  apply_parallel_for(aero_mass_concen, "update_aerosol_mass_and_concen", aero_mass_size);
}
