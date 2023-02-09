
#pragma once

#include <cassert>
#include <utility>

#include "invoke_kernel.hh"
#include "aerosol_physics_new.h"

template<typename... Args>
ACCELERATE
auto ELM::aero_impl::
aerosol_sum(const Args&... args)
{
  return (args + ...);
}

template<typename T, typename... Args>
ACCELERATE
T ELM::aero_impl::
aero_mass_deposited_during_dt(const T& dtime, Args&&... args)
{
  return dtime * aerosol_sum(std::forward<Args>(args)...);
}

template<typename T>
ACCELERATE
void ELM::aero_impl::
add_aero_deposition_to_snow(const T& dtime, const T& bcphi_in,
                            const T& bcpho_in, const T& bcdep_in,
                            const T& dst1_1_in, const T& dst1_2_in,
                            const T& dst2_1_in, const T& dst2_2_in,
                            const T& dst3_1_in, const T& dst3_2_in,
                            const T& dst4_1_in, const T& dst4_2_in,
                            T& mss_bcphi, T& mss_bcpho, T& mss_dst1,
                            T& mss_dst2, T& mss_dst3, T& mss_dst4)
{
  mss_bcphi += aero_mass_deposited_during_dt(dtime, bcphi_in);
  mss_bcpho += aero_mass_deposited_during_dt(dtime, bcpho_in, bcdep_in);
  mss_dst1  += aero_mass_deposited_during_dt(dtime, dst1_1_in, dst1_2_in);
  mss_dst2  += aero_mass_deposited_during_dt(dtime, dst2_1_in, dst2_2_in);
  mss_dst3  += aero_mass_deposited_during_dt(dtime, dst3_1_in, dst3_2_in);
  mss_dst4  += aero_mass_deposited_during_dt(dtime, dst4_1_in, dst4_2_in);
}

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

// this shouldn't be here
// move during refactor
template<typename T>
ACCELERATE
void ELM::aero_impl::
reset_snow_rds(const int& snow_idx, const int& snotop, T& snw_rds)
{
  if (snow_idx < snotop) {
    snw_rds = 0.0;
  }
}

template<typename T>
ACCELERATE
void ELM::aero_impl::
update_aerosol_mass(const T& snowcap_scl_fct,
                    T& mss_bcphi, T& mss_bcpho,
                    T& mss_dst1, T& mss_dst2,
                    T& mss_dst3, T& mss_dst4)
{
  mss_bcphi *= snowcap_scl_fct;
  mss_bcpho *= snowcap_scl_fct;
  mss_dst1  *= snowcap_scl_fct;
  mss_dst2  *= snowcap_scl_fct;
  mss_dst3  *= snowcap_scl_fct;
  mss_dst4  *= snowcap_scl_fct;
}

template<typename T>
ACCELERATE
void ELM::aero_impl::
update_aerosol_concen(const T& snowmass, const T& mss_bcphi, const T& mss_bcpho,
                      const T& mss_dst1, const T& mss_dst2, const T& mss_dst3,
                      const T& mss_dst4, T& cnc_bcphi, T& cnc_bcpho, T& cnc_dst1,
                      T& cnc_dst2, T& cnc_dst3, T& cnc_dst4)
{
  const double snwmss_inv = 1.0 / snowmass;
  cnc_bcphi = mss_bcphi * snwmss_inv;
  cnc_bcpho = mss_bcpho * snwmss_inv;
  cnc_dst1  = mss_dst1  * snwmss_inv;
  cnc_dst2  = mss_dst2  * snwmss_inv;
  cnc_dst3  = mss_dst3  * snwmss_inv;
  cnc_dst4  = mss_dst4  * snwmss_inv;
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
      aero_impl::add_aero_deposition_to_snow(dtime,
                                  aero_input.bcphi(i),  aero_input.bcpho(i),
                                  aero_input.bcdep(i),  aero_input.dst1_1(i),
                                  aero_input.dst1_2(i), aero_input.dst2_1(i),
                                  aero_input.dst2_2(i), aero_input.dst3_1(i),
                                  aero_input.dst3_2(i), aero_input.dst4_1(i),
                                  aero_input.dst4_2(i),
                                  aero_mass.mss_bcphi(i,j), aero_mass.mss_bcpho(i,j),
                                  aero_mass.mss_dst1(i,j), aero_mass.mss_dst2(i,j),
                                  aero_mass.mss_dst3(i,j), aero_mass.mss_dst4(i,j));

    }
  };
  apply_parallel_for(aero_deposition, "aero_deposition", aero_mss_size);
}


template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
void
ELM::update_aerosol_mass_and_concen(const double& dtime, const ArrayI1& snl,
                                    const ArrayI1& do_capsnow, const ArrayD1& qflx_snwcp_ice,
                                    const ArrayD2& h2osoi_ice, const ArrayD2& h2osoi_liq,
                                    ArrayD2& snw_rds,
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
    for (int sl = 0; sl < ELMdims::nlevsno(); ++sl) {
      // get mass of water in snow layer
      const double snowmass = get_snow_mass(sl, snotop, h2osoi_ice(i,sl), h2osoi_liq(i,sl));
      // get snowcap scaling factor
      const double snowcap_scl_fct = get_snowcap_scl_fct(sl, snotop, do_capsnow(i),
                                                         snowmass, qflx_snwcp_ice(i), dtime);
      // update aerosol masses
      aero_impl::update_aerosol_mass(snowcap_scl_fct,
                          aero_mass.mss_bcphi(i,sl), aero_mass.mss_bcpho(i,sl),
                          aero_mass.mss_dst1(i,sl),  aero_mass.mss_dst2(i,sl),
                          aero_mass.mss_dst3(i,sl),  aero_mass.mss_dst4(i,sl));
      // update aerosol concentration
      aero_impl::update_aerosol_concen(snowmass,
                          aero_mass.mss_bcphi(i,sl), aero_mass.mss_bcpho(i,sl),
                          aero_mass.mss_dst1(i,sl),  aero_mass.mss_dst2(i,sl),
                          aero_mass.mss_dst3(i,sl),  aero_mass.mss_dst4(i,sl),
                          aero_conc.cnc_bcphi(i,sl), aero_conc.cnc_bcpho(i,sl),
                          aero_conc.cnc_dst1(i,sl),  aero_conc.cnc_dst2(i,sl),
                          aero_conc.cnc_dst3(i,sl),  aero_conc.cnc_dst4(i,sl));
      // TODO move this somewhere else
      // this is where it gets invoked in ELM
      reset_snow_rds(sl, snotop, snw_rds(i,sl));
    }
  };
  apply_parallel_for(aero_mass_concen, "update_aerosol_mass_and_concen", aero_mass_size);
}
