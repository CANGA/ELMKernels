
/*
out:
hs_soil                   ! heat flux on soil [W/m2]
hs_h2osfc                 ! heat flux on standing water [W/m2]
dhsdT                     ! temperature derivative of "hs" [col]

hs_top                    ! net energy flux into surface layer (col) [W/m2]
hs_top_snow               ! heat flux on top snow layer [W/m2]
sabg_chk                  ! sum of soil/snow using current fsno, for balance check [W/m2]

*/

#pragma once

#include "elm_constants.h"

#include "kokkos_includes.hh"

namespace ELM::soil_temp {

/*
generic function to calc surface heat fluxes
can calc:
hs_soil: sabg_soil, t_soisno(nlevsno), eflx_sh_soil, qflx_ev_soil
hs_h2osfc: sabg_soil, t_h2osfc, eflx_sh_h2osfc, qflx_ev_h2osfc
hs_snow: sabg_snow, t_soisno(nlevsno - snl), eflx_sh_snow, qflx_ev_snow
hs_top: sabg_lyr(nlevsno - snl), t_grnd, eflx_sh_grnd, qflx_evap_soi
hs_top_snow: sabg_lyr(nlevsno - snl), t_soisno(nlevsno - snl), eflx_sh_snow, qflx_ev_snow
*/
  ACCELERATE
  double calc_surface_heat_flux(const int& frac_veg_nosno,
                                const double& dlrad,
                                const double& emg,
                                const double& forc_lwrad,
                                const double& htvp,
                                const double& solar_abg,
                                const double& temp,
                                const double& eflx_sh,
                                const double& qflx_ev);

  // calculate derivative of heat flux wrt temperature
  ACCELERATE
  double calc_dhsdT(const double& cgrnd, const double& emg, const double& t_grnd);

  // save sabg for balancecheck, in case frac_sno is set to zero later
  ACCELERATE
  double check_absorbed_solar(const double& frac_sno_eff, const double& sabg_snow, const double& sabg_soil);


  // calculates fn, the diffusive heat flux through layer interfaces
  template <typename ArrayD1>
  ACCELERATE
  void calc_diffusive_heat_flux(const int& snl,
                                const ArrayD1 tk,
                                const ArrayD1 t_soisno,
                                const ArrayD1 z,
                                ArrayD1 fn);


  template <typename ArrayD3, typename ArrayI1, typename ArrayD1, typename ArrayD2>
  void solve_temperature(const double& dtime,
                         const ArrayI1 snl,
                         const ArrayI1 frac_veg_nosno,
                         const ArrayD1 dlrad,
                         const ArrayD1 emg,
                         const ArrayD1 forc_lwrad,
                         const ArrayD1 htvp,
                         const ArrayD1 cgrnd,
                         const ArrayD1 eflx_sh_soil,
                         const ArrayD1 qflx_ev_soil,
                         const ArrayD1 eflx_sh_h2osfc,
                         const ArrayD1 qflx_ev_h2osfc,
                         const ArrayD1 eflx_sh_grnd,
                         const ArrayD1 qflx_evap_soi,
                         const ArrayD1 eflx_sh_snow,
                         const ArrayD1 qflx_ev_snow,
                         const ArrayD1 frac_sno_eff,
                         const ArrayD1 frac_sno,
                         const ArrayD1 frac_h2osfc,
                         const ArrayD1 h2osno,
                         const ArrayD1 h2osfc,
                         const ArrayD1 sabg_snow,
                         const ArrayD1 sabg_soil,
                         const ArrayD2 sabg_lyr,
                         const ArrayD2 h2osoi_liq,
                         const ArrayD2 h2osoi_ice,
                         const ArrayD2 watsat,
                         const ArrayD2 tkmg,
                         const ArrayD2 tkdry,
                         const ArrayD2 csol,
                         const ArrayD2 dz,
                         const ArrayD2 zsoi,
                         const ArrayD2 zisoi,
                         ArrayD1 t_h2osfc,
                         ArrayD1 t_grnd,
                         ArrayD2 t_soisno,
                         ArrayD2 fact);


} // namespace ELM::soil_temp


namespace ELM::soil_temp::detail {

  static constexpr double cnfac{0.5}; // Crank Nicholson factor between 0 and 1

  ACCELERATE
  double calc_lwrad_emit(const double& emg, const double& t_grnd);

  ACCELERATE
  double calc_dlwrad_emit(const double& emg, const double& t_grnd);


  // cv[nlevsno+nlevgrnd]     heat capacity [J/(m2 K)]
  // tk[nlevsno+nlevgrnd]     thermal conductivity [W/(m K)]
  // fn[nlevsno+nlevgrnd]     heat diffusion through the layer interface [W/m2]
  // fact[nlevsno+nlevgrnd]   used in computing matrix
  template <typename ArrayD1>
  ACCELERATE
  void calc_heat_flux_matrix_factor(const int& snl,
                                    const double& dtime,
                                    const ArrayD1 cv,
                                    const ArrayD1 dz,
                                    const ArrayD1 z,
                                    const ArrayD1 zi,
                                    ArrayD1 fact);

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void set_tvector(const int& c,
                   const ArrayI1 snl,
                   const ArrayD1 t_h2osfc,
                   const ArrayD2 t_soisno,
                   ArrayD2 tvector);


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void update_temperature(const int& c,
                          const ArrayI1 snl,
                          const ArrayD1 frac_h2osfc,
                          const ArrayD2 tvector,
                          ArrayD1 t_h2osfc,
                          ArrayD2 t_soisno);

} //namespace soil_temp::detail

#include "soil_temperature_impl.hh"
