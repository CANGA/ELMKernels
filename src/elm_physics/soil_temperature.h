
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


namespace ELM::soil_temperature::detail {

  ACCELERATE
  double calc_lwrad_emit(const double& emg, const double& t_grnd);
  
  ACCELERATE
  double calc_dlwrad_emit(const double& emg, const double& t_grnd);

} //namespace soil_temperature::detail


namespace ELM::soil_temperature {

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


/*
cv[nlevsno+nlevgrnd]     heat capacity [J/(m2 K)]
tk[nlevsno+nlevgrnd]     thermal conductivity [W/(m K)]
fn[nlevsno+nlevgrnd]     heat diffusion through the layer interface [W/m2]
fact[nlevsno+nlevgrnd]   used in computing tridiagonal matrix

*/
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


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  void solve_temperature(const ArrayI1 snl,
                         const ArrayD1 t_h2osfc,
                         const ArrayD2 t_soisno,
                         const ArrayD3 lhs_matrix,
                         ArrayD2 rhs_vector);


} // namespace ELM::soil_temperature

#include "soil_temperature_impl.hh"
