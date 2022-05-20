
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
  
  // t_soisno is t_soisno(nlevsno-snl)
  ACCELERATE
  double calc_lwrad_emit_snow(const double& emg, const double& t_soisno);
  
  // t_soisno is t_soisno(nlevsno)
  ACCELERATE
  double calc_lwrad_emit_soil(const double& emg, const double& t_soisno);
  
  ACCELERATE
  double calc_lwrad_emit_h2osfc(const double& emg, const double& t_h2osfc);

} //namespace soil_temperature::detail


namespace ELM::soil_temperature {

  // t_soisno is t_soisno(nlevsno)
  // return value is also ELM's elfx_gnet_soil without any pft weighting
  ACCELERATE
  double calc_hs_soil(const int& frac_veg_nosno,
                      const double& t_soisno,
                      const double& sabg_soil,
                      const double& dlrad,
                      const double& emg,
                      const double& forc_lwrad,
                      const double& eflx_sh_soil,
                      const double& qflx_ev_soil,
                      const double& htvp);

  // return value is also ELM's elfx_gnet_h2osfc without any pft weighting
  ACCELERATE
  double calc_hs_h2osfc(const int& frac_veg_nosno,
                        const double& t_h2osfc,
                        const double& sabg_soil,
                        const double& dlrad,
                        const double& emg,
                        const double& forc_lwrad,
                        const double& eflx_sh_soil,
                        const double& qflx_ev_soil,
                        const double& htvp);

  // sabg_lyr is sabg_lyr(nlevsno-snl)
  ACCELERATE
  double calc_hs_top(const int& frac_veg_nosno,
                     const double& t_grnd,
                     const double& sabg_lyr,
                     const double& dlrad,
                     const double& emg,
                     const double& forc_lwrad,
                     const double& eflx_sh_grnd,
                     const double& qflx_evap_soi,
                     const double& htvp);

  // t_soisno is t_soisno(nlevsno-snl)
  // sabg_lyr is sabg_lyr(nlevsno-snl)
  // return value is also ELM's elfx_gnet_snow without any pft weighting
  ACCELERATE
  double calc_hs_top_snow(const int& frac_veg_nosno,
                          const double& t_soisno,
                          const double& sabg_lyr,
                          const double& dlrad,
                          const double& emg,
                          const double& forc_lwrad,
                          const double& eflx_sh_snow,
                          const double& qflx_ev_snow,
                          const double& htvp);

  ACCELERATE
  double calc_dhsdT (const double& cgrnd, const double& emg, const double& t_grnd);

  ACCELERATE
  double check_absorbed_solar(const double& frac_sno_eff, const double& sabg_snow, const double& sabg_soil);


  template <typename ArrayI1, typename ArrayD1>
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
  template <typename ArrayI1, typename ArrayD1>
  ACCELERATE
  void calc_heat_flux_matrix_factor(const int& snl,
                                    const double& dtime,
                                    const ArrayD1 cv,
                                    const ArrayD1 dz,
                                    const ArrayD1 z,
                                    const ArrayD1 zi,
                                    ArrayD1 fact);

} // namespace ELM::soil_temperature

#include "soil_temperature_impl.hh"
