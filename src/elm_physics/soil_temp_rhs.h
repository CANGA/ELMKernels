

#pragma once

#include "elm_constants.h"

#include "kokkos_includes.hh"

namespace ELM::soil_temp_rhs {

  // rhs_vec(ncells, nlevgrnd+nlevsno+1) RHS vector for numerical solution of temperature
  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  void set_RHS(const double& dtime,
               const ArrayI1 snl,
               const ArrayD1 hs_top_snow,
               const ArrayD1 dhsdT,
               const ArrayD1 hs_soil,
               const ArrayD1 frac_sno_eff,
               const ArrayD2 t_soisno,
               const ArrayD2 fact,
               const ArrayD2 fn,
               const ArrayD2 sabg_lyr,
               const ArrayD2 z,
               const ArrayD1 tk_h2osfc,
               const ArrayD1 t_h2osfc,
               const ArrayD1 dz_h2osfc,
               const ArrayD1 c_h2osfc,
               const ArrayD1 hs_h2osfc,
               ArrayD2 rhs_vec);

} // namespace ELM::soil_temp_rhs


namespace ELM::soil_temp_rhs::detail {

  static constexpr double cnfac{0.5}; // Crank Nicholson factor between 0 and 1
  
  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_rhs_snow(const int& c,
                    const ArrayI1 snl,
                    const ArrayD1 hs_top_snow,
                    const ArrayD1 dhsdT,
                    const ArrayD2 t_soisno,
                    const ArrayD2 fact,
                    const ArrayD2 fn,
                    const ArrayD2 sabg_lyr,
                    ArrayD2 rt_snow);

  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_rhs_ssw(const int& c,
                   const double& dtime,
                   const ArrayD1 tk_h2osfc,
                   const ArrayD1 t_h2osfc,
                   const ArrayD1 dz_h2osfc,
                   const ArrayD1 c_h2osfc,
                   const ArrayD1 hs_h2osfc,
                   const ArrayD1 dhsdT,
                   const ArrayD2 t_soisno,
                   const ArrayD2 z,
                   ArrayD1 fn_h2osfc,
                   ArrayD1 rt_ssw);

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_rhs_soil(const int& c,
                    const ArrayI1 snl,
                    const ArrayD1 hs_soil,
                    const ArrayD1 hs_top_snow,
                    const ArrayD1 frac_sno_eff,
                    const ArrayD1 dhsdT,
                    const ArrayD2 t_soisno,
                    const ArrayD2 fact,
                    const ArrayD2 fn,
                    const ArrayD2 sabg_lyr,
                    ArrayD2 rt_soil);

  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void assemble_rhs(const int& c,
                    const ArrayD2 rt_snow,
                    const ArrayD1 rt_ssw,
                    const ArrayD2 rt_soil,
                    ArrayD2 rhs_vec);

} // namespace ELM::soil_temp_rhs::detail

#include "soil_temp_rhs_impl.hh"
