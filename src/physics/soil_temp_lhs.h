
#pragma once

#include "elm_constants.h"

#include "compile_options.hh"

namespace ELM::soil_temp {

  // lhs_matrix(ncells, nlevgrnd()+nlevsno()+1, nband()) LHS matrix for numerical solution of temperature
  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  void set_LHS(const double& dtime,
               const ArrayI1 snl,
               const ArrayD1 dz_h2osfc,
               const ArrayD1 c_h2osfc,
               const ArrayD1 tk_h2osfc,
               const ArrayD1 frac_h2osfc,
               const ArrayD1 frac_sno_eff,
               const ArrayD1 dhsdT,
               const ArrayD2 z,
               const ArrayD2 fact,
               const ArrayD2 tk,
               ArrayD3 lhs_matrix);

} // namespace ELM::soil_temp


namespace ELM::soil_temp::detail {

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void get_matrix_snow(const int& c,
                       const ArrayI1 snl,
                       const ArrayD1 dhsdT,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD3 bmatrix_snow);


  template <typename ArrayI1, typename ArrayD2>
  ACCELERATE
  void get_matrix_snow_soil(const int& c,
                       const ArrayI1 snl,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD2 bmatrix_snow_soil);


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void get_matrix_soil(const int& c,
                       const ArrayI1 snl,
                       const ArrayD1 dhsdT,
                       const ArrayD1 frac_sno_eff,
                       const ArrayD1 frac_h2osfc,
                       const ArrayD1 dz_h2osfc,
                       const ArrayD1 tk_h2osfc,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD3 bmatrix_soil);


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_soil_snow(const int& c,
                       const ArrayI1 snl,
                       const ArrayD1 frac_sno_eff,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD2 bmatrix_soil_snow);


  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_ssw(const int& c,
                      const double& dtime,
                      const ArrayD1 dz_h2osfc,
                      const ArrayD1 c_h2osfc,
                      const ArrayD1 tk_h2osfc,
                      const ArrayD1 dhsdT,
                      const ArrayD2 z,
                      ArrayD2 bmatrix_ssw);


  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_ssw_soil(const int& c,
                           const double& dtime,
                           const ArrayD1 dz_h2osfc,
                           const ArrayD1 c_h2osfc,
                           const ArrayD1 tk_h2osfc,
                           const ArrayD2 z,
                           ArrayD2 bmatrix_ssw_soil);


  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_soil_ssw(const int& c,
                           const double& dtime,
                           const ArrayD1 frac_h2osfc,
                           const ArrayD1 dz_h2osfc,
                           const ArrayD1 tk_h2osfc,
                           const ArrayD2 fact,
                           const ArrayD2 z,
                           ArrayD2 bmatrix_soil_ssw);


  template <typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void assemble_lhs(const int& c,
                    const ArrayD2 bmatrix_snow_soil,
                    const ArrayD2 bmatrix_ssw_soil,
                    const ArrayD2 bmatrix_soil_snow,
                    const ArrayD2 bmatrix_soil_ssw,
                    const ArrayD2 bmatrix_ssw,
                    const ArrayD3 bmatrix_snow,
                    const ArrayD3 bmatrix_soil,
                    ArrayD3 lhs_matrix);


} // namespace ELM::soil_temp::detail

#include "soil_temp_lhs_impl.hh"
