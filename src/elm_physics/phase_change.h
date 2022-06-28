
#pragma once

#include "elm_constants.h"

#include "compile_options.hh"

namespace ELM::soil_temp {

// h2osoi_ice, t_soisno, and fact are single elements of 
// larger column arrays - they are located at the level
// of the first snow layer (nlevsno - 1)
ACCELERATE
void phase_change_h2osfc(const int& snl,
                         const double& dtime,
                         const double& frac_sno,
                         const double& frac_h2osfc,
                         const double& dhsdT,
                         const double& c_h2osfc,
                         const double& fact_sl1,
                         double& t_h2osfc,
                         double& h2osfc,
                         double& xmf_h2osfc,
                         double& qflx_h2osfc_to_ice,
                         double& eflx_h2osfc_to_snow,
                         double& h2osno,
                         double& int_snow,
                         double& snow_depth,
                         double& h2osoi_ice_sl1,
                         double& t_soisno_sl1);


template<typename ArrayI1, typename ArrayD1>
ACCELERATE
void phase_change_soisno(const int& snl,
                         const int& ltype,
                         const double& dtime,
                         const double& dhsdT,
                         const double& frac_h2osfc,
                         const double& frac_sno_eff,
                         const ArrayD1 fact,
                         const ArrayD1 watsat,
                         const ArrayD1 sucsat,
                         const ArrayD1 bsw,
                         const ArrayD1 dz,
                         double& h2osno,
                         double& snow_depth,
                         double& xmf,
                         double& qflx_snofrz,
                         double& qflx_snow_melt,
                         double& qflx_snomelt,
                         double& eflx_snomelt,
                         ArrayI1 imelt,
                         ArrayD1 qflx_snofrz_lyr,
                         ArrayD1 h2osoi_ice,
                         ArrayD1 h2osoi_liq,
                         ArrayD1 t_soisno);


template<typename ArrayD1>
ACCELERATE
void phase_change_correction(const int& snl,
                             const ArrayD1 tk,
                             const ArrayD1 t_soisno,
                             const ArrayD1 z,
                             ArrayD1 fn1);


} // namespace ELM::soil_temp

#include "phase_change_impl.hh"
