#ifndef ELM_CANOPY_HYDROLOGY_INTERFACE_HH_
#define ELM_CANOPY_HYDROLOGY_INTERFACE_HH_

#include "CanopyHydrology_private.h"

namespace ELM {

inline void CanopyHydrology_Interception(double dtime,
        const double& forc_rain,
        const double& forc_snow,
        const double& irrig_rate,
        const int& ltype, const int& ctype,
        const bool& urbpoi, const bool& do_capsnow,
        const double& elai, const double& esai,
        const double& dewmx, const int& frac_veg_nosno,
        double& h2ocan,
        int& n_irrig_steps_left,
        double& qflx_prec_intr,
        double& qflx_irrig,
        double& qflx_prec_grnd,
        double& qflx_snwcp_liq,
        double& qflx_snwcp_ice,
        double& qflx_snow_grnd_patch,
        double& qflx_rain_grnd) {
  canopyhydrology_interception_(&dtime, &forc_rain, &forc_snow, &irrig_rate,
                        &ltype, &ctype, &urbpoi, &do_capsnow, &elai, &esai,
                        &dewmx, &frac_veg_nosno, &h2ocan, &n_irrig_steps_left,
                        &qflx_prec_intr, &qflx_irrig, &qflx_prec_grnd,
                        &qflx_snwcp_liq, &qflx_snwcp_ice, & qflx_snow_grnd_patch, &qflx_rain_grnd);
}



inline void CanopyHydrology_FracWet(const int& frac_veg_nosno,
        const double& h2ocan,
        const double& elai, 
        const double& esai,
        const double& dewmx,
        double& fwet,
        double& fdry) {
  canopyhydrology_fracwet_(&frac_veg_nosno, &h2ocan, &elai, &esai, &dewmx, &fwet, &fdry);
}


template<typename Array_d>
void CanopyHydrology_SnowWater(const double& dtime,
        const double& qflx_floodg,
        const int& ltype,
        const int& ctype,
        const bool& urbpoi,
        const bool& do_capsnow,                            
        const int& oldfflag,
        const double& forc_t,
        const double& t_grnd,
        const double& qflx_snow_grnd_col,
        const double& qflx_snow_melt,
        const double& n_melt,
        const double& frac_h2osfc,
        double& snow_depth,
        double& h2osno,
        double& int_snow,
        Array_d& swe_old,
        Array_d& h2osoi_liq,
        Array_d& h2osoi_ice,
        Array_d& t_soisno,
        Array_d& frac_iceold,
        int& snl,
        Array_d& dz,
        Array_d& z,
        Array_d& zi,
        int& newnode,
        double& qflx_floodc,
        double& qflx_snow_h2osfc,
        double& frac_sno_eff,
        double& frac_sno) {
  canopyhydrology_snowwater_(&dtime, &qflx_floodg, &ltype, &ctype, &urbpoi,
                             &do_capsnow, &oldfflag, &forc_t, &t_grnd,
                             &qflx_snow_grnd_col, &qflx_snow_melt, &n_melt,
                             &frac_h2osfc, &snow_depth, &h2osno, &int_snow,
                             &swe_old[0], &h2osoi_liq[0], &h2osoi_ice[0], &t_soisno[0],
                             &frac_iceold[0], &snl, &dz[0], &z[0], &zi[0], &newnode,
                             &qflx_floodc, &qflx_snow_h2osfc, &frac_sno_eff, &frac_sno);
}


inline void CanopyHydrology_FracH2OSfc(const double& dtime,
        const double& min_h2osfc,
        const int& ltype,
        const double& micro_sigma,
        const double& h2osno,
        double& h2osfc,
        double& h2osoi_liq,
        double& frac_sno,
        double& frac_sno_eff,
        double& qflx_h2osfc2topsoi,
        double& frac_h2osfc,
        bool no_update=false) {
  canopyhydrology_frach2osfc_(&dtime, &min_h2osfc, &ltype, &micro_sigma, &h2osno,
          &h2osfc, &h2osoi_liq, &frac_sno, &frac_sno_eff, &qflx_h2osfc2topsoi,
          &frac_h2osfc, &no_update);
}

  
} // namespace

#endif
