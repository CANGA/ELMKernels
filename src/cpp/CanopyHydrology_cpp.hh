#ifndef ELM_CANOPY_HYDROLOGY_CPP_HH_
#define ELM_CANOPY_HYDROLOGY_CPP_HH_

namespace ELM {

  void CanopyHydrology_Interception(double dtime,
        const double& forc_rain,
        const double& forc_snow,
        const double& irrig_rate,
        const int& ltype, const int& ctype,
        const bool& urbpoi, const bool& do_capsnow,
        const double& elai, const double& esai,
        const double& dewmx, const int& frac_veg_nosno,
        double& h2ocan,
        int n_irrig_steps_left, // need to be fixed
        double& qflx_prec_intr,
        double& qflx_irrig,
        double& qflx_prec_grnd,
        double& qflx_snwcp_liq,
        double& qflx_snwcp_ice,
        double& qflx_snow_grnd_patch,
        double& qflx_rain_grnd) ;



 void CanopyHydrology_FracWet(const int& frac_veg_nosno,
        const double& h2ocan,
        const double& elai, 
        const double& esai,
        const double& dewmx,
        double& fwet,
        double& fdry) ;


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
        double& frac_sno) ;


 void CanopyHydrology_FracH2OSfc(const double& dtime,
        const double& min_h2osfc,
        const int& ltype,
        const double& micro_sigma,
        double& h2osno,
        double& h2osfc,
        double& h2osoi_liq,
        double& frac_sno,
        double& frac_sno_eff,
        double& qflx_h2osfc2topsoi,
        double& frac_h2osfc) ;

  
} // namespace

#endif
