#ifndef ELM_CANOPY_HYDROLOGY_INTERFACE_HH_
#define ELM_CANOPY_HYDROLOGY_INTERFACE_HH_

namespace ELM {

void CanopyHydrologyKern1(double dtime,
                          const double& forc_rain, const double& forc_snow, const double& irrig_rate,
                          const int& ltype, const int& ctype,
                          const bool& urbpoi, const bool& do_capsnow,
                          const double& elai, const double& esai,
                          const double& dewmx, const double& frac_veg_nosno,
                          double& h2ocan,
                          int& n_irrig_steps_left,
                          double& qflx_prec_intr,
                          double& qflx_irrig,
                          double& qflx_prec_grnd,
                          double& qflx_snwcp_liq,
                          double& qflx_snwcp_ice,
                          double& qflx_snow_grnd_patch,
                          double& qflx_rain_grnd);
                            
  
} // namespace

#endif
