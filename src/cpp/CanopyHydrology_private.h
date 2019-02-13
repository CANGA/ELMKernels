#ifndef ELM_CANOPY_HYDROLOGY_INTERFACE_PRIVATE_HH_
#define ELM_CANOPY_HYDROLOGY_INTERFACE_PRIVATE_HH_


extern "C" {
  void canopyhydrology_interception_(const double* dtime,
                            const double* forc_rain,
                            const double* forc_snow,
                            const double* irrig_rate,
                            const int* ltype,
                            const int* ctype,
                            const bool* urbpoi,
                            const bool* do_capsnow,                            
                            const double* elai, 
                            const double* esai,
                            const double* dewmx,
                            const int* frac_veg_nosno,
                            const double* h2ocan,
                            const int* n_irrig_steps_left,
                            const double* qflx_prec_intr,
                            const double* qflx_irrig,
                            const double* qflx_prec_grnd,
                            const double* qflx_snwcp_liq,
                            const double* qflx_snwcp_ice,
                            const double* qflx_snow_grnd_patch,
                            const double* qflx_rain_grnd);

                            
  
}

#endif
