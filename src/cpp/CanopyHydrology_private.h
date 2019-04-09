#ifndef ELM_CANOPY_HYDROLOGY_FORT_PRIVATE_HH_
#define ELM_CANOPY_HYDROLOGY_FORT_PRIVATE_HH_


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

  void canopyhydrology_fracwet_(const int* frac_veg_nosno,
          const double* h2ocan,
          const double* elai, 
          const double* esai,
          const double* dewmx,
          const double* fwet,
          const double* fdry);

  void canopyhydrology_snowwater_(const double* dtime,
          const double* qflx_floodg,
          const int* ltype,
          const int* ctype,
          const bool* urbpoi,
          const bool* do_capsnow,                            
          const int* oldfflag,
          const double* forc_t,
          const double* t_grnd,
          const double* qflx_snow_grnd_col,
          const double* qflx_snow_melt,
          const double* n_melt,
          const double* frac_h2osfc,
          const double* snow_depth,
          const double* h2osno,
          const double* int_snow,
          const double* swe_old,
          const double* h2osoi_liq,
          const double* h2osoi_ice,
          const double* t_soisno,
          const double* frac_iceold,
          const int* snl,
          const double* dz,
          const double* z,
          const double* zi,
          const int* newnode,
          const double* qflx_floodc,
          const double* qflx_snow_h2osfc,
          const double* frac_sno_eff,
          const double* frac_sno);

  void canopyhydrology_frach2osfc_(const double* dtime,
          const double* min_h2osfc,
          const int* ltype,
          const double* micro_sigma,
          const double* h2osno,
          const double* h2osfc,
          const double* h2osoi_liq,
          const double* frac_sno,
          const double* frac_sno_eff,
          const double* qflx_h2osfc2topsoi,
          const double* frac_h2osfc,
          const bool* no_update);
  
}

#endif
