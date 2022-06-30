
#pragma once

#include "land_data.h"

namespace ELM {

  template
  <typename ArrayB1, typename ArrayI1, typename ArrayI2, typename ArrayD1,
   typename ArrayD2, typename ArrayD3, typename ArrayPSN1>
  struct ELMState {

      ELMState(const int ncells);
      ~ELMState() = default;

      // forcing data
      ArrayD1 forc_tbot, forc_thbot, forc_pbot, forc_qbot, forc_rh, forc_lwrad,
              forc_rain, forc_snow, forc_u, forc_v, forc_hgt, forc_hgt_u,
              forc_hgt_t, forc_hgt_q, forc_hgt_u_patch, forc_hgt_t_patch,
              forc_hgt_q_patch, forc_vp, forc_rho, forc_po2, forc_pco2;
      ArrayD2 forc_solai, forc_solad;

      // prescribed sat phenology
      ArrayD1 tlai, tsai, elai, esai, htop, hbot;

      // frac veg w and w/o albedo consideration
      ArrayI1 frac_veg_nosno_alb, frac_veg_nosno;

      // soil hydraulics
      ArrayD2 watsat, sucsat, bsw, watdry, watopt, watfc;

      // topo, microtopography
      ArrayD1 n_melt, micro_sigma, topo_slope, topo_std;

      // soil color and texture constants
      //int max_soil_color = 20; // largest option - can't know until NC file is read
      ArrayI1 isoicol;
      ArrayD2 albsat, albdry;

      // soil thermal constants
      ArrayD2 tkmg, tkdry, csol;

      // snow variables
      ArrayI1 snl;
      ArrayD1 snow_depth, frac_sno, int_snow;
      ArrayD2 snw_rds, swe_old;

      // soil/snow/surface h2o state
      ArrayD2 frac_iceold, h2osoi_liq, h2osoi_ice, h2osoi_vol;
      ArrayD1 h2ocan, h2osno, fwet, fdry, h2osfc, frac_h2osfc, frac_sno_eff;

      // water fluxes
      ArrayD1 qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd,
              qflx_rain_grnd, qflx_snow_melt;

      // soil and snow temp
      ArrayD2 t_soisno;
      ArrayD1 t_grnd;

      // for can_sun_shade
      ArrayI1 nrad;
      ArrayD1 laisun, laisha;
      ArrayD2 tlai_z, fsun_z, fabd_sun_z, fabd_sha_z, fabi_sun_z,
              fabi_sha_z, parsun_z, parsha_z, laisun_z, laisha_z;

      // for surface rad
      ArrayD1 sabg_soil, sabg_snow, sabg, sabv, fsa, fsr;
      ArrayD2 sabg_lyr, ftdd, ftid, ftii, fabd, fabi, albsod,
              albsoi, albsnd_hst, albsni_hst, albgrd, albgri,
              flx_absdv, flx_absdn, flx_absiv, flx_absin, albd, albi;

      // variables for CanopyTemperature
      ArrayD1 t_h2osfc, t_h2osfc_bef, soilalpha, soilbeta, qg_snow, qg_soil,
              qg, qg_h2osfc, dqgdT, htvp, emg, emv, z0mg, z0hg, z0qg, z0mv,
              z0hv, z0qv, thv, z0m, displa, thm, eflx_sh_tot, eflx_lh_tot,
              eflx_sh_veg, qflx_evap_tot, qflx_evap_veg, qflx_tran_veg;
      ArrayD2 tssbef;

      // bareground fluxes
      ArrayD1 dlrad, ulrad, eflx_sh_grnd, eflx_sh_snow, eflx_sh_soil, eflx_sh_h2osfc,
              qflx_evap_soi, qflx_ev_snow, qflx_ev_soil, qflx_ev_h2osfc, t_ref2m,
              q_ref2m, rh_ref2m, cgrnds, cgrndl, cgrnd;

      // canopy fluxes
      ArrayI1 altmax_indx, altmax_lastyear_indx;
      ArrayD1 t10, vcmaxcintsha, vcmaxcintsun, btran, t_veg;
      ArrayD2 rootfr, rootr, eff_porosity;

      // soil fluxes (outputs)
      ArrayD1 eflx_soil_grnd, qflx_evap_grnd, qflx_sub_snow, qflx_dew_snow,
              qflx_dew_grnd, eflx_lwrad_out, eflx_lwrad_net, soil_e_balance;
      // from soil temp - used in soil_e_balance 
      ArrayD2 fact;

      // surface albedo and snicar
      // required for SurfaceAlbedo kernels
      ArrayI1 snl_top, snl_btm, ncan, flg_nosnl;
      ArrayI2 snw_rds_lcl;
      ArrayD1 mu_not, coszen;
      ArrayD2 fabd_sun, fabd_sha, fabi_sun, fabi_sha, albsnd, albsni,
              tsai_z, albout_lcl, flx_slrd_lcl, flx_slri_lcl,
              h2osoi_ice_lcl, h2osoi_liq_lcl;
      ArrayD3 mss_cnc_aer_in_fdb, flx_absd_snw, flx_absi_snw, flx_abs_lcl,
              g_star, omega_star, tau_star;

      // outputs from soil temp/snow hydro
      ArrayI2 imelt;
      ArrayD1 snot_top, dTdz_top, snw_rds_top, sno_liq_top, qflx_sl_top_soil,
              qflx_snow2topsoi, mflx_snowlyr_col, qflx_top_soil, mflx_neg_snow,
              eflx_snomelt, qflx_snomelt, xmf_dummy, xmf_h2osfc_dummy,
              eflx_h2osfc_snow_dummy, qflx_h2osfc_ice_dummy,
              eflx_building_heat_dummy, qflx_snofrz;
      ArrayD2 qflx_snofrz_lyr;

      // main exchange variables
      // transpiration
      // vegetation/soil water exchange (m H2O/s) (+ = to atm)
      ArrayD2 qflx_rootsoi;

      // grid data
      // may not stay in ELM state
      ArrayD2 dz;
      ArrayD2 zsoi;
      ArrayD2 zisoi;

      // need to put away
      ArrayI1 vtype;
      // lat/lon
      double lat{0.0}, lon{0.0};
      ELM::LandType Land;
      
      // view of struct PSNVegData
      ArrayPSN1 psn_pft;
      
      // snow and veg indicators
      ArrayB1 veg_active, do_capsnow;
  };

} // namespace ELM

#include "elm_state_impl.hh"
