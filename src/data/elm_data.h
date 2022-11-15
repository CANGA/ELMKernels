
#pragma once

#include "land_data.h"

namespace ELM {

template<typename ArrayB1, typename ArrayI1, typename ArrayI2,
         typename ArrayD1, typename ArrayD2, typename ArrayD3>
struct ConstitutiveData {
  
  ConstitutiveData(const int ncells);
  ~ConstitutiveData() = default;

  // forcing data
  ArrayD1 forc_tbot, forc_thbot, forc_pbot, forc_qbot, forc_lwrad,
          forc_u, forc_v, forc_hgt, forc_hgt_u_patch, forc_hgt_t_patch,
          forc_hgt_q_patch;

  // prescribed sat phenology
  ArrayD1 tlai, tsai, elai, esai, htop, hbot;

  // frac veg w and w/o albedo consideration
  ArrayI1 frac_veg_nosno_alb, frac_veg_nosno;

  ArrayD2 swe_old, frac_iceold;
  
  ArrayD1 fwet, fdry, frac_h2osfc, frac_sno_eff;

  // for cansunshade
  ArrayD1 laisun, laisha;
  
  ArrayD2 parsun_z, parsha_z, laisun_z, laisha_z;

  // for surface rad
  ArrayD1 sabg_soil, sabg_snow, sabg, sabv, fsa, fsr;
  
  ArrayD2 sabg_lyr;

  // from albedo and snicar
  // recalculated every timestep
  ArrayD1 coszen, xmf, xmf_h2osfc;
  
  ArrayD2 albsnd, albsni, tlai_z, fsun_z, fabd_sun_z,
          fabd_sha_z, fabi_sun_z, fabi_sha_z, ftdd,
          ftid, ftii, fabd, fabi, albsod, albsoi,
          albgrd, albgri, flx_absdv, flx_absdn,
          flx_absiv, flx_absin, albd, albi;
 
  ArrayI2 imelt;

  ArrayD1 soilbeta, qg_snow, qg_soil, qg, qg_h2osfc,
          dqgdT, htvp, emg, emv, z0mg, z0hg, z0qg,
          z0mv, z0hv, z0qv, thv, z0m, displa, thm;

  ArrayD2 tssbef;

  // bareground fluxes
  ArrayD1 dlrad, ulrad, t_ref2m, q_ref2m, rh_ref2m,
          cgrnds, cgrndl, cgrnd;

  // canopy fluxes
  // need to get/calculate
  ArrayD1 altmax_indx, altmax_lastyear_indx, t10;

  ArrayD1 vcmaxcintsha, vcmaxcintsun, btran, t_veg;

  ArrayD2 rootr, eff_porosity;

  // from soil temp - used in soil_e_balance
  // could be moved into more local scope
  ArrayD3 fact;

  ArrayB1 do_capsnow;
};

template<typename ArrayB1, typename ArrayI1, typename ArrayD1,
         typename ArrayD2, typename ArrayPSN1>
struct ConstantData {

  ConstantData(const int ncells);
  ~ConstantData() = default;

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

  // root fraction
  // function of pft type
  ArrayD2 rootfr;

  // view of struct PSNVegData
  ArrayPSN1 psn_pft;

  // need to put away
  ArrayI1 vtype;
        
  // veg indicator
  ArrayB1 veg_active;
};


template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
struct StateData {
  
  StateData(const int ncell);
  ~StateData() = default;

  // snow
  // snow layers, depth, fraction SCA, integrated snow, snow radius(age proxy)
  ArrayI1 snl; 
  ArrayD1 snow_depth, frac_sno, int_snow;
  ArrayD2 snw_rds;

  // water
  // soil water and soil ice in [kg/m2] and volumetric soil water in [m3/m3]
  ArrayD2 h2osoi_liq, h2osoi_ice, h2osoi_vol;
  // canopy water, SWE, and standing surface water
  ArrayD1 h2ocan, h2osno, h2osfc;

  // temperature
  ArrayD2 t_soisno;
  // temp of ground interface, surface water at current and previous timestep
  // t_grnd should be calculated from incoming t_soisno - leave for now
  ArrayD1 t_grnd, t_h2osfc, t_h2osfc_bef;

  // number of canopy layers above snow for radiative transfer
  ArrayI1 nrad;
  
  // grid data
  // may not stay here
  // subsurface layer data is likely constant in time
  // snow data is variable in time
  ArrayD2 dz, zsoi, zisoi;
};


template<typename ArrayD1, typename ArrayD2>
struct FluxData {

  FluxData(const int ncells);
  ~FluxData() = default;

  // OUTPUTS
  // output at every step
  ArrayD1 eflx_sh_tot, eflx_lh_tot, eflx_sh_veg,
          qflx_evap_tot, qflx_evap_veg, qflx_tran_veg,
          eflx_sh_grnd, eflx_sh_snow, eflx_sh_soil,
          qflx_evap_soi, qflx_ev_snow, qflx_ev_soil,
          qflx_ev_h2osfc, qflx_snwcp_liq, qflx_snwcp_ice,
          qflx_snow_grnd, qflx_rain_grnd, qflx_snow_melt,
          eflx_soil_grnd, qflx_evap_grnd, qflx_sub_snow,
          qflx_dew_snow, qflx_dew_grnd, eflx_lwrad_out,
          eflx_lwrad_net, soil_e_balance, qflx_sl_top_soil,
          qflx_snow2topsoi, mflx_snowlyr_col, qflx_top_soil,
          mflx_neg_snow, eflx_snomelt, qflx_snomelt,
          eflx_h2osfc_snow, qflx_h2osfc_ice, qflx_snofrz;

  ArrayD2 qflx_snofrz_lyr, qflx_rootsoi;
};

} // namespace ELM
