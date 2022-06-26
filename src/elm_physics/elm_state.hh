
#pragma once

#include "land_data.h"

namespace ELM {

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayPSN1>
  struct ELMState {
  
      ELMState(const int ncells);
      ~ELMState() = default;
      
      // forcing data
      ArrayD1 forc_tbot, forc_thbot, forc_pbot, forc_qbot, forc_rh, forc_lwrad,
              forc_rain, forc_snow, forc_u, forc_v, forc_hgt, forc_hgt_u,
              forc_hgt_t, forc_hgt_q, forc_vp, forc_rho, forc_po2, forc_pco2;
      ArrayD2 forc_solai, forc_solad;
  
      // prescribed sat phenology
      ArrayI1 frac_veg_nosno_alb;
      ArrayD1 tlai, tsai, elai, esai, htop, hbot;
  
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
  
  
  
      // unchecked below this
      // may not be persistent
      
      // snow variables
      ArrayI1 snl;
      ArrayD1 snow_depth, frac_sno, int_snow;
      ArrayD2 snw_rds;
  
      // uncategorized
      ArrayD1 t_grnd, h2ocan, h2osno;
      ArrayI1 frac_veg_nosno;

      // soil h20 state
      ArrayD2 frac_iceold, h2osoi_liq, h2osoi_ice, h2osoi_vol;
  
      // for Canopy hydrology
      ArrayD1 qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd, qflx_rain_grnd,
              fwet, fdry, qflx_snow_melt, h2osfc, frac_h2osfc, frac_sno_eff;
      ArrayD2 swe_old;

      // soil temp
      ArrayD2 t_soisno;

      // grid data
      // may not stay in ELM state
      ArrayD2 dz;
      ArrayD2 zsoi;
      ArrayD2 zisoi;

      ArrayI1 vtype;


      // lat/lon
      double lat{0.0}, lon{0.0};

      ELM::LandType Land;



      // canopy fluxes
      ArrayD2 rootfr;


      // view of struct PSNVegData
      ArrayPSN1 psn_pft;

  
  };

} // namespace ELM

#include "elm_state_impl.hh"
