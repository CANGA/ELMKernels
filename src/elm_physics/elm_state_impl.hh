
#pragma once

#include "elm_constants.h"

namespace ELM {

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayPSN1>
  ELMState<ArrayI1, ArrayD1, ArrayD2, ArrayPSN1>::ELMState(const int ncells)
      :
        // forcing data
        forc_tbot("forc_tbot", ncells),
        forc_thbot("forc_thbot", ncells),
        forc_pbot("forc_pbot", ncells),
        forc_qbot("forc_qbot", ncells),
        forc_rh("forc_rh", ncells),
        forc_lwrad("forc_lwrad", ncells),
        forc_solai("forc_solai", ncells, 2),
        forc_solad("forc_solad", ncells, 2),
        forc_rain("forc_rain", ncells),
        forc_snow("forc_snow", ncells),
        forc_u("forc_u", ncells),
        forc_v("forc_v", ncells),
        forc_hgt("forc_hgt", ncells),
        forc_hgt_u("forc_hgt_u", ncells),
        forc_hgt_t("forc_hgt_t", ncells),
        forc_hgt_q("forc_hgt_q", ncells),
        forc_vp("forc_vp", ncells),
        forc_rho("forc_rho", ncells),
        forc_po2("forc_po2", ncells),
        forc_pco2("forc_pco2", ncells),

        // prescribed sat phenology
        frac_veg_nosno_alb("frac_veg_nosno_alb", ncells),
        tlai("tlai", ncells),
        tsai("tsai", ncells),
        elai("elai", ncells),
        esai("esai", ncells),
        htop("htop", ncells),
        hbot("hbot", ncells),

        // soil hydraulics
        watsat("watsat", ncells, ELMdims::nlevgrnd),
        sucsat("sucsat", ncells, ELMdims::nlevgrnd),
        bsw("bsw", ncells, ELMdims::nlevgrnd),
        watdry("watdry", ncells, ELMdims::nlevgrnd),
        watopt("watopt", ncells, ELMdims::nlevgrnd),
        watfc("watfc", ncells, ELMdims::nlevgrnd),
    
        // topo, microtopography
        n_melt("n_melt", ncells),
        micro_sigma("micro_sigma", ncells),
        topo_slope("topo_slope", ncells),
        topo_std("topo_std", ncells),

        // soil color and texture constants
        //int max_soil_color = 20; // largest option - can't know until NC file is read
        isoicol("isoicol", ncells),
        albsat("albsat", 20, 2),
        albdry("albdry", 20, 2),

        // soil thermal constants
        tkmg("tkmg", ncells, ELMdims::nlevgrnd),
        tkdry("tkdry", ncells, ELMdims::nlevgrnd),
        csol("csol", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno),



        // unchecked below this
        // may not be persistent
    
        // snow variables
        snl("snl", ncells),
        snow_depth("snow_depth", ncells),
        frac_sno("frac_sno", ncells),
        int_snow("int_snow", ncells),
        snw_rds("snw_rds", ncells, ELMdims::nlevsno),

        // uncategorized
        t_grnd("t_grnd", ncells),
        h2ocan("h2ocan", ncells),
        h2osno("h2osno", ncells),
        frac_veg_nosno("frac_veg_nosno", ncells),
        frac_iceold("frac_iceold", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_liq("h2osoi_liq", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_ice("h2osoi_ice", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_vol("h2osoi_vol", ncells, ELMdims::nlevgrnd),

        // for Canopy hydrology
        qflx_snwcp_liq("qflx_snwcp_liq", ncells),
        qflx_snwcp_ice("qflx_snwcp_ice", ncells),
        qflx_snow_grnd("qflx_snow_grnd", ncells),
        qflx_rain_grnd("qflx_rain_grnd", ncells),
        fwet("fwet", ncells),
        fdry("fdry", ncells),
        qflx_snow_melt("qflx_snow_melt", ncells),
        h2osfc("h2osfc", ncells),
        frac_h2osfc("frac_h2osfc", ncells),
        frac_sno_eff("frac_sno_eff", ncells),
        swe_old("swe_old", ncells, ELMdims::nlevsno),
    
        // soil and snow temp
        t_soisno("t_soisno", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),

        dz("dz", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        zsoi("zsoi", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        zisoi("zisoi", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd + 1),

        vtype("vtype", ncells),

        rootfr("rootfr", ncells, ELMdims::nlevgrnd),


        psn_pft("psn_pft", ncells)

{}


} // namespace ELM
