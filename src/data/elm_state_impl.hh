
#pragma once

#include "elm_constants.h"

namespace ELM {

  template
  <typename ArrayB1, typename ArrayI1, typename ArrayI2, typename ArrayD1,
   typename ArrayD2, typename ArrayD3, typename ArrayPSN1>
  ELMState<ArrayB1, ArrayI1, ArrayI2, ArrayD1, ArrayD2, ArrayD3, ArrayPSN1>::ELMState(const int ncells)
      :
        // forcing data
        forc_tbot("forc_tbot", ncells),
        forc_thbot("forc_thbot", ncells),
        forc_pbot("forc_pbot", ncells),
        forc_qbot("forc_qbot", ncells),
        forc_rh("forc_rh", ncells),
        forc_lwrad("forc_lwrad", ncells),
        forc_rain("forc_rain", ncells),
        forc_snow("forc_snow", ncells),
        forc_u("forc_u", ncells),
        forc_v("forc_v", ncells),
        forc_hgt("forc_hgt", ncells),
        forc_hgt_u("forc_hgt_u", ncells),
        forc_hgt_t("forc_hgt_t", ncells),
        forc_hgt_q("forc_hgt_q", ncells),
        forc_hgt_u_patch("forc_hgt_u_patch", ncells),
        forc_hgt_t_patch("forc_hgt_t_patch", ncells),
        forc_hgt_q_patch("forc_hgt_q_patch", ncells),
        forc_vp("forc_vp", ncells),
        forc_rho("forc_rho", ncells),
        forc_po2("forc_po2", ncells),
        forc_pco2("forc_pco2", ncells),
        forc_solai("forc_solai", ncells, 2),
        forc_solad("forc_solad", ncells, 2),

        // prescribed sat phenology
        tlai("tlai", ncells),
        tsai("tsai", ncells),
        elai("elai", ncells),
        esai("esai", ncells),
        htop("htop", ncells),
        hbot("hbot", ncells),

        // frac veg w and w/o albedo consideration
        frac_veg_nosno_alb("frac_veg_nosno_alb", ncells),
        frac_veg_nosno("frac_veg_nosno", ncells),

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
        swe_old("swe_old", ncells, ELMdims::nlevsno),

        // soil/snow/surface h2o state
        frac_iceold("frac_iceold", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_liq("h2osoi_liq", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_ice("h2osoi_ice", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_vol("h2osoi_vol", ncells, ELMdims::nlevgrnd),
        h2ocan("h2ocan", ncells),
        h2osno("h2osno", ncells),
        fwet("fwet", ncells),
        fdry("fdry", ncells),
        h2osfc("h2osfc", ncells),
        frac_h2osfc("frac_h2osfc", ncells),
        frac_sno_eff("frac_sno_eff", ncells),

        // // water fluxes
        qflx_snwcp_liq("qflx_snwcp_liq", ncells),
        qflx_snwcp_ice("qflx_snwcp_ice", ncells),
        qflx_snow_grnd("qflx_snow_grnd", ncells),
        qflx_rain_grnd("qflx_rain_grnd", ncells),
        qflx_snow_melt("qflx_snow_melt", ncells),
    
        // soil and snow temp
        t_soisno("t_soisno", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        t_grnd("t_grnd", ncells),

        // for can_sun_shade
        nrad("nrad", ncells),
        laisun("laisun", ncells),
        laisha("laisha", ncells),
        tlai_z("tlai_z", ncells, ELMdims::nlevcan),
        fsun_z("fsun_z", ncells, ELMdims::nlevcan),
        fabd_sun_z("fabd_sun_z", ncells, ELMdims::nlevcan),
        fabd_sha_z("fabd_sha_z", ncells, ELMdims::nlevcan),
        fabi_sun_z("fabi_sun_z", ncells, ELMdims::nlevcan),
        fabi_sha_z("fabi_sha_z", ncells, ELMdims::nlevcan),
        parsun_z("parsun_z", ncells, ELMdims::nlevcan),
        parsha_z("parsha_z", ncells, ELMdims::nlevcan),
        laisun_z("laisun_z", ncells, ELMdims::nlevcan),
        laisha_z("laisha_z", ncells, ELMdims::nlevcan),

        // for surface rad
        sabg_soil("sabg_soil", ncells),
        sabg_snow("sabg_snow", ncells),
        sabg("sabg", ncells),
        sabv("sabv", ncells),
        fsa("fsa", ncells),
        fsr("fsr", ncells),
        sabg_lyr("sabg_lyr", ncells, ELMdims::nlevsno + 1),
        ftdd("ftdd", ncells, ELMdims::numrad),
        ftid("ftid", ncells, ELMdims::numrad),
        ftii("ftii", ncells, ELMdims::numrad),
        fabd("fabd", ncells, ELMdims::numrad),
        fabi("fabi", ncells, ELMdims::numrad),
        albsod("albsod", ncells, ELMdims::numrad),
        albsoi("albsoi", ncells, ELMdims::numrad),
        albsnd_hst("albsnd_hst", ncells, ELMdims::numrad),
        albsni_hst("albsni_hst", ncells, ELMdims::numrad),
        albgrd("albgrd", ncells, ELMdims::numrad),
        albgri("albgri", ncells, ELMdims::numrad),
        flx_absdv("flx_absdv", ncells, ELMdims::nlevsno + 1),
        flx_absdn("flx_absdn", ncells, ELMdims::nlevsno + 1),
        flx_absiv("flx_absiv", ncells, ELMdims::nlevsno + 1),
        flx_absin("flx_absin", ncells, ELMdims::nlevsno + 1),
        albd("albd", ncells, ELMdims::numrad),
        albi("albi", ncells, ELMdims::numrad),

        // variables for CanopyTemperature
        t_h2osfc("t_h2osfc", ncells),
        t_h2osfc_bef("t_h2osfc_bef", ncells),
        soilalpha("soilalpha", ncells),
        soilbeta("soilbeta", ncells), 
        qg_snow("qg_snow", ncells),
        qg_soil("qg_soil", ncells),
        qg("qg", ncells),
        qg_h2osfc("qg_h2osfc", ncells),
        dqgdT("dqgdT", ncells),
        htvp("htvp", ncells),
        emg("emg", ncells),
        emv("emv", ncells),
        z0mg("z0mg", ncells),
        z0hg("z0hg", ncells),
        z0qg("z0qg", ncells),
        z0mv("z0mv", ncells),
        z0hv("z0hv", ncells),
        z0qv("z0qv", ncells),
        thv("thv", ncells),
        z0m("z0m", ncells),
        displa("displa", ncells),
        thm("thm", ncells),
        eflx_sh_tot("eflx_sh_tot", ncells),
        eflx_lh_tot("eflx_lh_tot", ncells),
        eflx_sh_veg("eflx_sh_veg", ncells),
        qflx_evap_tot("qflx_evap_tot", ncells),
        qflx_evap_veg("qflx_evap_veg", ncells),
        qflx_tran_veg("qflx_tran_veg", ncells),
        tssbef("tssbef", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno),

        // bareground fluxes
        dlrad("dlrad", ncells),
        ulrad("ulrad", ncells),
        eflx_sh_grnd("eflx_sh_grnd", ncells),
        eflx_sh_snow("eflx_sh_snow", ncells),
        eflx_sh_soil("eflx_sh_soil", ncells),
        eflx_sh_h2osfc("eflx_sh_h2osfc", ncells),
        qflx_evap_soi("qflx_evap_soi", ncells),
        qflx_ev_snow("qflx_ev_snow", ncells),
        qflx_ev_soil("qflx_ev_soil", ncells),
        qflx_ev_h2osfc("qflx_ev_h2osfc", ncells),
        t_ref2m("t_ref2m", ncells),
        q_ref2m("q_ref2m", ncells),
        rh_ref2m("rh_ref2m", ncells),
        cgrnds("cgrnds", ncells),
        cgrndl("cgrndl", ncells),
        cgrnd("cgrnd", ncells),

        // canopy fluxes
        altmax_indx("altmax_indx", ncells), // NEED!! - active layer thickness index into subsurface
        altmax_lastyear_indx("altmax_lastyear_indx", ncells), // NEED!! - active layer thickness index from previous year
        t10("t10", ncells), // NEED!! - running 10-day mean 2m air temperature 
        vcmaxcintsha("vcmaxcintsha", ncells),
        vcmaxcintsun("vcmaxcintsun", ncells),
        btran("btran", ncells),
        t_veg("t_veg", ncells),
        rootfr("rootfr", ncells, ELMdims::nlevgrnd),
        rootr("rootr", ncells, ELMdims::nlevgrnd),
        eff_porosity("eff_porosity", ncells, ELMdims::nlevgrnd),

        // soil fluxes (outputs)
        eflx_soil_grnd("eflx_soil_grnd", ncells),
        qflx_evap_grnd("qflx_evap_grnd", ncells),
        qflx_sub_snow("qflx_sub_snow", ncells),
        qflx_dew_snow("qflx_dew_snow", ncells),
        qflx_dew_grnd("qflx_dew_grnd", ncells),
        eflx_lwrad_out("eflx_lwrad_out", ncells),
        eflx_lwrad_net("eflx_lwrad_net", ncells),
        soil_e_balance("soil_e_balance", ncells),

        // from soil temp - used in soil_e_balance 
        fact("fact", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno), // factors used in computing tridiagonal matrix

        // surface albedo and snicar
        // required for SurfaceAlbedo kernels
        coszen("coszen", ncells),
        fabd_sun("fabd_sun", ncells, ELMdims::numrad),
        fabd_sha("fabd_sha", ncells, ELMdims::numrad),
        fabi_sun("fabi_sun", ncells, ELMdims::numrad),
        fabi_sha("fabi_sha", ncells, ELMdims::numrad),
        albsnd("albsnd", ncells, ELMdims::numrad),
        albsni("albsni", ncells, ELMdims::numrad),
        tsai_z("tsai_z", ncells, ELMdims::nlevcan),

        // outputs from soil temp/snow hydro
        // many of these could be removed if not desired for output
        imelt("imelt", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno),
        qflx_sl_top_soil("qflx_sl_top_soil", ncells),
        qflx_snow2topsoi("qflx_snow2topsoi", ncells),
        mflx_snowlyr_col("mflx_snowlyr_col", ncells),
        qflx_top_soil("qflx_top_soil", ncells),
        mflx_neg_snow("mflx_neg_snow", ncells),
        eflx_snomelt("eflx_snomelt", ncells),
        qflx_snomelt("qflx_snomelt", ncells),
        xmf("xmf", ncells),
        xmf_h2osfc("xmf_h2osfc", ncells),
        eflx_h2osfc_snow("eflx_h2osfc_to_snow", ncells),
        qflx_h2osfc_ice("qflx_h2osfc_to_ice", ncells),
        qflx_snofrz("qflx_snofrz", ncells),
        qflx_snofrz_lyr("qflx_snofrz_lyr", ncells, ELMdims::nlevsno),

        // main exchange variables
        // transpiration
        // vegetation/soil water exchange (m H2O/s) (+ = to atm)
        qflx_rootsoi("qflx_rootsoi", ncells, ELMdims::nlevgrnd),

        // grid data
        // may not stay in ELM state
        dz("dz", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        zsoi("zsoi", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        zisoi("zisoi", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd + 1),

        // hardwired pft type
        // this is the format phenology reader wants 
        vtype("vtype", ncells),


        psn_pft("psn_pft", ncells),

        veg_active("veg_active", ncells),
        do_capsnow("do_capsnow", ncells)

{}


} // namespace ELM
