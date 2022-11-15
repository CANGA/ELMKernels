
#pragma once

#include "elm_constants.h"

namespace ELM {

  template
  <typename ArrayB1, typename ArrayI1, typename ArrayI2, typename ArrayD1,
   typename ArrayD2, typename ArrayD3, typename ArrayPSN1>
  ELMState<ArrayB1, ArrayI1, ArrayI2, ArrayD1, ArrayD2, ArrayD3, ArrayPSN1>::ELMState(const int ncells)
      :
        // recalculated every dt
        // forcing data
        forc_tbot("forc_tbot", ncells),
        forc_thbot("forc_thbot", ncells),
        forc_pbot("forc_pbot", ncells),
        forc_qbot("forc_qbot", ncells),
        forc_lwrad("forc_lwrad", ncells),
        forc_u("forc_u", ncells), //in canflux and baregroundflux - could be calculated on the fly
        forc_v("forc_v", ncells), //in canflux and baregroundflux - could be calculated on the fly
        forc_hgt("forc_hgt", ncells), // not used now because forc_hgt hardwired as 30.0 - will be useful if forc_hgt read from file
        forc_hgt_u_patch("forc_hgt_u_patch", ncells),
        forc_hgt_t_patch("forc_hgt_t_patch", ncells),
        forc_hgt_q_patch("forc_hgt_q_patch", ncells),
        //forc_rho("forc_rho", ncells), // in canflux and baregroundflux - could be calculated on the fly

        // recalculated every dt
        // prescribed sat phenology
        tlai("tlai", ncells),
        tsai("tsai", ncells),
        elai("elai", ncells),
        esai("esai", ncells),
        htop("htop", ncells),
        hbot("hbot", ncells),

        // recalculated every dt
        // frac veg w and w/o albedo consideration
        frac_veg_nosno_alb("frac_veg_nosno_alb", ncells),
        frac_veg_nosno("frac_veg_nosno", ncells),



        // these are constant in time
        // based on soil texture
        // make congruent w Van Genuchten 
        // soil hydraulics
        watsat("watsat", ncells, ELMdims::nlevgrnd),
        sucsat("sucsat", ncells, ELMdims::nlevgrnd),
        bsw("bsw", ncells, ELMdims::nlevgrnd),
        watdry("watdry", ncells, ELMdims::nlevgrnd), // only used in canopy_temperature
        watopt("watopt", ncells, ELMdims::nlevgrnd), // only used in canopy_temperature
        watfc("watfc", ncells, ELMdims::nlevgrnd),

        // these are constant in time
        // topo, microtopography
        n_melt("n_melt", ncells),
        micro_sigma("micro_sigma", ncells),
        topo_slope("topo_slope", ncells),
        topo_std("topo_std", ncells),

        // these are constant in time
        // soil color and texture constants
        //int max_soil_color = 20; // largest option - can't know until NC file is read
        // should be placed somewhere else
        isoicol("isoicol", ncells), // only used in surface_albedo
        albsat("albsat", 20, 2), // only used in surface_albedo
        albdry("albdry", 20, 2), // only used in surface_albedo

        // these are constant in time
        // soil thermal constants
        tkmg("tkmg", ncells, ELMdims::nlevgrnd), // created in soil_texture_hydraulic, only used in soil_thermal and soil_temp
        tkdry("tkdry", ncells, ELMdims::nlevgrnd),
        csol("csol", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno),

        // unchecked below this
        // may not be persistent
    
        // snow variables
        // time-variable state
        // need to retain value between timesteps
        snl("snl", ncells),
        snow_depth("snow_depth", ncells),
        frac_sno("frac_sno", ncells),
        int_snow("int_snow", ncells), // calced in can_hydro, used in snow_hydro
        snw_rds("snw_rds", ncells, ELMdims::nlevsno),

        // recalculated in can_hydro every dt
        swe_old("swe_old", ncells, ELMdims::nlevsno), // calced in can_hydro, used in snow_hydro

        // soil/snow/surface h2o state

        // recalculated in init_timestep every dt
        frac_iceold("frac_iceold", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),

        // exchange state
        // exchange variables: soil water and soil ice in [kg/m2] and volumetric soil water in [m3/m3]
        h2osoi_liq("h2osoi_liq", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_ice("h2osoi_ice", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        h2osoi_vol("h2osoi_vol", ncells, ELMdims::nlevgrnd),

        // time-variable state
        // need to retain value between timesteps
        h2ocan("h2ocan", ncells),
        h2osno("h2osno", ncells),

        // recalced, but albedo needs this - save as time-variable state
        // or call before/from albedo
        fwet("fwet", ncells),
        fdry("fdry", ncells),
        // need to retain value between timesteps
        h2osfc("h2osfc", ncells),
        // recalculated every timestep
        frac_h2osfc("frac_h2osfc", ncells),
        frac_sno_eff("frac_sno_eff", ncells),

        // // water fluxes
        // output at every step
        qflx_snwcp_liq("qflx_snwcp_liq", ncells),
        qflx_snwcp_ice("qflx_snwcp_ice", ncells),
        qflx_snow_grnd("qflx_snow_grnd", ncells),
        qflx_rain_grnd("qflx_rain_grnd", ncells),
        qflx_snow_melt("qflx_snow_melt", ncells),
    
        // soil and snow temp
        // time-variable state
        // need to retain value between timesteps or exchange
        t_soisno("t_soisno", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        // should be calced from t_soisno
        t_grnd("t_grnd", ncells),

        // for can_sun_shade

        // time-variable state
        nrad("nrad", ncells),

        // recalculated every timestep
        laisun("laisun", ncells), // in cansunshade
        laisha("laisha", ncells),
        parsun_z("parsun_z", ncells, ELMdims::nlevcan), // in cansunshade
        parsha_z("parsha_z", ncells, ELMdims::nlevcan),
        laisun_z("laisun_z", ncells, ELMdims::nlevcan),
        laisha_z("laisha_z", ncells, ELMdims::nlevcan),

        // from surface rad
        // recalculated every timestep
        sabg_soil("sabg_soil", ncells),
        sabg_snow("sabg_snow", ncells),
        sabg("sabg", ncells),
        sabv("sabv", ncells),
        fsa("fsa", ncells),
        fsr("fsr", ncells),
        sabg_lyr("sabg_lyr", ncells, ELMdims::nlevsno + 1),

        // from albedo and snicar
        // recalculated every timestep
        tlai_z("tlai_z", ncells, ELMdims::nlevcan), // in albedo
        fsun_z("fsun_z", ncells, ELMdims::nlevcan),
        fabd_sun_z("fabd_sun_z", ncells, ELMdims::nlevcan),
        fabd_sha_z("fabd_sha_z", ncells, ELMdims::nlevcan),
        fabi_sun_z("fabi_sun_z", ncells, ELMdims::nlevcan),
        fabi_sha_z("fabi_sha_z", ncells, ELMdims::nlevcan),
        ftdd("ftdd", ncells, ELMdims::numrad),
        ftid("ftid", ncells, ELMdims::numrad),
        ftii("ftii", ncells, ELMdims::numrad),
        fabd("fabd", ncells, ELMdims::numrad),
        fabi("fabi", ncells, ELMdims::numrad),
        albsod("albsod", ncells, ELMdims::numrad),
        albsoi("albsoi", ncells, ELMdims::numrad),
        albgrd("albgrd", ncells, ELMdims::numrad),
        albgri("albgri", ncells, ELMdims::numrad),
        flx_absdv("flx_absdv", ncells, ELMdims::nlevsno + 1),
        flx_absdn("flx_absdn", ncells, ELMdims::nlevsno + 1),
        flx_absiv("flx_absiv", ncells, ELMdims::nlevsno + 1),
        flx_absin("flx_absin", ncells, ELMdims::nlevsno + 1),
        albd("albd", ncells, ELMdims::numrad),
        albi("albi", ncells, ELMdims::numrad),

        // variables for CanopyTemperature
        // time-variable state
        t_h2osfc("t_h2osfc", ncells),
        t_h2osfc_bef("t_h2osfc_bef", ncells),

        // recalculated every dt
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

        // recalculated every dt
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
        // need to get/calculate
        altmax_indx("altmax_indx", ncells), // NEED!! - active layer thickness index into subsurface
        altmax_lastyear_indx("altmax_lastyear_indx", ncells), // NEED!! - active layer thickness index from previous year
        t10("t10", ncells), // NEED!! - running 10-day mean 2m air temperature

        // recalculated every dt
        vcmaxcintsha("vcmaxcintsha", ncells),
        vcmaxcintsun("vcmaxcintsun", ncells),
        btran("btran", ncells),
        t_veg("t_veg", ncells),

        // constant value, initialized once
        rootfr("rootfr", ncells, ELMdims::nlevgrnd),

        // recalculated every dt
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
        // could be moved into more local scope
        // recalculate every dt
        fact("fact", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno), // factors used in computing tridiagonal matrix

        // surface albedo and snicar
        // required for SurfaceAlbedo kernels
        // recalculate every dt
        coszen("coszen", ncells),
        fabd_sun("fabd_sun", ncells, ELMdims::numrad),
        fabd_sha("fabd_sha", ncells, ELMdims::numrad),
        fabi_sun("fabi_sun", ncells, ELMdims::numrad),
        fabi_sha("fabi_sha", ncells, ELMdims::numrad),
        albsnd("albsnd", ncells, ELMdims::numrad),
        albsni("albsni", ncells, ELMdims::numrad),

        // outputs from soil temp/snow hydro
        // many of these could be removed if not desired for output

        // recalculate every dt
        imelt("imelt", ncells, ELMdims::nlevgrnd + ELMdims::nlevsno),
        xmf("xmf", ncells),
        xmf_h2osfc("xmf_h2osfc", ncells),

        // outputs from surface fluxes
        qflx_sl_top_soil("qflx_sl_top_soil", ncells),
        qflx_snow2topsoi("qflx_snow2topsoi", ncells),
        mflx_snowlyr_col("mflx_snowlyr_col", ncells),
        qflx_top_soil("qflx_top_soil", ncells),
        mflx_neg_snow("mflx_neg_snow", ncells),
        eflx_snomelt("eflx_snomelt", ncells),
        qflx_snomelt("qflx_snomelt", ncells),
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
        // subsurface layer data is likely constant in time
        // snow data is variable in time
        dz("dz", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        zsoi("zsoi", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd),
        zisoi("zisoi", ncells, ELMdims::nlevsno + ELMdims::nlevgrnd + 1),

        // hardwired pft type
        // this is the format phenology reader wants 
        vtype("vtype", ncells),

        // constant 
        psn_pft("psn_pft", ncells),

        veg_active("veg_active", ncells), // not sure this is needed - investigate - set as true for now
        do_capsnow("do_capsnow", ncells) // calced in init_timestep

{}
  

} // namespace ELM
