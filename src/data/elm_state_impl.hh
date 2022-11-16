
#pragma once

#include "elm_constants.h"

namespace ELM {

template<typename ArrayB1, typename ArrayI1, typename ArrayI2, typename ArrayD1,
         typename ArrayD2, typename ArrayD3, typename ArrayPSN1>
ELMState<ArrayB1, ArrayI1, ArrayI2, ArrayD1, ArrayD2, ArrayD3, ArrayPSN1>::
ELMState(size_t ncols,
        const Utils::DomainDecomposition<2>& domain,
        const std::string& filename,
        const ELM::Utils::Date &file_start_time,
        int atm_nsteps) :

    // recalculated every dt
    // forcing data
    forc_tbot("forc_tbot", ncols),
    forc_thbot("forc_thbot", ncols),
    forc_pbot("forc_pbot", ncols),
    forc_qbot("forc_qbot", ncols),
    forc_lwrad("forc_lwrad", ncols),
    forc_u("forc_u", ncols), //in canflux and baregroundflux - could be calculated on the fly
    forc_v("forc_v", ncols), //in canflux and baregroundflux - could be calculated on the fly
    forc_hgt("forc_hgt", ncols), // not used now because forc_hgt hardwired as 30.0 - will be useful if forc_hgt read from file
    forc_hgt_u_patch("forc_hgt_u_patch", ncols),
    forc_hgt_t_patch("forc_hgt_t_patch", ncols),
    forc_hgt_q_patch("forc_hgt_q_patch", ncols),
    //forc_rho("forc_rho", ncols), // in canflux and baregroundflux - could be calculated on the fly

    // recalculated every dt
    // prescribed sat phenology
    tlai("tlai", ncols),
    tsai("tsai", ncols),
    elai("elai", ncols),
    esai("esai", ncols),
    htop("htop", ncols),
    hbot("hbot", ncols),

    // recalculated every dt
    // frac veg w and w/o albedo consideration
    frac_veg_nosno_alb("frac_veg_nosno_alb", ncols),
    frac_veg_nosno("frac_veg_nosno", ncols),

    // these are constant in time
    // based on soil texture
    // make congruent w Van Genuchten 
    // soil hydraulics
    watsat("watsat", ncols, ELMdims::nlevgrnd),
    sucsat("sucsat", ncols, ELMdims::nlevgrnd),
    bsw("bsw", ncols, ELMdims::nlevgrnd),
    watdry("watdry", ncols, ELMdims::nlevgrnd), // only used in canopy_temperature
    watopt("watopt", ncols, ELMdims::nlevgrnd), // only used in canopy_temperature
    watfc("watfc", ncols, ELMdims::nlevgrnd),

    // these are constant in time
    // topo, microtopography
    n_melt("n_melt", ncols),
    micro_sigma("micro_sigma", ncols),
    topo_slope("topo_slope", ncols),
    topo_std("topo_std", ncols),

    // these are constant in time
    // soil color and texture constants
    //int max_soil_color = 20; // largest option - can't know until NC file is read
    // should be placed somewhere else
    isoicol("isoicol", ncols), // only used in surface_albedo
    albsat("albsat", 20, 2), // only used in surface_albedo
    albdry("albdry", 20, 2), // only used in surface_albedo

    // these are constant in time
    // soil thermal constants
    tkmg("tkmg", ncols, ELMdims::nlevgrnd), // created in soil_texture_hydraulic, only used in soil_thermal and soil_temp
    tkdry("tkdry", ncols, ELMdims::nlevgrnd),
    csol("csol", ncols, ELMdims::nlevgrnd + ELMdims::nlevsno),

    // unchecked below this
    // may not be persistent
  
    // snow variables
    // time-variable state
    // need to retain value between timesteps
    snl("snl", ncols),
    snow_depth("snow_depth", ncols),
    frac_sno("frac_sno", ncols),
    int_snow("int_snow", ncols), // calced in can_hydro, used in snow_hydro
    snw_rds("snw_rds", ncols, ELMdims::nlevsno),

    // recalculated in can_hydro every dt
    swe_old("swe_old", ncols, ELMdims::nlevsno), // calced in can_hydro, used in snow_hydro

    // soil/snow/surface h2o state

    // recalculated in init_timestep every dt
    frac_iceold("frac_iceold", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd),

    // exchange state
    // exchange variables: soil water and soil ice in [kg/m2] and volumetric soil water in [m3/m3]
    h2osoi_liq("h2osoi_liq", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd),
    h2osoi_ice("h2osoi_ice", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd),
    h2osoi_vol("h2osoi_vol", ncols, ELMdims::nlevgrnd),

    // time-variable state
    // need to retain value between timesteps
    h2ocan("h2ocan", ncols),
    h2osno("h2osno", ncols),

    // recalced, but albedo needs this - save as time-variable state
    // or call before/from albedo
    fwet("fwet", ncols),
    fdry("fdry", ncols),
    // need to retain value between timesteps
    h2osfc("h2osfc", ncols),
    // recalculated every timestep
    frac_h2osfc("frac_h2osfc", ncols),
    frac_sno_eff("frac_sno_eff", ncols),

    // // water fluxes
    // output at every step
    qflx_snwcp_liq("qflx_snwcp_liq", ncols),
    qflx_snwcp_ice("qflx_snwcp_ice", ncols),
    qflx_snow_grnd("qflx_snow_grnd", ncols),
    qflx_rain_grnd("qflx_rain_grnd", ncols),
    qflx_snow_melt("qflx_snow_melt", ncols),
  
    // soil and snow temp
    // time-variable state
    // need to retain value between timesteps or exchange
    t_soisno("t_soisno", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd),
    // should be calced from t_soisno
    t_grnd("t_grnd", ncols),

    // for can_sun_shade

    // time-variable state
    nrad("nrad", ncols),

    // recalculated every timestep
    laisun("laisun", ncols), // in cansunshade
    laisha("laisha", ncols),
    parsun_z("parsun_z", ncols, ELMdims::nlevcan), // in cansunshade
    parsha_z("parsha_z", ncols, ELMdims::nlevcan),
    laisun_z("laisun_z", ncols, ELMdims::nlevcan),
    laisha_z("laisha_z", ncols, ELMdims::nlevcan),

    // from surface rad
    // recalculated every timestep
    sabg_soil("sabg_soil", ncols),
    sabg_snow("sabg_snow", ncols),
    sabg("sabg", ncols),
    sabv("sabv", ncols),
    fsa("fsa", ncols),
    fsr("fsr", ncols),
    sabg_lyr("sabg_lyr", ncols, ELMdims::nlevsno + 1),

    // from albedo and snicar
    // recalculated every timestep
    tlai_z("tlai_z", ncols, ELMdims::nlevcan), // in albedo
    fsun_z("fsun_z", ncols, ELMdims::nlevcan),
    fabd_sun_z("fabd_sun_z", ncols, ELMdims::nlevcan),
    fabd_sha_z("fabd_sha_z", ncols, ELMdims::nlevcan),
    fabi_sun_z("fabi_sun_z", ncols, ELMdims::nlevcan),
    fabi_sha_z("fabi_sha_z", ncols, ELMdims::nlevcan),
    ftdd("ftdd", ncols, ELMdims::numrad),
    ftid("ftid", ncols, ELMdims::numrad),
    ftii("ftii", ncols, ELMdims::numrad),
    fabd("fabd", ncols, ELMdims::numrad),
    fabi("fabi", ncols, ELMdims::numrad),
    albsod("albsod", ncols, ELMdims::numrad),
    albsoi("albsoi", ncols, ELMdims::numrad),
    albgrd("albgrd", ncols, ELMdims::numrad),
    albgri("albgri", ncols, ELMdims::numrad),
    flx_absdv("flx_absdv", ncols, ELMdims::nlevsno + 1),
    flx_absdn("flx_absdn", ncols, ELMdims::nlevsno + 1),
    flx_absiv("flx_absiv", ncols, ELMdims::nlevsno + 1),
    flx_absin("flx_absin", ncols, ELMdims::nlevsno + 1),
    albd("albd", ncols, ELMdims::numrad),
    albi("albi", ncols, ELMdims::numrad),

    // variables for CanopyTemperature
    // time-variable state
    t_h2osfc("t_h2osfc", ncols),
    t_h2osfc_bef("t_h2osfc_bef", ncols),

    // recalculated every dt
    soilbeta("soilbeta", ncols),
    qg_snow("qg_snow", ncols),
    qg_soil("qg_soil", ncols),
    qg("qg", ncols),
    qg_h2osfc("qg_h2osfc", ncols),
    dqgdT("dqgdT", ncols),
    htvp("htvp", ncols),
    emg("emg", ncols),
    emv("emv", ncols),
    z0mg("z0mg", ncols),
    z0hg("z0hg", ncols),
    z0qg("z0qg", ncols),
    z0mv("z0mv", ncols),
    z0hv("z0hv", ncols),
    z0qv("z0qv", ncols),
    thv("thv", ncols),
    z0m("z0m", ncols),
    displa("displa", ncols),
    thm("thm", ncols),
    eflx_sh_tot("eflx_sh_tot", ncols),
    eflx_lh_tot("eflx_lh_tot", ncols),
    eflx_sh_veg("eflx_sh_veg", ncols),
    qflx_evap_tot("qflx_evap_tot", ncols),
    qflx_evap_veg("qflx_evap_veg", ncols),
    qflx_tran_veg("qflx_tran_veg", ncols),
    tssbef("tssbef", ncols, ELMdims::nlevgrnd + ELMdims::nlevsno),

    // recalculated every dt
    // bareground fluxes
    dlrad("dlrad", ncols),
    ulrad("ulrad", ncols),
    eflx_sh_grnd("eflx_sh_grnd", ncols),
    eflx_sh_snow("eflx_sh_snow", ncols),
    eflx_sh_soil("eflx_sh_soil", ncols),
    eflx_sh_h2osfc("eflx_sh_h2osfc", ncols),
    qflx_evap_soi("qflx_evap_soi", ncols),
    qflx_ev_snow("qflx_ev_snow", ncols),
    qflx_ev_soil("qflx_ev_soil", ncols),
    qflx_ev_h2osfc("qflx_ev_h2osfc", ncols),
    t_ref2m("t_ref2m", ncols),
    q_ref2m("q_ref2m", ncols),
    rh_ref2m("rh_ref2m", ncols),
    cgrnds("cgrnds", ncols),
    cgrndl("cgrndl", ncols),
    cgrnd("cgrnd", ncols),

    // canopy fluxes
    // need to get/calculate
    altmax_indx("altmax_indx", ncols), // NEED!! - active layer thickness index into subsurface
    altmax_lastyear_indx("altmax_lastyear_indx", ncols), // NEED!! - active layer thickness index from previous year
    t10("t10", ncols), // NEED!! - running 10-day mean 2m air temperature

    // recalculated every dt
    vcmaxcintsha("vcmaxcintsha", ncols),
    vcmaxcintsun("vcmaxcintsun", ncols),
    btran("btran", ncols),
    t_veg("t_veg", ncols),

    // constant value, initialized once
    rootfr("rootfr", ncols, ELMdims::nlevgrnd),

    // recalculated every dt
    rootr("rootr", ncols, ELMdims::nlevgrnd),
    eff_porosity("eff_porosity", ncols, ELMdims::nlevgrnd),

    // soil fluxes (outputs)
    eflx_soil_grnd("eflx_soil_grnd", ncols),
    qflx_evap_grnd("qflx_evap_grnd", ncols),
    qflx_sub_snow("qflx_sub_snow", ncols),
    qflx_dew_snow("qflx_dew_snow", ncols),
    qflx_dew_grnd("qflx_dew_grnd", ncols),
    eflx_lwrad_out("eflx_lwrad_out", ncols),
    eflx_lwrad_net("eflx_lwrad_net", ncols),
    soil_e_balance("soil_e_balance", ncols),

    // from soil temp - used in soil_e_balance
    // could be moved into more local scope
    // recalculate every dt
    fact("fact", ncols, ELMdims::nlevgrnd + ELMdims::nlevsno), // factors used in computing tridiagonal matrix

    // surface albedo and snicar
    // required for SurfaceAlbedo kernels
    // recalculate every dt
    coszen("coszen", ncols),
    fabd_sun("fabd_sun", ncols, ELMdims::numrad),
    fabd_sha("fabd_sha", ncols, ELMdims::numrad),
    fabi_sun("fabi_sun", ncols, ELMdims::numrad),
    fabi_sha("fabi_sha", ncols, ELMdims::numrad),
    albsnd("albsnd", ncols, ELMdims::numrad),
    albsni("albsni", ncols, ELMdims::numrad),

    // outputs from soil temp/snow hydro
    // many of these could be removed if not desired for output

    // recalculate every dt
    imelt("imelt", ncols, ELMdims::nlevgrnd + ELMdims::nlevsno),
    xmf("xmf", ncols),
    xmf_h2osfc("xmf_h2osfc", ncols),

    // outputs from surface fluxes
    qflx_sl_top_soil("qflx_sl_top_soil", ncols),
    qflx_snow2topsoi("qflx_snow2topsoi", ncols),
    mflx_snowlyr_col("mflx_snowlyr_col", ncols),
    qflx_top_soil("qflx_top_soil", ncols),
    mflx_neg_snow("mflx_neg_snow", ncols),
    eflx_snomelt("eflx_snomelt", ncols),
    qflx_snomelt("qflx_snomelt", ncols),
    eflx_h2osfc_snow("eflx_h2osfc_to_snow", ncols),
    qflx_h2osfc_ice("qflx_h2osfc_to_ice", ncols),
    qflx_snofrz("qflx_snofrz", ncols),
    qflx_snofrz_lyr("qflx_snofrz_lyr", ncols, ELMdims::nlevsno),

    // main exchange variables
    // transpiration
    // vegetation/soil water exchange (m H2O/s) (+ = to atm)
    qflx_rootsoi("qflx_rootsoi", ncols, ELMdims::nlevgrnd),

    // grid data
    // may not stay in ELM state
    // subsurface layer data is likely constant in time
    // snow data is variable in time
    dz("dz", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd),
    zsoi("zsoi", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd),
    zisoi("zisoi", ncols, ELMdims::nlevsno + ELMdims::nlevgrnd + 1),

    // hardwired pft type
    // this is the format phenology reader wants 
    vtype("vtype", ncols),

    // constant 
    psn_pft("psn_pft", ncols),

    veg_active("veg_active", ncols), // not sure this is needed - investigate - set as true for now
    do_capsnow("do_capsnow", ncols), // calced in init_timestep

    // const in time
    // containers for SNICAR model parameters
    snicar_data(std::make_shared<ELM::SnicarData<ArrayD1, ArrayD2, ArrayD3>>()),
    snw_rds_table(std::make_shared<ELM::SnwRdsTable<ArrayD3>>()),
    // container for PFT type parameters
    pft_data(std::make_shared<ELM::PFTData<ArrayD1>>()),
    // container for aerosol deposition forcing - 1 year repeating
    aerosol_data(std::make_shared<ELM::AerosolDataManager<ArrayD1>>()),

    // time variable
    // containers for aerosol deposition and concentration within snowpack layers
    aerosol_masses(std::make_shared<ELM::AerosolMasses<ArrayD2>>(ncols)),
    aerosol_concentrations(std::make_shared<ELM::AerosolConcentrations<ArrayD2>>(ncols)),
    // container for satellite phenology forcing data
    phen_data(std::make_shared<ELM::PhenologyDataManager<ArrayD2>>(domain, ncols, ELMdims::numveg)),
    // container for atmospheric forcing data
    atm_forcing(std::make_shared<ELM::AtmForcObjects<ArrayD1, ArrayD2>>(filename, file_start_time, atm_nsteps, ncols))
{}
  

} // namespace ELM
