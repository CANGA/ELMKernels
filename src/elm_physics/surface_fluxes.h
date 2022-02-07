
#pragma once

#include "elm_constants.h"
#include <cmath>

namespace ELM::surface_fluxes {
// calculate previous t_grnd for flux correction
//tssbef_snotop = tssbef(ELM::nlevsno-snl)
//tssbef_soitop = tssbef(EML::nlevsno)
double prev_tgrnd(const int& snl, const double& frac_sno_eff, const double& frac_h2osfc, const double& t_h2osfc_bef, const double& tssbef_snotop, const double& tssbef_soitop);

// calculate difference in soil temperature from last time step, for flux corrections
double delta_t(const double& t_grnd, const double& t_grnd0);

// get the ratio of evaporative demand to h2o availability
// if demand is greater than availability, or 1.0 otherwise.
// !!!! ATTENTION: call AFTER qflx_evap_soi is updated by initial_flux_calc()
//h2osoi_ice_snotop = (ELM::nlevsno-snl)
//h2osoi_liq_snotop = (ELM::nlevsno-snl)
double evap_ratio(const double& h2osoi_ice_snotop, const double& h2osoi_liq_snotop, const double& dtime, const double& qflx_evap_soi);

/*
IN:
urbpoi
snl
frac_sno_eff
frac_h2osfc
t_h2osfc_bef
tssbef_snotop
tssbef_soitop
t_grnd
cgrnds
cgrndl


OUT:
eflx_sh_grnd
qflx_evap_soi
qflx_ev_snow
qflx_ev_soil
qflx_ev_h2osfc
*/
// corrct fluxes for changing temperature during timestep
// must CALL BEFORE any other surface flux functions
//tssbef_snotop = tssbef(ELM::nlevsno-snl)
//tssbef_soitop = tssbef(EML::nlevsno)
void initial_flux_calc(const bool& urbpoi, const int& snl, const double& frac_sno_eff, const double& frac_h2osfc, const double& t_h2osfc_bef, const double& tssbef_snotop, const double& tssbef_soitop,
  const double& t_grnd, const double& cgrnds, const double& cgrndl, double& eflx_sh_grnd, double& qflx_evap_soi, double& qflx_ev_snow, double& qflx_ev_soil,
  double& qflx_ev_h2osfc);


/*
IN:
urbpoi
do_capsnow
snl
dtime
t_grnd
htvp
frac_sno_eff
frac_h2osfc
t_h2osfc_bef
sabg_soil
sabg_snow
dlrad
frac_veg_nosno
emg
forc_lwrad
tssbef_snotop
tssbef_soitop
h2osoi_ice_snotop
h2osoi_liq_snotop
eflx_sh_veg
qflx_evap_veg


OUT:
qflx_evap_soi
eflx_sh_grnd
qflx_ev_snow
qflx_ev_soil
qflx_ev_h2osfc
lw_grnd
eflx_soil_grnd
eflx_sh_tot
qflx_evap_tot
eflx_lh_tot
qflx_evap_grnd
qflx_sub_snow
qflx_dew_snow
qflx_dew_grnd
qflx_snwcp_liq
qflx_snwcp_ice

*/
// calculate new surface fluxes for current timestep
//h2osoi_ice_snotop = (ELM::nlevsno-snl)
//h2osoi_liq_snotop = (ELM::nlevsno-snl)
//tssbef_snotop = tssbef(ELM::nlevsno-snl)
//tssbef_soitop = tssbef(EML::nlevsno)
void update_surface_fluxes(const bool& urbpoi, const bool& do_capsnow, const int& snl, const double& dtime, const double& t_grnd, const double& htvp, const double& frac_sno_eff, const double& frac_h2osfc, const double& t_h2osfc_bef, const double& sabg_soil,
  const double& sabg_snow, const double& dlrad, const double& frac_veg_nosno, const double& emg, const double& forc_lwrad, const double& tssbef_snotop, const double& tssbef_soitop,
  const double& h2osoi_ice_snotop, const double& h2osoi_liq_snotop, 
  const double& eflx_sh_veg,
  const double& qflx_evap_veg, double& qflx_evap_soi, double& eflx_sh_grnd, double& qflx_ev_snow, double& qflx_ev_soil, double& qflx_ev_h2osfc,
  double& eflx_soil_grnd, double& eflx_sh_tot, double& qflx_evap_tot, double& eflx_lh_tot, double& qflx_evap_grnd, double& qflx_sub_snow, double& qflx_dew_snow, double& qflx_dew_grnd, double& qflx_snwcp_liq,
  double& qflx_snwcp_ice);

// Soil Energy balance check
template<typename ArrayD1>
double soil_energy_balance(const bool& ctype, const int& snl, const double& eflx_soil_grnd, const double& xmf, const double& xmf_h2osfc, const double& frac_h2osfc, const double& t_h2osfc,
  const double& t_h2osfc_bef, const double& dtime, const double& eflx_h2osfc_to_snow, const double& eflx_building_heat, const double&frac_sno_eff, const ArrayD1 t_soisno, const ArrayD1 tssbef, const ArrayD1 fact);

// Outgoing long-wave radiation from vegetation + ground
// For conservation we put the increase of ground longwave to outgoing
// For urban patches, ulrad=0 and (1-fracveg_nosno)=1, and eflx_lwrad_out and eflx_lwrad_net 
// are calculated in UrbanRadiation. The increase of ground longwave is added directly 
// to the outgoing longwave and the net longwave.
//tssbef_snotop = tssbef(ELM::nlevsno-snl)
//tssbef_soitop = tssbef(EML::nlevsno)
void lwrad_outgoing(const bool& urbpoi, const int& snl, const int& frac_veg_nosno, const double& forc_lwrad, const double& frac_sno_eff, const double& tssbef_snotop, const double& tssbef_soitop,
  const double& frac_h2osfc, const double& t_h2osfc_bef, const double& t_grnd, const double& ulrad, const double& emg, double& eflx_lwrad_out, double& eflx_lwrad_net);

} // namespace ELM::surface_fluxes

#include "surface_fluxes_impl.hh"
