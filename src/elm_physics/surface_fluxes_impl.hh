
#pragma once

namespace ELM::surface_fluxes {

// calculate previous t_grnd for flux correction
// tssbef_snotop = tssbef(ELM::nlevsno-snl)
// tssbef_soitop = tssbef(EML::nlevsno)
ACCELERATE
double prev_tgrnd(const int& snl, const double& frac_sno_eff,
                  const double& frac_h2osfc, const double& t_h2osfc_bef,
                  const double& tssbef_snotop, const double& tssbef_soitop)
{
  if (snl > 0) {
    return frac_sno_eff * tssbef_snotop + (1.0 - frac_sno_eff - frac_h2osfc) * tssbef_soitop +
           frac_h2osfc * t_h2osfc_bef;
  } else {
    return (1.0 - frac_h2osfc) * tssbef_soitop + frac_h2osfc * t_h2osfc_bef;
  }
}

// calculate difference in soil temperature from last time step, for flux corrections
ACCELERATE
double delta_t(const double& t_grnd, const double& t_grnd0) { return t_grnd - t_grnd0; }

// get the ratio of evaporative demand to h2o availability
// if demand is greater than availability, or 1.0 otherwise.
// !!!! ATTENTION: call AFTER qflx_evap_soi is updated by initial_flux_calc()
// h2osoi_ice_snotop = (ELM::nlevsno-snl)
// h2osoi_liq_snotop = (ELM::nlevsno-snl)
ACCELERATE
double evap_ratio(const double& h2osoi_ice_snotop, const double& h2osoi_liq_snotop,
                  const double& dtime, const double& qflx_evap_soi)
{
  // Determine ratio of topsoil_evap_tot
  double egsmax = (h2osoi_ice_snotop + h2osoi_liq_snotop) / dtime;
  // added to trap very small negative soil water,ice
  if (egsmax < 0.0) {
    egsmax = 0.0;
  }
  if (qflx_evap_soi > egsmax) {
    return egsmax / qflx_evap_soi;
  } else {
    return 1.0;
  }
}

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
// tssbef_snotop = tssbef(ELM::nlevsno-snl)
// tssbef_soitop = tssbef(EML::nlevsno)
ACCELERATE
void initial_flux_calc(const bool& urbpoi, const int& snl, const double& frac_sno_eff,
                       const double& frac_h2osfc, const double& t_h2osfc_bef,
                       const double& tssbef_snotop, const double& tssbef_soitop,
                       const double& t_grnd, const double& cgrnds, const double& cgrndl,
                       double& eflx_sh_grnd, double& qflx_evap_soi, double& qflx_ev_snow,
                       double& qflx_ev_soil, double& qflx_ev_h2osfc)
{
  const double t_grnd0{prev_tgrnd(snl, frac_sno_eff, frac_h2osfc, t_h2osfc_bef, tssbef_snotop, tssbef_soitop)};
  const double tinc{delta_t(t_grnd, t_grnd0)};
  eflx_sh_grnd += tinc * cgrnds;
  qflx_evap_soi += tinc * cgrndl;
  // set ev_snow, ev_soil for urban landunits here
  if (!urbpoi) {
    qflx_ev_snow += tinc * cgrndl;
    qflx_ev_soil += tinc * cgrndl;
    qflx_ev_h2osfc += tinc * cgrndl;
  } else {
    qflx_ev_snow = qflx_evap_soi;
    qflx_ev_soil = 0.0;
    qflx_ev_h2osfc = 0.0;
  }
}

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
// h2osoi_ice_snotop = (ELM::nlevsno-snl)
// h2osoi_liq_snotop = (ELM::nlevsno-snl)
// tssbef_snotop = tssbef(ELM::nlevsno-snl)
// tssbef_soitop = tssbef(EML::nlevsno)
ACCELERATE
void update_surface_fluxes(const bool& urbpoi, const bool& do_capsnow, const int& snl, const double& dtime,
                           const double& t_grnd, const double& htvp, const double& frac_sno_eff,
                           const double& frac_h2osfc, const double& t_h2osfc_bef, const double& sabg_soil,
                           const double& sabg_snow, const double& dlrad, const double& frac_veg_nosno,
                           const double& emg, const double& forc_lwrad, const double& tssbef_snotop,
                           const double& tssbef_soitop, const double& h2osoi_ice_snotop,
                           const double& h2osoi_liq_snotop, const double& eflx_sh_veg, const double& qflx_evap_veg,
                           double& qflx_evap_soi, double& eflx_sh_grnd, double& qflx_ev_snow, double& qflx_ev_soil,
                           double& qflx_ev_h2osfc, double& eflx_soil_grnd, double& eflx_sh_tot,
                           double& qflx_evap_tot, double& eflx_lh_tot, double& qflx_evap_grnd,
                           double& qflx_sub_snow, double& qflx_dew_snow, double& qflx_dew_grnd,
                           double& qflx_snwcp_liq, double& qflx_snwcp_ice)
{
  const double egirat{evap_ratio(h2osoi_ice_snotop, h2osoi_liq_snotop, dtime, qflx_evap_soi)};
  // Correct soil fluxes for possible evaporation in excess of top layer water
  // excess energy is added to the sensible heat flux from soil
  if (egirat < 1.0) {
    const double save_qflx_evap_soi{qflx_evap_soi};
    qflx_evap_soi *= egirat;
    eflx_sh_grnd += (save_qflx_evap_soi - qflx_evap_soi) * htvp;
    qflx_ev_snow *= egirat;
    qflx_ev_soil *= egirat;
    qflx_ev_h2osfc *= egirat;
  }

  // Ground heat flux
  if (!urbpoi) {
    const double lw_grnd =
        (frac_sno_eff * pow(tssbef_snotop, 4.0) + (1.0 - frac_sno_eff - frac_h2osfc) * pow(tssbef_soitop, 4.0) +
         frac_h2osfc * pow(t_h2osfc_bef, 40));

    const double t_grnd0{prev_tgrnd(snl, frac_sno_eff, frac_h2osfc, t_h2osfc_bef, tssbef_snotop, tssbef_soitop)};
    const double tinc{delta_t(t_grnd, t_grnd0)};
    eflx_soil_grnd = ((1.0 - frac_sno_eff) * sabg_soil + frac_sno_eff * sabg_snow) + dlrad +
                     (1.0 - frac_veg_nosno) * emg * forc_lwrad - emg * ELMconst::STEBOL * lw_grnd -
                     pow(emg * ELMconst::STEBOL * t_grnd0, 3.0) * (4.0 * tinc) - (eflx_sh_grnd + qflx_evap_soi * htvp);
  } else {
    // For all urban columns we use the net longwave radiation (eflx_lwrad_net) since
    // the term (emg*sb*tssbef(col_pp%snl+1)**4) is not the upward longwave flux because of
    // interactions between urban columns.
    // eflx_lwrad_del = 4.0 * emg * ELMconst::STEBOL * pow(t_grnd0, 3.0) * tinc;
    // ! Include transpiration term because needed for pervious road
    // ! and wasteheat and traffic flux
    // eflx_soil_grnd = sabg + dlrad - eflx_lwrad_net - eflx_lwrad_del
    //      - (eflx_sh_grnd + qflx_evap_soi * htvp + qflx_tran_veg * ELMconst::HVAP)
    //      + eflx_wasteheat_patch + eflx_heat_from_ac_patch + eflx_traffic_patch;
    // eflx_soil_grnd_u = eflx_soil_grnd;
  }

  // Total fluxes (vegetation + ground)
  eflx_sh_tot = eflx_sh_veg + eflx_sh_grnd;
  qflx_evap_tot = qflx_evap_veg + qflx_evap_soi;
  eflx_lh_tot = ELMconst::HVAP * qflx_evap_veg + htvp * qflx_evap_soi;

  // Assign ground evaporation to sublimation from soil ice or to dew
  // on snow or ground
  qflx_evap_grnd = 0.0;
  qflx_sub_snow = 0.0;
  qflx_dew_snow = 0.0;
  qflx_dew_grnd = 0.0;

  if (qflx_ev_snow >= 0.0) {
    // for evaporation partitioning between liquid evap and ice sublimation,
    // use the ratio of liquid to (liquid+ice) in the top layer to determine split
    if ((h2osoi_liq_snotop + h2osoi_ice_snotop) > 0.0) {
      qflx_evap_grnd = std::max(qflx_ev_snow * (h2osoi_liq_snotop / (h2osoi_liq_snotop + h2osoi_ice_snotop)), 0.0);
    } else {
      qflx_evap_grnd = 0.0;
    }
    qflx_sub_snow = qflx_ev_snow - qflx_evap_grnd;
  } else {
    if (t_grnd < ELMconst::TFRZ) {
      qflx_dew_snow = std::abs(qflx_ev_snow);
    } else {
      qflx_dew_grnd = std::abs(qflx_ev_snow);
    }
  }

  // Update the pft-level qflx_snwcp
  if (snl > 0 && do_capsnow) {
    qflx_snwcp_liq = qflx_snwcp_liq + frac_sno_eff * qflx_dew_grnd;
    qflx_snwcp_ice = qflx_snwcp_ice + frac_sno_eff * qflx_dew_snow;
  }
}

// Outgoing long-wave radiation from vegetation + ground
// For conservation we put the increase of ground longwave to outgoing
// For urban patches, ulrad=0 and (1-fracveg_nosno)=1, and eflx_lwrad_out and eflx_lwrad_net
// are calculated in UrbanRadiation. The increase of ground longwave is added directly
// to the outgoing longwave and the net longwave.
// tssbef_snotop = tssbef(ELM::nlevsno-snl)
// tssbef_soitop = tssbef(EML::nlevsno)
ACCELERATE
void lwrad_outgoing(const bool& urbpoi, const int& snl, const int& frac_veg_nosno, const double& forc_lwrad,
                        const double& frac_sno_eff, const double& tssbef_snotop, const double& tssbef_soitop,
                        const double& frac_h2osfc, const double& t_h2osfc_bef, const double& t_grnd,
                        const double& ulrad, const double& emg, double& eflx_lwrad_out, double& eflx_lwrad_net) {

  if (!urbpoi) {
    const double lw_grnd =
        (frac_sno_eff * pow(tssbef_snotop, 4.0) + (1.0 - frac_sno_eff - frac_h2osfc) * pow(tssbef_soitop, 4.0) +
         frac_h2osfc * pow(t_h2osfc_bef, 4.0));
    const double t_grnd0{prev_tgrnd(snl, frac_sno_eff, frac_h2osfc, t_h2osfc_bef, tssbef_snotop, tssbef_soitop)};
    const double tinc{delta_t(t_grnd, t_grnd0)};
    eflx_lwrad_out = ulrad + (1 - frac_veg_nosno) * (1.0 - emg) * forc_lwrad +
                     (1 - frac_veg_nosno) * emg * ELMconst::STEBOL * lw_grnd +
                     4.0 * emg * ELMconst::STEBOL * pow(t_grnd0, 3.0) * tinc;
    eflx_lwrad_net = eflx_lwrad_out - forc_lwrad;
  } else {
    //  eflx_lwrad_out = eflx_lwrad_out + eflx_lwrad_del;
    //  eflx_lwrad_net = eflx_lwrad_net + eflx_lwrad_del;
  }
}

// Soil Energy balance check
template <typename ArrayD1>
ACCELERATE
double soil_energy_balance(const int& ctype, const int& snl, const double& eflx_soil_grnd, const double& xmf,
                           const double& xmf_h2osfc, const double& frac_h2osfc, const double& t_h2osfc,
                           const double& t_h2osfc_bef, const double& dtime, const double& eflx_h2osfc_to_snow,
                           const double& eflx_building_heat, const double& frac_sno_eff, const ArrayD1 t_soisno,
                           const ArrayD1 tssbef, const ArrayD1 fact) {

  double errsoi = eflx_soil_grnd - xmf - xmf_h2osfc - frac_h2osfc * (t_h2osfc - t_h2osfc_bef) * (t_h2osfc / dtime);
  errsoi += eflx_h2osfc_to_snow;
  // For urban sunwall, shadewall, and roof columns, the "soil" energy balance check
  // must include the heat flux from the interior of the building.
  if (ctype == LND::icol_sunwall || ctype == LND::icol_shadewall || ctype == LND::icol_roof) {
    errsoi += eflx_building_heat;
  }
  for (int j = 0; j < ELM::nlevgrnd + ELM::nlevsno; ++j) {
    if ((ctype != LND::icol_sunwall && ctype != LND::icol_shadewall && ctype != LND::icol_roof) || (j < ELM::nlevurb)) {
      // area weight heat absorbed by snow layers
      if (j >= ELM::nlevsno - snl && j < ELM::nlevsno) {
        errsoi -= frac_sno_eff * (t_soisno(j) - tssbef(j)) / fact(j);
      }
      if (j >= ELM::nlevsno) {
        errsoi -= (t_soisno(j) - tssbef(j)) / fact(j);
      }
    }
  }
  return errsoi;
}

} // namespace ELM::surface_fluxes
