// functions derived from CanopyTemperatureMod.F90

#pragma once

namespace ELM::canopy_temperature {

template <typename ArrayD1>
ACCELERATE
void old_ground_temp(const LandType& Land, const double& t_h2osfc, const ArrayD1 t_soisno,
                     double& t_h2osfc_bef, ArrayD1 tssbef)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevgrnd;
  using ELMdims::nlevurb;
  using ELMconst::SPVAL;
  
  if (!Land.lakpoi) {
    for (int i = 0; i < nlevgrnd + nlevsno; i++) {
      if ((Land.ctype == LND::icol_sunwall || Land.ctype == LND::icol_shadewall || Land.ctype == LND::icol_roof) && i > nlevurb) {
        tssbef(i) = SPVAL;
      } else {
        tssbef(i) = t_soisno(i);
      }
      // record t_h2osfc prior to updating
      t_h2osfc_bef = t_h2osfc;
    }
  }
} // old_ground_temp

template <typename ArrayD1>
ACCELERATE
void ground_temp(const LandType& Land, const int& snl, const double& frac_sno_eff,
                 const double& frac_h2osfc, const double& t_h2osfc,
                 const ArrayD1 t_soisno, double& t_grnd)
{
  using ELMdims::nlevsno;

  // ground temperature is weighted average of exposed soil, snow, and h2osfc
  if (!Land.lakpoi) {
    if (snl > 0) {
      t_grnd = frac_sno_eff * t_soisno(nlevsno - snl) + (1.0 - frac_sno_eff - frac_h2osfc) * t_soisno(nlevsno) +
               frac_h2osfc * t_h2osfc;
    } else {
      t_grnd = (1.0 - frac_h2osfc) * t_soisno(nlevsno) + frac_h2osfc * t_h2osfc;
    }
  }
} // ground_temp

template <typename ArrayD1>
ACCELERATE
void calc_soilalpha(const LandType& Land, const double& frac_sno, const double& frac_h2osfc,
                    const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 dz, const ArrayD1 t_soisno,
                    const ArrayD1 watsat, const ArrayD1 sucsat, const ArrayD1 bsw, const ArrayD1 watdry,
                    const ArrayD1 watopt, double& qred, double& hr, double& soilalpha)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevbed;
  using ELMdims::nlevsoi;
  using ELMconst::SPVAL;

  static constexpr double smpmin{-1.e8}; // restriction for min of soil potential (mm)
  
  qred = 1.0; // soil surface relative humidity

  if (!Land.lakpoi) {

    // urban off for now
    //double hr_road_perv; // relative humidity for urban pervious road
    // rootfr_road_perv[nlevgrnd]   [double] fraction of roots in each soil layer for urban pervious road
    // rootr_road_perv[nlevgrnd]    [double] effective fraction of roots in each soil layer for urban pervious road
    double fac;          // soil wetness of surface layer
    double wx;           // partial volume of ice and water of surface layer
    double psit;         // negative potential of soil
    double eff_porosity; //  effective porosity in layer
    double vol_ice;      //  partial volume of ice lens in layer
    double vol_liq;      //  partial volume of liquid water in layer

    //if (Land.ctype == LND::icol_road_perv) {
    //  hr_road_perv = 0.0;
    //}
    if (Land.ltype != LND::istwet && Land.ltype != LND::istice && Land.ltype != LND::istice_mec) {
      if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
        wx = (h2osoi_liq(nlevsno) / ELMconst::DENH2O + h2osoi_ice(nlevsno) / ELMconst::DENICE) / dz(nlevsno);
        fac = std::min(1.0, wx / watsat(0));
        fac = std::max(fac, 0.01);
        psit = -sucsat(0) * pow(fac, (-bsw(0)));
        psit = std::max(smpmin, psit);
        // modify qred to account for h2osfc
        hr = exp(psit / ELMconst::ROVERG / t_soisno(nlevsno));
        qred = (1.0 - frac_sno - frac_h2osfc) * hr + frac_sno + frac_h2osfc;
        soilalpha = qred;
      //} else if (Land.ctype == LND::icol_road_perv) {

      //  // Pervious road depends on water in total soil column
      //  for (int j = 0; j < nlevbed; j++) {
      //    if (t_soisno(j + nlevsno) >= ELMconst::TFRZ) {
      //      vol_ice = std::min(watsat(j), h2osoi_ice(j + nlevsno) / (dz(j + nlevsno) * ELMconst::DENICE));
      //      eff_porosity = watsat(j) - vol_ice;
      //      vol_liq = std::min(eff_porosity, h2osoi_liq(j + nlevsno) / (dz(j + nlevsno) * ELMconst::DENH2O));
      //      fac = std::min(std::max(vol_liq - watdry(j), 0.0) / (watopt(j) - watdry(j)), 1.0);
      //    } else {
      //      fac = 0.0;
      //    }
      //    rootr_road_perv(j) = rootfr_road_perv(j) * fac;
      //    hr_road_perv = hr_road_perv + rootr_road_perv(j);
      //  }
      //  // Allows for sublimation of snow or dew on snow
      //  qred = (1.0 - frac_sno) * hr_road_perv + frac_sno;
      //  // Normalize root resistances to get layer contribution to total ET
      //  if (hr_road_perv > 0.0) {
      //    for (int j = 0; j < nlevsoi; j++) {
      //      rootr_road_perv(j) = rootr_road_perv(j) / hr_road_perv;
      //    }
      //  }
        // soilalpha_u [double] Urban factor that reduces ground saturated specific humidity (-)
        // off until urban code ready
        //soilalpha_u = qred;
      } else if (Land.ctype == LND::icol_sunwall || Land.ctype == LND::icol_shadewall) {
        qred = 0.0;
        //soilalpha_u = SPVAL;
      } else if (Land.ctype == LND::icol_roof || Land.ctype == LND::icol_road_imperv) {
        qred = 1.0;
        //soilalpha_u = SPVAL;
      }
    } else {
      soilalpha = SPVAL;
    }
  }
} // calc_soilalpha

template <typename ArrayD1>
ACCELERATE
void calc_soilbeta(const LandType& Land, const double& frac_sno, const double& frac_h2osfc, const ArrayD1 watsat,
                   const ArrayD1 watfc, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 dz,
                   double& soilbeta)
{
  surface_resistance::calc_soilevap_stress(Land, frac_sno, frac_h2osfc, watsat, watfc, h2osoi_liq, h2osoi_ice, dz,
                                           soilbeta);
} // calc_soilbeta()

template <typename ArrayD1>
ACCELERATE
void humidities(const LandType& Land, const int& snl, const double& forc_q, const double& forc_pbot,
                const double& t_h2osfc, const double& t_grnd, const double& frac_sno, const double& frac_sno_eff,
                const double& frac_h2osfc, const double& qred, const double& hr, const ArrayD1 t_soisno,
                double& qg_snow, double& qg_soil, double& qg, double& qg_h2osfc, double& dqgdT)
{
  using ELMdims::nlevsno;

  if (!Land.lakpoi) {

    double eg;      // water vapor pressure at temperature T [pa]
    double qsatg;   // saturated humidity [kg/kg]
    double degdT;   // d(eg)/dT
    double qsatgdT; // d(qsatg)/dT

    if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
      qsat(t_soisno(nlevsno - snl), forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_snow = qsatg;
      dqgdT = frac_sno * qsatgdT;
      qsat(t_soisno(nlevsno), forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > hr * qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_soil = hr * qsatg;
      dqgdT = dqgdT + (1.0 - frac_sno - frac_h2osfc) * hr * qsatgdT;

      // to be consistent with hs_top values in SoilTemp, set qg_snow to qg_soil for snl = 0 case
      // this ensures hs_top_snow will equal hs_top_soil
      if (snl == 0) {
        qg_snow = qg_soil;
        dqgdT = (1.0 - frac_h2osfc) * hr * dqgdT;
      }
      qsat(t_h2osfc, forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_h2osfc = qsatg;
      dqgdT = dqgdT + frac_h2osfc * qsatgdT;
      qg = frac_sno_eff * qg_snow + (1.0 - frac_sno_eff - frac_h2osfc) * qg_soil + frac_h2osfc * qg_h2osfc;
    } else {
      qsat(t_grnd, forc_pbot, eg, degdT, qsatg, qsatgdT);
      qg = qred * qsatg;
      dqgdT = qred * qsatgdT;

      if (qsatg > forc_q && forc_q > qred * qsatg) {
        qg = forc_q;
        dqgdT = 0.0;
      }
      qg_snow = qg;
      qg_soil = qg;
      qg_h2osfc = qg;
    }
  }
} // humidities

template <typename ArrayD1, typename SubviewD1>
ACCELERATE
void ground_properties(const LandType& Land, const int& snl, const double& frac_sno, const double& forc_th,
                       const double& forc_q, const double& elai, const double& esai, const double& htop,
                       const SubviewD1 displar, const SubviewD1 z0mr, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice,
                       double& emg, double& emv, double& htvp, double& z0mg, double& z0hg, double& z0qg, double& z0mv,
                       double& z0hv, double& z0qv, double& thv, double& z0m, double& displa)
{
  using ELMdims::nlevsno;

  if (!Land.lakpoi) {
    double avmuir; // ir inverse optical depth per unit leaf area

    // Ground emissivity - only calculate for non-urban landunits
    // Urban emissivities are currently read in from data file
    if (!Land.urbpoi) {
      if (Land.ltype == LND::istice || Land.ltype == LND::istice_mec) {
        emg = 0.97;
      } else {
        emg = (1.0 - frac_sno) * 0.96 + frac_sno * 0.97;
      }
    }

    // Vegetation Emissivity
    avmuir = 1.0;
    emv = 1.0 - exp(-(elai + esai) / avmuir);

    // Latent heat. We arbitrarily assume that the sublimation occurs
    // only as h2osoi_liq = 0
    htvp = ELMconst::HVAP;
    if (h2osoi_liq(nlevsno - snl) <= 00 && h2osoi_ice(nlevsno - snl) > 0.0) {
      htvp = ELMconst::HSUB;
    }

    // Ground roughness lengths over non-lake columns (includes bare ground, ground
    // underneath canopy, wetlands, etc.)
    if (frac_sno > 0.0) {
      z0mg = ELMconst::ZSNO;
    } else {
      z0mg = ELMconst::ZLND;
    }
    z0hg = z0mg; // initial set only
    z0qg = z0mg; // initial set only
    z0m = z0mr(Land.vtype) * htop;
    displa = displar(Land.vtype) * htop;

    // vegetation roughness lengths
    z0mv = z0m;
    z0hv = z0mv;
    z0qv = z0mv;

    // Virtual potential temperature
    thv = forc_th * (1.0 + 0.61 * forc_q);
  }
} // ground_properties

ACCELERATE
void forcing_height(const LandType& Land, const bool& veg_active, const int& frac_veg_nosno,
                    const double& forc_hgt_u, const double& forc_hgt_t, const double& forc_hgt_q,
                    const double& z0m, const double& z0mg, const double& forc_t, const double& displa,
                    double& forc_hgt_u_patch, double& forc_hgt_t_patch,
                    double& forc_hgt_q_patch, double& thm)
{
  // Make forcing height a pft-level quantity that is the atmospheric forcing
  // height plus each pft's z0m+displa
  if (veg_active) {
    if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
      if (frac_veg_nosno == 0) {
        forc_hgt_u_patch = forc_hgt_u + z0mg + displa;
        forc_hgt_t_patch = forc_hgt_t + z0mg + displa;
        forc_hgt_q_patch = forc_hgt_q + z0mg + displa;
      } else {
        forc_hgt_u_patch = forc_hgt_u + z0m + displa;
        forc_hgt_t_patch = forc_hgt_t + z0m + displa;
        forc_hgt_q_patch = forc_hgt_q + z0m + displa;
      }
    } else if (Land.ltype == LND::istwet || Land.ltype == LND::istice || Land.ltype == LND::istice_mec) {
      forc_hgt_u_patch = forc_hgt_u + z0mg;
      forc_hgt_t_patch = forc_hgt_t + z0mg;
      forc_hgt_q_patch = forc_hgt_q + z0mg;
    } else if (Land.ltype == LND::istdlak) {
      forc_hgt_u_patch = forc_hgt_u;
      forc_hgt_t_patch = forc_hgt_t;
      forc_hgt_q_patch = forc_hgt_q;
    } else if (Land.urbpoi) {
      //  z_0_town         [double] momentum roughness length of urban landunit (m)
      //  z_d_town         [double] displacement height of urban landunit (m)
      // hardwired until urban area code completed
      const double z_0_town{0.0}, z_d_town{0.0};
      forc_hgt_u_patch = forc_hgt_u + z_0_town + z_d_town;
      forc_hgt_t_patch = forc_hgt_t + z_0_town + z_d_town;
      forc_hgt_q_patch = forc_hgt_q + z_0_town + z_d_town;
    }
  }

  thm = forc_t + 0.0098 * forc_hgt_t_patch;
} // forcing_height

ACCELERATE
void init_energy_fluxes(const LandType& Land, double& eflx_sh_tot, double& eflx_lh_tot,
                       double& eflx_sh_veg, double& qflx_evap_tot, double& qflx_evap_veg,
                       double& qflx_tran_veg)
{
  // Initial set (needed for history tape fields)
  eflx_sh_tot = 0.0;

  // off for now
  // eflx_sh_tot_u  [double] urban total sensible heat flux (W/m**2) [+ to atm]
  // eflx_sh_tot_r  [double] rural total sensible heat flux (W/m**2) [+ to atm]
  // eflx_lh_tot_u  [double] urban total latent heat flux (W/m**2)  [+ to atm]
  // eflx_lh_tot_r  [double] rural total latent heat flux (W/m**2)  [+ to atm]

  //if (Land.urbpoi) {
  //  eflx_sh_tot_u = 0.0;
  //} else if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
  //  eflx_sh_tot_r = 0.0;
  //}
  eflx_lh_tot = 0.0;
  //if (Land.urbpoi) {
  //  eflx_lh_tot_u = 0.0;
  //} else if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
  //  eflx_lh_tot_r = 0.0;
  //}
  eflx_sh_veg = 0.0;
  qflx_evap_tot = 0.0;
  qflx_evap_veg = 0.0;
  qflx_tran_veg = 0.0;
} // init_energy_fluxes

} // namespace ELM::canopy_temperature
