/* physics functions from CanopyTemperatureMod.F90



*/
#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

/* SavePriorGroundTemp()
record t_h2osfc and t_soisno prior to updating

*/
void SavePriorGroundTemp(
  const int& ctype,
  const bool& lakpoi,
  const int& nlevsno,
  const int& nlevgrnd,
  const double* t_soisno_,

  double* tssbef_,
  double* t_h2osfc_bef)
{
  // reference variables with pointer offset to allow use of CLM/ELM index without modification
  const double* tssbef = &tssbef_[nlevsno-1];
  double* t_soisno = &t_soisno_[nlevsno-1];

  if (!lakpoi) {
    for (int i = -nlevsno+1; i <= nlevgrnd; i++) {
      if ((ctype == icol_sunwall || ctype == icol_shadewall || ctype == icol_roof) && i > nlevurb) {
        tssbef[i] = spval;
      } else {
        tssbef[i] = t_soisno[i];
      }
      // record t_h2osfc prior to updating
      t_h2osfc_bef = t_h2osfc;   
    }
  }
}


/* CalculateGroundTemp()

*/
void CalculateGroundTemp(
  const bool& lakpoi,
  const int& nlevsno,
  const int& snl,
  const double& frac_sno_eff, 
  const double& frac_h2osfc,
  const double& t_h2osfc,
  const double* t_soisno_,

  double& t_grnd)
{
  // reference variables with pointer offset to allow use of CLM/ELM index without modification
  const double* t_soisno = &t_soisno_[nlevsno-1];

  if (!lakpoi) {
    if (snl < 0) {
      t_grnd = frac_sno_eff * t_soisno[snl+1] + (1.0 - frac_sno_eff - frac_h2osfc) * t_soisno[1] + frac_h2osfc * t_h2osfc;
    } else {
      t_grnd = (1.0 - frac_h2osfc) * t_soisno[1] + frac_h2osfc * t_h2osfc;
    }
  }
}


/*

*/
void CalculateSoilAlpha(
  )
{
  if (!lakpoi) {
    
    if (ctype == icol_road_perv) {
      hr_road_perv = 0.0
    }

    double qred = 1.0;

    if (ltype != istwet && ltype != istice && ltype != istice_mec) {
      if (ltype == istsoil || ltype == istcrop) {
        wx   = (h2osoi_liq[1] / denh2o + h2osoi_ice[1] / denice) / dz[1];
        fac  = std::min(1.0, wx / watsat[1]);
        fac  = std::max(fac, 0.01);
        psit = -sucsat[1] * pow(fac, (-bsw[1]));
        psit = std::max(smpmin, psit);
        // modify qred to account for h2osfc
        hr   = exp(psit / roverg / t_soisno[1]);
        qred = (1.0 - frac_sno - frac_h2osfc) * hr + frac_sno + frac_h2osfc;
        soilalpha = qred;
      } else if (ctype == icol_road_perv) {
        // Pervious road depends on water in total soil column
        nlevbed = nlev2bed;

        for (int j = 1; j <= nlevbed; j++) {
          if (t_soisno[j] >= tfrz) {
            vol_ice = std::min(watsat[j], h2osoi_ice[j] / (dz[j] * denice));
            eff_porosity = watsat[j] - vol_ice;
            vol_liq = std::min(eff_porosity, h2osoi_liq[j] / (dz[j] * denh2o));
            fac = std::min(std::max(vol_liq - watdry[j], 0.0) / (watopt[j] - watdry[j]), 1.0);
          } else {
            fac = 0.0;
          }
          rootr_road_perv[j] = rootfr_road_perv[j] * fac;
          hr_road_perv = hr_road_perv + rootr_road_perv[j];
        }
        // Allows for sublimation of snow or dew on snow
        qred = (1.0 - frac_sno) * hr_road_perv + frac_sno;
        // Normalize root resistances to get layer contribution to total ET
        if (hr_road_perv > 0.0) {
          for (int j = 1; j <= nlevsoi; j++) {
            rootr_road_perv[j] = rootr_road_perv[j] / hr_road_perv;
          }
        }
        soilalpha_u = qred;
      } else if (ctype == icol_sunwall || ctype == icol_shadewall) {
        qred = 0.0;
        soilalpha_u = spval;
      } else if (ctype == icol_roof || ctype == icol_road_imperv) {
        qred = 1.0;
        soilalpha_u = spval;
      }
    } else {
      soilalpha = spval;
    }
  }
}

/*

*/
void CalculateHumidities(
  )
{
  if (!lakpoi) {
    if (ltype == istsoil || ltype == istcrop) {
      QSat(t_soisno[snl+1], forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_snow = qsatg;
      dqgdT = frac_sno * qsatgdT;
      QSat(t_soisno[1] , forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > hr * qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_soil = hr * qsatg;
      dqgdT = dqgdT + (1.0 - frac_sno - frac_h2osfc) * hr * qsatgdT;

      // to be consistent with hs_top values in SoilTemp, set qg_snow to qg_soil for snl = 0 case
      // this ensures hs_top_snow will equal hs_top_soil
      if (snl >= 0) {
        qg_snow; = qg_soil;
        dqgdT = (1.0 - frac_h2osfc) * hr * dqgdT;
      }
      QSat(t_h2osfc, forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_h2osfc = qsatg;
      dqgdT = dqgdT + frac_h2osfc * qsatgdT;
      qg = frac_sno_eff * qg_snow + (1.0 - frac_sno_eff - frac_h2osfc) * qg_soil + frac_h2osfc * qg_h2osfc;
    } else {
      QSat(t_grnd, forc_pbot, eg, degdT, qsatg, qsatgdT);
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
}
