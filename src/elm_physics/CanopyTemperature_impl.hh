// functions derived from CanopyTemperatureMod.F90

#pragma once

#include "QSat.h"
#include "SurfaceResistance.h"
#include "clm_constants.h"
#include "landtype.h"
#include <algorithm>
#include <cmath>

namespace ELM {

template <class ArrayD1>
void SaveGroundTemp(const LandType &Land, const double &t_h2osfc, const ArrayD1 t_soisno, double &t_h2osfc_bef,
                    ArrayD1 tssbef) {

  if (!Land.lakpoi) {
    for (int i = 0; i < nlevgrnd + nlevsno; i++) {
      if ((Land.ctype == icol_sunwall || Land.ctype == icol_shadewall || Land.ctype == icol_roof) && i > nlevurb) {
        tssbef[i] = spval;
      } else {
        tssbef[i] = t_soisno[i];
      }
      // record t_h2osfc prior to updating
      t_h2osfc_bef = t_h2osfc;
    }
  }
} // SaveGroundTemp

template <class ArrayD1>
void CalculateGroundTemp(const LandType &Land, const int &snl, const double &frac_sno_eff, const double &frac_h2osfc,
                         const double &t_h2osfc, const ArrayD1 t_soisno, double &t_grnd) {
  // ground temperature is weighted average of exposed soil, snow, and h2osfc
  if (!Land.lakpoi) {
    if (snl > 0) {
      t_grnd = frac_sno_eff * t_soisno[nlevsno - snl] + (1.0 - frac_sno_eff - frac_h2osfc) * t_soisno[nlevsno] +
               frac_h2osfc * t_h2osfc;
    } else {
      t_grnd = (1.0 - frac_h2osfc) * t_soisno[nlevsno] + frac_h2osfc * t_h2osfc;
    }
  }
} // CalculateGroundTemp

template <class ArrayD1>
void CalculateSoilAlpha(const LandType &Land, const double &frac_sno, const double &frac_h2osfc, const double &smpmin,
                        const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 dz,
                        const ArrayD1 t_soisno, const ArrayD1 watsat, const ArrayD1 sucsat,
                        const ArrayD1 bsw, const ArrayD1 watdry, const ArrayD1 watopt,
                        const ArrayD1 rootfr_road_perv, ArrayD1 rootr_road_perv, double &qred, double &hr,
                        double &soilalpha, double &soilalpha_u) {
  qred = 1.0; // soil surface relative humidity

  if (!Land.lakpoi) {

    double hr_road_perv; // relative humidity for urban pervious road
    double fac;          // soil wetness of surface layer
    double wx;           // partial volume of ice and water of surface layer
    double psit;         // negative potential of soil
    double eff_porosity; //  effective porosity in layer
    double vol_ice;      //  partial volume of ice lens in layer
    double vol_liq;      //  partial volume of liquid water in layer

    if (Land.ctype == icol_road_perv) {
      hr_road_perv = 0.0;
    }
    if (Land.ltype != istwet && Land.ltype != istice && Land.ltype != istice_mec) {
      if (Land.ltype == istsoil || Land.ltype == istcrop) {
        wx = (h2osoi_liq[nlevsno] / denh2o + h2osoi_ice[nlevsno] / denice) / dz[nlevsno];
        fac = std::min(1.0, wx / watsat[0]);
        fac = std::max(fac, 0.01);
        psit = -sucsat[0] * pow(fac, (-bsw[0]));
        psit = std::max(smpmin, psit);
        // modify qred to account for h2osfc
        hr = exp(psit / roverg / t_soisno[nlevsno]);
        qred = (1.0 - frac_sno - frac_h2osfc) * hr + frac_sno + frac_h2osfc;
        soilalpha = qred;
      } else if (Land.ctype == icol_road_perv) {

        // Pervious road depends on water in total soil column
        for (int j = 0; j < nlevbed; j++) {
          if (t_soisno[j + nlevsno] >= tfrz) {
            vol_ice = std::min(watsat[j], h2osoi_ice[j + nlevsno] / (dz[j + nlevsno] * denice));
            eff_porosity = watsat[j] - vol_ice;
            vol_liq = std::min(eff_porosity, h2osoi_liq[j + nlevsno] / (dz[j + nlevsno] * denh2o));
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
          for (int j = 0; j < nlevsoi; j++) {
            rootr_road_perv[j] = rootr_road_perv[j] / hr_road_perv;
          }
        }
        soilalpha_u = qred;
      } else if (Land.ctype == icol_sunwall || Land.ctype == icol_shadewall) {
        qred = 0.0;
        soilalpha_u = spval;
      } else if (Land.ctype == icol_roof || Land.ctype == icol_road_imperv) {
        qred = 1.0;
        soilalpha_u = spval;
      }
    } else {
      soilalpha = spval;
    }
  }
} // CalculateSoilAlpha

template <class ArrayD1>
void CalculateSoilBeta(const LandType &Land, const double &frac_sno, const double &frac_h2osfc,
                       const ArrayD1 watsat, const ArrayD1 watfc, const ArrayD1 h2osoi_liq,
                       const ArrayD1 h2osoi_ice, const ArrayD1 dz, double &soilbeta) {

  calc_soilevap_stress(Land, frac_sno, frac_h2osfc, watsat, watfc, h2osoi_liq, h2osoi_ice, dz, soilbeta);
} // CalculateSoilBeta()

template <class ArrayD1>
void CalculateHumidities(const LandType &Land, const int &snl, const double &forc_q, const double &forc_pbot,
                         const double &t_h2osfc, const double &t_grnd, const double &frac_sno,
                         const double &frac_sno_eff, const double &frac_h2osfc, const double &qred, const double &hr,
                         const ArrayD1 t_soisno, double &qg_snow, double &qg_soil, double &qg, double &qg_h2osfc,
                         double &dqgdT) {
  if (!Land.lakpoi) {

    double eg;      // water vapor pressure at temperature T [pa]
    double qsatg;   // saturated humidity [kg/kg]
    double degdT;   // d(eg)/dT
    double qsatgdT; // d(qsatg)/dT

    if (Land.ltype == istsoil || Land.ltype == istcrop) {
      QSat(t_soisno[nlevsno - snl], forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_snow = qsatg;
      dqgdT = frac_sno * qsatgdT;
      QSat(t_soisno[nlevsno], forc_pbot, eg, degdT, qsatg, qsatgdT);
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
} // CalculateHumidities

template <class ArrayD1>
void GroundProperties(const LandType &Land, const int &snl, const double &frac_sno, const double &forc_th,
                      const double &forc_q, const double &elai, const double &esai, const double &htop,
                      const ArrayD1 displar, const ArrayD1 z0mr, const ArrayD1 h2osoi_liq,
                      const ArrayD1 h2osoi_ice, double &emg, double &emv, double &htvp, double &z0mg, double &z0hg,
                      double &z0qg, double &z0mv, double &z0hv, double &z0qv, double &beta, double &zii, double &thv,
                      double &z0m, double &displa, double &cgrnd, double &cgrnds, double &cgrndl) {
  if (!Land.lakpoi) {
    double avmuir; // ir inverse optical depth per unit leaf area

    // Ground emissivity - only calculate for non-urban landunits
    // Urban emissivities are currently read in from data file
    if (!Land.urbpoi) {
      if (Land.ltype == istice || Land.ltype == istice_mec) {
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
    htvp = hvap;
    if (h2osoi_liq[nlevsno - snl] <= 00 && h2osoi_ice[nlevsno - snl] > 0.0) {
      htvp = hsub;
    }

    // Ground roughness lengths over non-lake columns (includes bare ground, ground
    // underneath canopy, wetlands, etc.)
    if (frac_sno > 0.0) {
      z0mg = zsno;
    } else {
      z0mg = zlnd;
    }
    z0hg = z0mg; // initial set only
    z0qg = z0mg; // initial set only
    z0m = z0mr[Land.vtype] * htop;
    displa = displar[Land.vtype] * htop;

    // vegetation roughness lengths
    z0mv = z0m;
    z0hv = z0mv;
    z0qv = z0mv;

    // Potential, virtual potential temperature, and wind speed at the reference height
    beta = 1.0;
    zii = 1000.0;
    thv = forc_th * (1.0 + 0.61 * forc_q);

    // Initial set for calculation
    cgrnd = 0.0;
    cgrnds = 0.0;
    cgrndl = 0.0;
  }
} // GroundProperties

} // namespace ELM