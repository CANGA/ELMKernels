/* functions derived from CanopyTemperatureMod.F90

*/
#pragma once

#include "QSat.h"
#include "SurfaceResistance.h"
#include "clm_constants.h"
#include "landtype.h"
#include <algorithm>
#include <cmath>

namespace ELM {

/* SaveGroundTemp()
DESCRIPTION: Record t_h2osfc and t_soisno prior to updating.

INPUTS:
Land                       [LandType] struct containing information about landtype
t_h2osfc                   [double] surface water temperature
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)

OUTPUTS:
t_h2osfc_bef               [double] saved surface water temperature
tssbef[nlevgrnd+nlevsno]   [double] soil/snow temperature before update
*/
template <class dArray_type>
void SaveGroundTemp(const LandType &Land, const double &t_h2osfc, const dArray_type t_soisno, double &t_h2osfc_bef,
                    dArray_type tssbef) {

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

/* CalculateGroundTemp()
DESCRIPTION: Calculate average ground temp.

INPUTS:
Land                       [LandType] struct containing information about landtype
snl                        [int] number of snow layers
frac_sno_eff               [double] eff. fraction of ground covered by snow (0 to 1)
frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
t_h2osfc                   [double] surface water temperature
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)

OUTPUTS:
t_grnd                     [double] ground temperature (Kelvin)
*/
template <class dArray_type>
void CalculateGroundTemp(const LandType &Land, const int &snl, const double &frac_sno_eff, const double &frac_h2osfc,
                         const double &t_h2osfc, const dArray_type t_soisno, double &t_grnd) {
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

/* CalculateSoilAlpha()
DESCRIPTION: Calculate soilalpha factor that reduces ground saturated specific humidity.
It looks like soilalpha doesn't get used in maint-1.0 branch, but both qred and hr do.

INPUTS:
Land                       [LandType] struct containing information about landtype
frac_sno                   [double] fraction of ground covered by snow (0 to 1)
frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
smpmin [double] restriction for min of soil potential (mm)
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)
t_soisno[nlevgrnd+nlevsno]   [double] col soil temperature (Kelvin)
watsat[nlevgrnd]             [double] volumetric soil water at saturation (porosity)
sucsat[nlevgrnd]             [double] minimum soil suction (mm)
bsw[nlevgrnd]                [double] Clapp and Hornberger "b"
watdry[nlevgrnd]             [double] btran parameter for btran = 0
watopt[nlevgrnd]             [double] btran parameter for btran = 1
rootfr_road_perv[nlevgrnd]   [double] fraction of roots in each soil layer for urban pervious road

OUTPUTS:
rootr_road_perv[nlevgrnd]    [double] effective fraction of roots in each soil layer for urban pervious road
qred                         [double] soil surface relative humidity
hr                           [double] relative humidity
soilalpha                    [double] factor that reduces ground saturated specific humidity (-)
soilalpha_u                  [double] Urban factor that reduces ground saturated specific humidity (-)
*/
template <class dArray_type>
void CalculateSoilAlpha(const LandType &Land, const double &frac_sno, const double &frac_h2osfc, const double &smpmin,
                        const dArray_type h2osoi_liq, const dArray_type h2osoi_ice, const dArray_type dz,
                        const dArray_type t_soisno, const dArray_type watsat, const dArray_type sucsat,
                        const dArray_type bsw, const dArray_type watdry, const dArray_type watopt,
                        const dArray_type rootfr_road_perv, dArray_type rootr_road_perv, double &qred, double &hr,
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

/* CalculateSoilBeta()
INPUTS:
Land                         [LandType] struct containing information about landtypet
frac_sno                     [double] fraction of ground covered by snow (0 to 1)
frac_h2osfc                  [double] fraction of ground covered by surface water (0 to 1)
watsat[nlevgrnd]             [double] volumetric soil water at saturation (porosity)
watfc[nlevgrnd]              [double] volumetric soil water at field capacity
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
soilbeta [double] factor that reduces ground evaporation
*/
template <class dArray_type>
void CalculateSoilBeta(const LandType &Land, const double &frac_sno, const double &frac_h2osfc,
                       const dArray_type watsat, const dArray_type watfc, const dArray_type h2osoi_liq,
                       const dArray_type h2osoi_ice, const dArray_type dz, double &soilbeta) {

  calc_soilevap_stress(Land, frac_sno, frac_h2osfc, watsat, watfc, h2osoi_liq, h2osoi_ice, dz, soilbeta);
} // CalculateSoilBeta()

/* CalculateHumidities()
DESCRIPTION: Saturated vapor pressure, specific humidity and their derivatives at ground surface.
Compute humidities individually for snow, soil, h2osfc for vegetated landunits.

INPUTS:
Land                       [LandType] struct containing information about landtype
snl                        [int] number of snow layers
forc_q                     [double] atmospheric specific humidity (kg/kg)
forc_pbot                  [double] atmospheric pressure (Pa)
t_h2osfc                   [double] surface water temperature
t_grnd                     [double] ground temperature (Kelvin)
frac_sno                   [double] fraction of ground covered by snow (0 to 1)
frac_sno_eff               [double] eff. fraction of ground covered by snow (0 to 1)
frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
qred                       [double] soil surface relative humidity
hr                         [double] relative humidity
t_soisno[nlevgrnd+nlevsno] [double] soil temperature (Kelvin)

OUTPUTS:
qg_snow                    [double] specific humidity at snow surface [kg/kg]
qg_soil                    [double] specific humidity at soil surface [kg/kg]
qg                         [double] ground specific humidity [kg/kg]
qg_h2osfc                  [double] specific humidity at h2osfc surface [kg/kg]
dqgdT                      [double] d(qg)/dT
*/
template <class dArray_type>
void CalculateHumidities(const LandType &Land, const int &snl, const double &forc_q, const double &forc_pbot,
                         const double &t_h2osfc, const double &t_grnd, const double &frac_sno,
                         const double &frac_sno_eff, const double &frac_h2osfc, const double &qred, const double &hr,
                         const dArray_type t_soisno, double &qg_snow, double &qg_soil, double &qg, double &qg_h2osfc,
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

/* GroundProperties()
DESCRIPTION: Calculate ground emissivity, latent heat constant, roughness lengths,
potential temp and wind speed.

INPUTS:
Land                         [LandType] struct containing information about landtype
snl                          [int] number of snow layers
frac_sno                     [double] fraction of ground covered by snow (0 to 1)
forc_th                      [double] atmospheric potential temperature (Kelvin)
forc_q                       [double] atmospheric specific humidity (kg/kg)
elai                         [double] one-sided leaf area index with burying by snow
esai                         [double] one-sided stem area index with burying by snow
htop                         [double] canopy top (m)
displar[numpft]              [double] ratio of displacement height to canopy top height (-)
z0mr[numpft]                 [double] ratio of momentum roughness length to canopy top height (-)
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)

OUTPUTS:
emg                          [double] ground emissivity
emv                          [double] vegetation emissivity
htvp                         [double] latent heat of vapor of water (or sublimation) [j/kg]
z0mg                         [double] roughness length over ground, momentum [m]
z0hg                         [double] roughness length over ground, sensible heat [m]
z0qg                         [double] roughness length over ground, latent heat [m]
z0mv                         [double] roughness length over vegetation, momentum [m]
z0hv                         [double] roughness length over vegetation, sensible heat [m]
z0qv                         [double] roughness length over vegetation, latent heat [m]
beta                         [double] coefficient of convective velocity [-]
zii                          [double] convective boundary height [m]
thv                          [double] virtual potential temperature (kelvin)
z0m                          [double] momentum roughness length (m)
displa                       [double] displacement height (m)
cgrnd                        [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
cgrnds                       [double] deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
cgrndl                       [double] deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
*/
template <class dArray_type>
void GroundProperties(const LandType &Land, const int &snl, const double &frac_sno, const double &forc_th,
                      const double &forc_q, const double &elai, const double &esai, const double &htop,
                      const dArray_type displar, const dArray_type z0mr, const dArray_type h2osoi_liq,
                      const dArray_type h2osoi_ice, double &emg, double &emv, double &htvp, double &z0mg, double &z0hg,
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