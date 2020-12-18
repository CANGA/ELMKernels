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
  const double* t_soisno,

  double* tssbef,
  double* t_h2osfc_bef)
{
  if (!lakpoi) {
    for (int i = 0; i < nlevgrnd + nlevsno; i++) {
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
{ // ground temperature is weighted average of exposed soil, snow, and h2osfc
  if (!lakpoi) {
    if (snl > 0) {
      t_grnd = frac_sno_eff * t_soisno[nlevsno - snl] + (1.0 - frac_sno_eff - frac_h2osfc) * 
      t_soisno[nlevsno] + frac_h2osfc * t_h2osfc;
    } else {
      t_grnd = (1.0 - frac_h2osfc) * t_soisno[nlevsno] + frac_h2osfc * t_h2osfc;
    }
  }
}


/* CalculateSoilAlpha()
DESCRIPTION: Calculate soilalpha factor that reduces ground saturated specific humidity. 
It looks like soilalpha doesn't get used in maint-1.0 branch.

INPUTS:
smpmin [double] restriction for min of soil potential (mm)
h2osoi_liq[nlengrnd+nlevsno-1] [double] liquid water (kg/m2)
h2osoi_ice[nlengrnd+nlevsno-1] [double] ice lens (kg/m2)
t_soisno[nlengrnd+nlevsno-1] [double]
dz[nlengrnd+nlevsno-1] [double] layer thickness (m) 
watsat[nlevgrnd]  [double]  volumetric soil water at saturation (porosity)
bsw[nlevgrnd]            [double] Clapp and Hornberger "b" 
sucsat[nlevgrnd] [double]  minimum soil suction (mm)
watdry[nlevgrnd] [double]  btran parameter for btran = 0
watopt[nlevgrnd] [double]  btran parameter for btran = 1
rootr_road_perv[nlevgrnd]  [double] effective fraction of roots in each soil layer for urban pervious road
rootfr_road_perv[nlevgrnd] [double] fraction of roots in each soil layer for urban pervious road


*/

void CalculateSoilAlpha(
  const bool& lakpoi,
  const int& ctype,
  const int& ltype,
  const int& nlevbed,
  const int& nlevsoi,
  const int& nlevsno,
  const double& frac_sno,
  const double& frac_h2osfc,
  const double& smpmin,

  const double* h2osoi_liq,
  const double* h2osoi_ice,
  const double* dz,
  const double* t_soisno,

  const double* watsat,
  const double* sucsat,
  const double* bsw,
  const double* watdry,
  const double* watopt,
  const double* rootr_road_perv,
  const double* rootfr_road_perv,

  double& qred,
  double& soilalpha,
  double& soilalpha_u)
{
  qred = 1.0; // soil surface relative humidity

  if (!lakpoi) {

    double hr_road_perv; // relative humidity for urban pervious road
    double fac; // soil wetness of surface layer
    double wx; // partial volume of ice and water of surface layer
    double psit; // negative potential of soil
    double hr; // relative humidity
    double eff_porosity; //  effective porosity in layer
    double vol_ice;      //  partial volume of ice lens in layer
    double vol_liq;      //  partial volume of liquid water in layer

    if (ctype == icol_road_perv) { hr_road_perv = 0.0; }
    if (ltype != istwet && ltype != istice && ltype != istice_mec) {
      if (ltype == istsoil || ltype == istcrop) {
        wx   = (h2osoi_liq[nlevsno] / denh2o + h2osoi_ice[nlevsno] / denice) / dz[nlevsno];
        fac  = std::min(1.0, wx / watsat[0]);
        fac  = std::max(fac, 0.01);
        psit = -sucsat[0] * pow(fac, (-bsw[0]));
        psit = std::max(smpmin, psit);
        // modify qred to account for h2osfc
        hr   = exp(psit / roverg / t_soisno[nlevsno]);
        qred = (1.0 - frac_sno - frac_h2osfc) * hr + frac_sno + frac_h2osfc;
        soilalpha = qred;
      } else if (ctype == icol_road_perv) {

        // Pervious road depends on water in total soil column
        for (int j = 0; j < nlevbed; j++) {
          if (t_soisno[j+nlevsno] >= tfrz) {
            vol_ice = std::min(watsat[j], h2osoi_ice[j+nlevsno] / (dz[j+nlevsno] * denice));
            eff_porosity = watsat[j] - vol_ice;
            vol_liq = std::min(eff_porosity, h2osoi_liq[j+nlevsno] / (dz[j+nlevsno] * denh2o));
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
DESCRIPTION:
Saturated vapor pressure, specific humidity and their derivatives at ground surface. 
Compute humidities individually for snow, soil, h2osfc for vegetated landunits.

INPUTS:
forc_q atmospheric specific humidity (kg/kg)    
forc_pbot atmospheric pressure (Pa) 

frac_sno_eff eff. fraction of ground covered by snow (0 to 1)




qg_snow   specific humidity at snow surface [kg/kg]
qg_soil   specific humidity at soil surface [kg/kg]
qg        ground specific humidity [kg/kg]         
qg_h2osfc specific humidity at h2osfc surface [kg/kg]
dqgdT     d(qg)/dT                                 
*/
void CalculateHumidities(
  const bool& lakpoi,
  const int& ltype,
  const int& nlevsno,
  const int& snl,
  const double& forc_q,
  const double& forc_pbot,
  const double& t_h2osfc,
  const double& t_grnd,
  const double& frac_sno,
  const double& frac_sno_eff,
  const double& qred,

  const double* t_soisno,



  double& qg_snow,
  double& qg_soil,
  double& qg,
  double& qg_h2osfc,
  double& dqgdT,


  )
{
  if (!lakpoi) {

    double eg;           // water vapor pressure at temperature T [pa]
    double qsatg;        // saturated humidity [kg/kg]
    double degdT;        // d(eg)/dT
    double qsatgdT;      // d(qsatg)/dT

    if (ltype == istsoil || ltype == istcrop) {
      QSat(t_soisno[nlevsno-snl], forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_snow = qsatg;
      dqgdT = frac_sno * qsatgdT;
      QSat(t_soisno[nlevsno] , forc_pbot, eg, degdT, qsatg, qsatgdT);
      if (qsatg > forc_q && forc_q > hr * qsatg) {
        qsatg = forc_q;
        qsatgdT = 0.0;
      }
      qg_soil = hr * qsatg;
      dqgdT = dqgdT + (1.0 - frac_sno - frac_h2osfc) * hr * qsatgdT;

      // to be consistent with hs_top values in SoilTemp, set qg_snow to qg_soil for snl = 0 case
      // this ensures hs_top_snow will equal hs_top_soil
      if (snl == 0) {
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



/* GroundProperties()
DESCRIPTION:
Calculate ground emissivity, latent heat constant, roughness lengths, 
potential temp and wind speed.

*/
void GroundProperties(
  const bool& lakpoi,
  const bool& urbpoi,
  const double& frac_sno,
  const double* h2osoi_liq,
  const double* h2osoi_ice,


  double& emg,
  double& htvp,
  double& z0hg,
  double& z0mg,
  double& z0qg,
  double& beta,
  double& zii,
  double& thv)
{
  if (!lakpoi) {

    // Ground emissivity - only calculate for non-urban landunits 
    // Urban emissivities are currently read in from data file
    if (!urbpoi) {
      if (ltype == istice || ltype == istice_mec) {
        emg = 0.97;
      } else {
        emg = (1.0 - frac_sno) * 0.96 + frac_sno * 0.97;
      }
    }

    // Latent heat. We arbitrarily assume that the sublimation occurs
    // only as h2osoi_liq = 0
    htvp = hvap;
    if (h2osoi_liq[nlevsno-snl] <= 00 && h2osoi_ice[nlevsno-snl] > 0.0) { htvp = hsub; }

    // Ground roughness lengths over non-lake columns (includes bare ground, ground
    // underneath canopy, wetlands, etc.)
    if (frac_sno > 0.0) {
      z0mg = zsno;
    } else {
      z0mg = zlnd;
    }
    z0hg = z0mg;            // initial set only
    z0qg = z0mg;            // initial set only

    // Potential, virtual potential temperature, and wind speed at the reference height
    beta = 1.0;
    zii  = 1000.0;
    thv  = forc_th * (1.0 + 0.61 * forc_q);
  }
}
