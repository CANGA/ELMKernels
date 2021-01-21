/* functions derived from CanopyTemperatureMod.F90

*/
#include <algorithm>
#include <cmath>
#include "clm_constants.h"
#include "qsat.h"
#include "landtype.h"

namespace ELM {

class CanopyTemperature {

private:
  double qred; // soil surface relative humidity
  double hr; // relative humidity

public:

/* CanopyTemperature::SaveGroundTemp()
DESCRIPTION: Record t_h2osfc and t_soisno prior to updating.

INPUTS:
Land                       [LandType] struct containing information about landtype
t_h2osfc                   [double] surface water temperature 
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)

OUTPUTS:
t_h2osfc_bef               [double] saved surface water temperature         
tssbef[nlevgrnd+nlevsno]   [double] soil/snow temperature before update    
*/
  void SaveGroundTemp(
    const LandType& Land,
    const double& t_h2osfc,
    const double* t_soisno,
    double& t_h2osfc_bef,
    double* tssbef) {

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


/* CanopyTemperature::CalculateGroundTemp()
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
  void CalculateGroundTemp(
    const LandType& Land,
    const int& snl,
    const double& frac_sno_eff, 
    const double& frac_h2osfc,
    const double& t_h2osfc,
    const double* t_soisno,
    double& t_grnd) { 
    // ground temperature is weighted average of exposed soil, snow, and h2osfc
    if (!Land.lakpoi) {
      if (snl > 0) {
        t_grnd = frac_sno_eff * t_soisno[nlevsno - snl] + (1.0 - frac_sno_eff - frac_h2osfc) * 
        t_soisno[nlevsno] + frac_h2osfc * t_h2osfc;
      } else {
        t_grnd = (1.0 - frac_h2osfc) * t_soisno[nlevsno] + frac_h2osfc * t_h2osfc;
      }
    }
  } // CalculateGroundTemp


/* CanopyTemperature::CalculateSoilAlpha()
DESCRIPTION: Calculate soilalpha factor that reduces ground saturated specific humidity. 
It looks like soilalpha doesn't get used in maint-1.0 branch, but both qred and hr do.

INPUTS:
Land                       [LandType] struct containing information about landtype
nlevbed                    [int] number of layers to bedrock
nlevsoi                    [int] number of hydrologically active soil layers
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

  void CalculateSoilAlpha(
    const LandType& Land,
    const int& nlevbed,
    const int& nlevsoi,
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
    const double* rootfr_road_perv,
  
    double* rootr_road_perv,
    double& qred,
    double& hr,
    double& soilalpha,
    double& soilalpha_u)
  {
    qred = 1.0; // soil surface relative humidity
  
    if (!Land.lakpoi) {
  
      double hr_road_perv; // relative humidity for urban pervious road
      double fac; // soil wetness of surface layer
      double wx; // partial volume of ice and water of surface layer
      double psit; // negative potential of soil
      double eff_porosity; //  effective porosity in layer
      double vol_ice;      //  partial volume of ice lens in layer
      double vol_liq;      //  partial volume of liquid water in layer
  
      if (Land.ctype == icol_road_perv) { hr_road_perv = 0.0; }
      if (Land.ltype != istwet && Land.ltype != istice && Land.ltype != istice_mec) {
        if (Land.ltype == istsoil || Land.ltype == istcrop) {
          wx   = (h2osoi_liq[nlevsno] / denh2o + h2osoi_ice[nlevsno] / denice) / dz[nlevsno];
          fac  = std::min(1.0, wx / watsat[0]);
          fac  = std::max(fac, 0.01);
          psit = -sucsat[0] * pow(fac, (-bsw[0]));
          psit = std::max(smpmin, psit);
          // modify qred to account for h2osfc
          hr   = exp(psit / roverg / t_soisno[nlevsno]);
          qred = (1.0 - frac_sno - frac_h2osfc) * hr + frac_sno + frac_h2osfc;
          soilalpha = qred;
        } else if (Land.ctype == icol_road_perv) {
  
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

/* CanopyTemperature::CalculateHumidities()
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
  void CalculateHumidities(
    const LandType& Land,
    const int& snl,
    const double& forc_q,
    const double& forc_pbot,
    const double& t_h2osfc,
    const double& t_grnd,
    const double& frac_sno,
    const double& frac_sno_eff,
    const double& frac_h2osfc,
    const double& qred,
    const double& hr,
    const double* t_soisno,
  
    double& qg_snow,
    double& qg_soil,
    double& qg,
    double& qg_h2osfc,
    double& dqgdT)
  {
    if (!Land.lakpoi) {
  
      double eg;           // water vapor pressure at temperature T [pa]
      double qsatg;        // saturated humidity [kg/kg]
      double degdT;        // d(eg)/dT
      double qsatgdT;      // d(qsatg)/dT
  
      if (Land.ltype == istsoil || Land.ltype == istcrop) {
        QSat(t_soisno[nlevsno-snl], forc_pbot, eg, degdT, qsatg, qsatgdT);
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



/* CanopyTemperature::GroundProperties()
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
  void GroundProperties(
    const LandType& Land,
    const int& snl,
    const double& frac_sno,
    const double& forc_th,
    const double& forc_q,
    const double& elai,
    const double& esai,
    const double& htop,
    const double* displar,
    const double* z0mr,
    const double* h2osoi_liq,
    const double* h2osoi_ice,
  
    double& emg,
    double& emv,
    double& htvp,
    double& z0mg,
    double& z0hg,
    double& z0qg,
    double& z0mv,
    double& z0hv,
    double& z0qv,
    double& beta,
    double& zii,
    double& thv,
    double& z0m,
    double& displa,
    double& cgrnd,
    double& cgrnds,
    double& cgrndl)
  {
    if (!Land.lakpoi) {
      double avmuir; //ir inverse optical depth per unit leaf area
  
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
      z0m = z0mr[Land.vtype] * htop;
      displa = displar[Land.vtype] * htop;
  
      // vegetation roughness lengths
      z0mv = z0m;
      z0hv = z0mv;
      z0qv = z0mv;
  
      // Potential, virtual potential temperature, and wind speed at the reference height
      beta = 1.0;
      zii  = 1000.0;
      thv  = forc_th * (1.0 + 0.61 * forc_q);
  
      //Initial set for calculation
      cgrnd = 0.0;
      cgrnds = 0.0;
      cgrndl = 0.0;
    }
  } // GroundProperties


/* CanopyTemperature::CalculateForcingHeight()
DESCRIPTION: Calculate roughness length, displacement height, and forcing height for wind , temp, and humidity.

INPUTS:
Land             [LandType] struct containing information about landtype
veg_active       [bool] true => landunit is an urban point
frac_veg_nosno   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
forc_hgt_u       [double] observational height of wind [m]         
forc_hgt_t       [double] observational height of temperature [m]          
forc_hgt_q       [double] observational height of specific humidity [m]
z0m              [double] momentum roughness length (m)
z0mg             [double] roughness length over ground, momentum [m]
z_0_town         [double] momentum roughness length of urban landunit (
z_d_town         [double] displacement height of urban landunit (m)
forc_t           [double] atmospheric temperature (Kelvin)         
displa           [double] displacement height (m)

OUTPUTS:
forc_hgt_u_patch [double] observational height of wind at pft level [m]
forc_hgt_t_patch [double] observational height of temperature at pft level [m]
forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
thm              [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
*/
  void CalculateForcingHeight(
    const LandType& Land,
    const bool& veg_active,
    const int& frac_veg_nosno,
    const double& forc_hgt_u,
    const double& forc_hgt_t,
    const double& forc_hgt_q,
    const double& z0m,
    const double& z0mg,
    const double& z_0_town,
    const double& z_d_town,
    const double& forc_t,
    const double& displa,
  
    double& forc_hgt_u_patch,
    double& forc_hgt_t_patch,
    double& forc_hgt_q_patch,
    double& thm)
  {
    // Make forcing height a pft-level quantity that is the atmospheric forcing 
    // height plus each pft's z0m+displa
    if (veg_active) {
      if (Land.ltype == istsoil || Land.ltype == istcrop) {
        if (frac_veg_nosno == 0) {
          forc_hgt_u_patch = forc_hgt_u + z0mg + displa;
          forc_hgt_t_patch = forc_hgt_t + z0mg + displa;
          forc_hgt_q_patch = forc_hgt_q + z0mg + displa;
        } else {
          forc_hgt_u_patch = forc_hgt_u + z0m + displa;
          forc_hgt_t_patch = forc_hgt_t + z0m + displa;
          forc_hgt_q_patch = forc_hgt_q + z0m + displa;
        }
      } else if (Land.ltype == istwet || Land.ltype == istice || Land.ltype == istice_mec) {
        forc_hgt_u_patch = forc_hgt_u + z0mg;
        forc_hgt_t_patch = forc_hgt_t + z0mg;
        forc_hgt_q_patch = forc_hgt_q + z0mg;
      } else if (Land.ltype == istdlak) {
        forc_hgt_u_patch = forc_hgt_u;
        forc_hgt_t_patch = forc_hgt_t;
        forc_hgt_q_patch = forc_hgt_q;
      } else if (Land.urbpoi) {
        forc_hgt_u_patch = forc_hgt_u + z_0_town + z_d_town;
        forc_hgt_t_patch = forc_hgt_t + z_0_town + z_d_town;
        forc_hgt_q_patch = forc_hgt_q + z_0_town + z_d_town;
      }
    }
  
    thm  = forc_t + 0.0098 * forc_hgt_t_patch;
  } // CalculateForcingHeight


/* CanopyTemperature::InitializeEnergyFluxes()
DESCRIPTION: Set energy flux terms to 0.0 before calculation.

INPUTS:
Land             [LandType] struct containing information about landtype

OUTPUTS:
eflx_sh_tot      [double] total sensible heat flux (W/m**2) [+ to atm]
eflx_sh_tot_u    [double] urban total sensible heat flux (W/m**2) [+ to atm]
eflx_sh_tot_r    [double] rural total sensible heat flux (W/m**2) [+ to atm]
eflx_lh_tot      [double] total latent heat flux (W/m**2)  [+ to atm]
eflx_lh_tot_u    [double] urban total latent heat flux (W/m**2)  [+ to atm]
eflx_lh_tot_r    [double] rural total latent heat flux (W/m**2)  [+ to atm]
eflx_sh_veg      [double] sensible heat flux from leaves (W/m**2) [+ to atm]
qflx_evap_tot    [double] qflx_evap_soi + qflx_evap_can + qflx_tran_veg
qflx_evap_veg    [double] vegetation evaporation (mm H2O/s) (+ = to atm)
qflx_tran_veg    [double] vegetation transpiration (mm H2O/s) (+ = to atm)
*/
  void InitializeEnergyFluxes(
    const LandType& Land,
    double& eflx_sh_tot,
    double& eflx_sh_tot_u,
    double& eflx_sh_tot_r,
    double& eflx_lh_tot,
    double& eflx_lh_tot_u,
    double& eflx_lh_tot_r,
    double& eflx_sh_veg,
    double& qflx_evap_tot,
    double& qflx_evap_veg,
    double& qflx_tran_veg)
  {
    // Initial set (needed for history tape fields)
    eflx_sh_tot = 0.0;
    if (Land.urbpoi) {
      eflx_sh_tot_u = 0.0;
    } else if (Land.ltype == istsoil || Land.ltype == istcrop) {
      eflx_sh_tot_r = 0.0;
    }
    eflx_lh_tot = 0.0;
    if (Land.urbpoi) {
      eflx_lh_tot_u = 0.0;
    } else if (Land.ltype == istsoil || Land.ltype == istcrop) {
      eflx_lh_tot_r = 0.0;
    }
    eflx_sh_veg = 0.0;
    qflx_evap_tot = 0.0;
    qflx_evap_veg = 0.0;
    qflx_tran_veg = 0.0;
  } // InitializeEnergyFluxes

};

} // namespace ELM
