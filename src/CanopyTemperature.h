/* functions derived from CanopyTemperatureMod.F90

*/
#pragma once

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
                    dArray_type tssbef);

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
                         const double &t_h2osfc, const dArray_type t_soisno, double &t_grnd);

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
                        double &soilalpha, double &soilalpha_u);

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
                       const dArray_type h2osoi_ice, const dArray_type dz, double &soilbeta);

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
                         double &dqgdT);

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
                      double &z0m, double &displa, double &cgrnd, double &cgrnds, double &cgrndl);

/* CalculateForcingHeight()
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
void CalculateForcingHeight(const LandType &Land, const bool &veg_active, const int &frac_veg_nosno,
                            const double &forc_hgt_u, const double &forc_hgt_t, const double &forc_hgt_q,
                            const double &z0m, const double &z0mg, const double &z_0_town, const double &z_d_town,
                            const double &forc_t, const double &displa, double &forc_hgt_u_patch,
                            double &forc_hgt_t_patch, double &forc_hgt_q_patch, double &thm);

/* InitializeEnergyFluxes()
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
void InitializeEnergyFluxes(const LandType &Land, double &eflx_sh_tot, double &eflx_sh_tot_u, double &eflx_sh_tot_r,
                            double &eflx_lh_tot, double &eflx_lh_tot_u, double &eflx_lh_tot_r, double &eflx_sh_veg,
                            double &qflx_evap_tot, double &qflx_evap_veg, double &qflx_tran_veg);

} // namespace ELM

#include "CanopyTemperature_impl.hh"
