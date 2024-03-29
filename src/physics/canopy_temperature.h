/*! \file canopy_temperature.h
\brief Functions derived from CanopyTemperatureMod.F90

Calculates ground temperature, emissivity, humidity, roughness lengths, forcing height, etc.

Call sequence:
old_ground_temp() -> ground_temp() -> calc_soilalpha() -> calc_soilbeta() -> humidities() ->
ground_properties() -> forcing_height() -> init_energy_fluxes()
*/
#pragma once

#include "elm_constants.h"
#include "land_data.h"
#include "qsat.h"
#include "surface_resistance.h"

#include <algorithm>
#include <cmath>

#include "compile_options.hh"

namespace ELM::canopy_temperature {

/*! Record t_h2osfc and t_soisno prior to updating.

\param[in]  Land                       [LandType] struct containing information about landtype
\param[in]  t_h2osfc                   [double] surface water temperature
\param[in]  t_soisno[nlevgrnd()+nlevsno()] [double] col soil temperature (Kelvin)
\param[out] t_h2osfc_bef               [double] saved surface water temperature
\param[out] tssbef[nlevgrnd()+nlevsno()]   [double] soil/snow temperature before update
*/
template <typename ArrayD1>
ACCELERATE
void old_ground_temp(const LandType& Land, const double& t_h2osfc, const ArrayD1 t_soisno, double& t_h2osfc_bef,
                     ArrayD1 tssbef);

/*! Calculate average ground temp.

\param[in]  Land                       [LandType] struct containing information about landtype
\param[in]  snl                        [int] number of snow layers
\param[in]  frac_sno_eff               [double] eff. fraction of ground covered by snow (0 to 1)
\param[in]  frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
\param[in]  t_h2osfc                   [double] surface water temperature
\param[in]  t_soisno[nlevgrnd()+nlevsno()] [double] col soil temperature (Kelvin)
\param[out] t_grnd                     [double] ground temperature (Kelvin)
*/
template <typename ArrayD1>
ACCELERATE
void ground_temp(const LandType& Land, const int& snl, const double& frac_sno_eff, const double& frac_h2osfc,
                 const double& t_h2osfc, const ArrayD1 t_soisno, double& t_grnd);

/*! Calculate soilalpha factor that reduces ground saturated specific humidity.
It looks like soilalpha doesn't get used in maint-1.2 branch, but both qred and hr do.

\param[in]  Land                         [LandType] struct containing information about landtype
\param[in]  frac_sno                     [double] fraction of ground covered by snow (0 to 1)
\param[in]  frac_h2osfc                  [double] fraction of ground covered by surface water (0 to 1)
\param[in]  h2osoi_liq[nlevgrnd()+nlevsno()] [double] liquid water (kg/m2)
\param[in]  h2osoi_ice[nlevgrnd()+nlevsno()] [double] ice lens (kg/m2)
\param[in]  dz[nlevgrnd()+nlevsno()]         [double] layer thickness (m)
\param[in]  t_soisno[nlevgrnd()+nlevsno()]   [double] col soil temperature (Kelvin)
\param[in]  watsat[nlevgrnd()]             [double] volumetric soil water at saturation (porosity)
\param[in]  sucsat[nlevgrnd()]             [double] minimum soil suction (mm)
\param[in]  bsw[nlevgrnd()]                [double] Clapp and Hornberger "b"
\param[in]  watdry[nlevgrnd()]             [double] btran parameter for btran = 0
\param[in]  watopt[nlevgrnd()]             [double] btran parameter for btran = 1
\param[out] qred                         [double] soil surface relative humidity
\param[out] hr                           [double] relative humidity
\param[out] soilalpha                    [double] factor that reduces ground saturated specific humidity (-)
*/
template <typename ArrayD1>
ACCELERATE
void calc_soilalpha(const LandType& Land, const double& frac_sno, const double& frac_h2osfc,
                    const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 dz, const ArrayD1 t_soisno,
                    const ArrayD1 watsat, const ArrayD1 sucsat, const ArrayD1 bsw, const ArrayD1 watdry,
                    const ArrayD1 watopt, double& qred, double& hr, double& soilalpha);

/*! Calculate soilbeta parameter.

\param[in]  Land                         [LandType] struct containing information about landtypet
\param[in]  frac_sno                     [double] fraction of ground covered by snow (0 to 1)
\param[in]  frac_h2osfc                  [double] fraction of ground covered by surface water (0 to 1)
\param[in]  watsat[nlevgrnd()]             [double] volumetric soil water at saturation (porosity)
\param[in]  watfc[nlevgrnd()]              [double] volumetric soil water at field capacity
\param[in]  h2osoi_liq[nlevgrnd()+nlevsno()] [double] liquid water (kg/m2)
\param[in]  h2osoi_ice[nlevgrnd()+nlevsno()] [double] ice lens (kg/m2)
\param[in]  dz[nlevgrnd()+nlevsno()]         [double] layer thickness (m)
\param[out] soilbeta                     [double] factor that reduces ground evaporation
*/
template <typename ArrayD1>
ACCELERATE
void calc_soilbeta(const LandType& Land, const double& frac_sno, const double& frac_h2osfc, const ArrayD1 watsat,
                   const ArrayD1 watfc, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 dz,
                   double& soilbeta);

/*! Calculate saturated vapor pressure, specific humidity and their derivatives at ground surface.
Compute humidities individually for snow, soil, h2osfc for vegetated landunits.

\param[in]  Land                       [LandType] struct containing information about landtype
\param[in]  snl                        [int] number of snow layers
\param[in]  forc_q                     [double] atmospheric specific humidity (kg/kg)
\param[in]  forc_pbot                  [double] atmospheric pressure (Pa)
\param[in]  t_h2osfc                   [double] surface water temperature
\param[in]  t_grnd                     [double] ground temperature (Kelvin)
\param[in]  frac_sno                   [double] fraction of ground covered by snow (0 to 1)
\param[in]  frac_sno_eff               [double] eff. fraction of ground covered by snow (0 to 1)
\param[in]  frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
\param[in]  qred                       [double] soil surface relative humidity
\param[in]  hr                         [double] relative humidity
\param[in]  t_soisno[nlevgrnd()+nlevsno()] [double] soil temperature (Kelvin)
\param[out] qg_snow                    [double] specific humidity at snow surface [kg/kg]
\param[out] qg_soil                    [double] specific humidity at soil surface [kg/kg]
\param[out] qg                         [double] ground specific humidity [kg/kg]
\param[out] qg_h2osfc                  [double] specific humidity at h2osfc surface [kg/kg]
\param[out] dqgdT                      [double] d(qg)/dT
*/
template <typename ArrayD1>
ACCELERATE
void humidities(const LandType& Land, const int& snl, const double& forc_q, const double& forc_pbot,
                const double& t_h2osfc, const double& t_grnd, const double& frac_sno, const double& frac_sno_eff,
                const double& frac_h2osfc, const double& qred, const double& hr, const ArrayD1 t_soisno,
                double& qg_snow, double& qg_soil, double& qg, double& qg_h2osfc, double& dqgdT);

/*! Calculate ground emissivity, latent heat constant, roughness lengths,
potential temp and wind speed.

\param[in]  Land                         [LandType] struct containing information about landtype
\param[in]  snl                          [int] number of snow layers
\param[in]  frac_sno                     [double] fraction of ground covered by snow (0 to 1)
\param[in]  forc_th                      [double] atmospheric potential temperature (Kelvin)
\param[in]  forc_q                       [double] atmospheric specific humidity (kg/kg)
\param[in]  elai                         [double] one-sided leaf area index with burying by snow
\param[in]  esai                         [double] one-sided stem area index with burying by snow
\param[in]  htop                         [double] canopy top (m)
\param[in]  displar[numpft()]              [double] ratio of displacement height to canopy top height (-)
\param[in]  z0mr[numpft()]                 [double] ratio of momentum roughness length to canopy top height (-)
\param[in]  h2osoi_liq[nlevgrnd()+nlevsno()] [double] liquid water (kg/m2)
\param[in]  h2osoi_ice[nlevgrnd()+nlevsno()] [double] ice lens (kg/m2)
\param[out] emg                          [double] ground emissivity
\param[out] emv                          [double] vegetation emissivity
\param[out] htvp                         [double] latent heat of vapor of water (or sublimation) [j/kg]
\param[out] z0mg                         [double] roughness length over ground, momentum [m]
\param[out] z0hg                         [double] roughness length over ground, sensible heat [m]
\param[out] z0qg                         [double] roughness length over ground, latent heat [m]
\param[out] z0mv                         [double] roughness length over vegetation, momentum [m]
\param[out] z0hv                         [double] roughness length over vegetation, sensible heat [m]
\param[out] z0qv                         [double] roughness length over vegetation, latent heat [m]
\param[out] thv                          [double] virtual potential temperature (kelvin)
\param[out] z0m                          [double] momentum roughness length (m)
\param[out] displa                       [double] displacement height (m)
*/
template <typename ArrayD1, typename SubviewD1>
ACCELERATE
void ground_properties(const LandType& Land, const int& snl, const double& frac_sno, const double& forc_th,
                       const double& forc_q, const double& elai, const double& esai, const double& htop,
                       const SubviewD1 displar, const SubviewD1 z0mr, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice,
                       double& emg, double& emv, double& htvp, double& z0mg, double& z0hg, double& z0qg, double& z0mv,
                       double& z0hv, double& z0qv, double& thv, double& z0m, double& displa);

/*! Calculate roughness length, displacement height, and forcing height for wind, temp, and humidity.

\param[in]  Land             [LandType] struct containing information about landtype
\param[in]  veg_active       [bool] true => cell contains active vegetation
\param[in]  frac_veg_nosno   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  z0m              [double] momentum roughness length (m)
\param[in]  z0mg             [double] roughness length over ground, momentum [m]
\param[in]  forc_t           [double] atmospheric temperature (Kelvin)
\param[in]  displa           [double] displacement height (m)
\param[out] forc_hgt_u_patch [double] observational height of wind at pft level [m]
\param[out] forc_hgt_t_patch [double] observational height of temperature at pft level [m]
\param[out] forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
\param[out] thm              [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
*/
ACCELERATE
void forcing_height(const LandType& Land, const bool& veg_active, const int& frac_veg_nosno,
                    const double& z0m, const double& z0mg, const double& forc_t, const double& displa,
                    double& forc_hgt_u_patch, double& forc_hgt_t_patch,
                    double& forc_hgt_q_patch, double& thm);

/*! Set energy flux terms to 0.0 before calculation.

\param[in]  Land             [LandType] struct containing information about landtype
\param[out] eflx_sh_tot      [double] total sensible heat flux (W/m**2) [+ to atm]
\param[out] eflx_lh_tot      [double] total latent heat flux (W/m**2)  [+ to atm]
\param[out] eflx_sh_veg      [double] sensible heat flux from leaves (W/m**2) [+ to atm]
\param[out] qflx_evap_tot    [double] qflx_evap_soi + qflx_evap_can + qflx_tran_veg
\param[out] qflx_evap_veg    [double] vegetation evaporation (mm H2O/s) (+ = to atm)
\param[out] qflx_tran_veg    [double] vegetation transpiration (mm H2O/s) (+ = to atm)
*/
ACCELERATE
void init_energy_fluxes(const LandType& Land, double& eflx_sh_tot, double& eflx_lh_tot,
                       double& eflx_sh_veg, double& qflx_evap_tot, double& qflx_evap_veg,
                       double& qflx_tran_veg);

} // namespace ELM::canopy_temperature

#include "canopy_temperature_impl.hh"
