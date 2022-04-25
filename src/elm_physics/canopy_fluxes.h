/*! \file canopy_fluxes.h
\brief Functions derived from CanopyFluxesMod.F90

Calculate leaf temperature, leaf fluxes, transpiration, photosynthesis and update the dew accumulation due to
evaporation. Also calculates the irrigation rate (not currently active). Irrigation() can be called anytime after
InitializeFlux()

Call sequence: initialize_flux() -> stability_iteration() -> compute_flux()
*/
#pragma once

#include "elm_constants.h"
#include "friction_velocity.h"
#include "land_data.h"
#include "pft_data.h"
#include "photosynthesis.h"
#include "qsat.h"
#include "soil_moist_stress.h"
#include <algorithm>
#include <assert.h>
#include <cmath>

#include "kokkos_includes.hh"

namespace ELM::canopy_fluxes {

/*! Initialize variables for photosynthesis and call monin_obukhov_length() for vegetated cells.

\param[in]  Land                            [LandType] struct containing information about landtype
\param[in]  snl                             [int] number of snow layers
\param[in]  frac_veg_nosno                  [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  frac_sno                        [double] fraction of ground covered by snow (0 to 1)
\param[in]  forc_hgt_u_patch                [double] observational height of wind at pft level [m]
\param[in]  thm                             [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
\param[in]  thv                             [double] virtual potential temperature (kelvin)
\param[in]  max_dayl                        [double] maximum daylength for this grid cell (s)
\param[in]  dayl                            [double] daylength (seconds)
\param[in]  altmax_indx                     [int] index corresponding to maximum active layer depth from current year
\param[in]  altmax_lastyear_indx            [int] index corresponding to maximum active layer depth from prior year
\param[in]  t_soisno[nlevgrnd+nlevsno]      [double] col soil temperature (Kelvin)
\param[in]  h2osoi_ice[nlevgrnd+nlevsno]    [double] ice lens (kg/m2)
\param[in]  h2osoi_liq[nlevgrnd+nlevsno]    [double] liquid water (kg/m2)
\param[in]  dz[nlevgrnd+nlevsno]            [double] layer thickness (m)
\param[in]  rootfr[nlevgrnd]                [double] fraction of roots in each soil layer
\param[in]  tc_stress                       [double] critical soil temperature for soil water stress (C)
\param[in]  sucsat[nlevgrnd]                [double] minimum soil suction (mm)
\param[in]  watsat[nlevgrnd]                [double] volumetric soil water at saturation (porosity)
\param[in]  bsw[nlevgrnd]                   [double] Clapp and Hornberger "b
\param[in]  smpso                           [double] soil water potential at full stomatal opening (mm)
\param[in]  smpsc                           [double] soil water potential at full stomatal closure (mm)
\param[in]  elai                            [double] one-sided leaf area index with burying by snow
\param[in]  esai                            [double] one-sided stem area index with burying by snow
\param[in]  emv                             [double] vegetation emissivity
\param[in]  emg                             [double] ground emissivity
\param[in]  qg                              [double] ground specific humidity [kg/kg]
\param[in]  t_grnd                          [double] ground temperature (Kelvin)
\param[in]  forc_t                          [double] atmospheric temperature (Kelvin)
\param[in]  forc_pbot                       [double] atmospheric pressure (Pa)
\param[in]  forc_lwrad                      [double] downward infrared (longwave) radiation (W/m**2)
\param[in]  forc_u                          [double] atmospheric wind speed in east direction (m/s)
\param[in]  forc_v                          [double] atmospheric wind speed in north direction (m/s)
\param[in]  forc_q                          [double] atmospheric specific humidity (kg/kg)
\param[in]  forc_th                         [double] atmospheric potential temperature (Kelvin)
\param[in]  z0mg                            [double] roughness length over ground, momentum [m]
\param[out] btran                           [double]  transpiration wetness factor (0 to 1)
\param[out] z0mv                            [double] roughness length over vegetation, momentum [m]
\param[out] z0hv                            [double] roughness length over vegetation, sensible heat [m]
\param[out] z0qv                            [double] roughness length over vegetation, latent heat [m]
\param[out] displa                          [double]  displacement height (m)
\param[out] rootr[nlevgrnd]                 [double] effective fraction of roots in each soil layer
\param[out] eff_porosity[nlevgrnd]          [double] effective soil porosity
\param[out] dayl_factor                     [double] scalar (0-1) for daylength effect on Vcmax
\param[out] air                             [double] atmos. radiation temporay set
\param[out] bir                             [double] atmos. radiation temporay set
\param[out] cir                             [double] atmos. radiation temporay set
\param[out] el                              [double] vapor pressure on leaf surface [pa]
\param[out] qsatl                           [double] leaf specific humidity [kg/kg]
\param[out] qsatldT                         [double] derivative of "qsatl" on "t_veg"
\param[out] taf                             [double] air temperature within canopy space [K]
\param[out] qaf                             [double] humidity of canopy air [kg/kg]
\param[out] um                              [double] wind speed including the stablity effect [m/s]
\param[out] ur                              [double] wind speed at reference height [m/s]
\param[out] obu                             [double] Monin-Obukhov length (m)
\param[out] zldis                           [double] reference height "minus" zero displacement height [m]
\param[out] delq                            [double] temporary
\param[out] t_veg                           [double]  vegetation temperature (Kelvin)
*/
template <class ArrayD1>
ACCELERATED
void initialize_flux(const LandType& Land, const int& snl, const int& frac_veg_nosno, const double& frac_sno,
                     const double& forc_hgt_u_patch, const double& thm, const double& thv, const double& max_dayl,
                     const double& dayl, const int& altmax_indx, const int& altmax_lastyear_indx,
                     const ArrayD1 t_soisno, const ArrayD1 h2osoi_ice, const ArrayD1 h2osoi_liq, const ArrayD1 dz,
                     const ArrayD1 rootfr, const double& tc_stress, const ArrayD1 sucsat, const ArrayD1 watsat,
                     const ArrayD1 bsw, const double& smpso, const double& smpsc, const double& elai,
                     const double& esai, const double& emv, const double& emg, const double& qg, const double& t_grnd,
                     const double& forc_t, const double& forc_pbot, const double& forc_lwrad, const double& forc_u,
                     const double& forc_v, const double& forc_q, const double& forc_th, const double& z0mg,
                     double& btran, double& displa, double& z0mv, double& z0hv, double& z0qv, ArrayD1 rootr,
                     ArrayD1 eff_porosity, double& dayl_factor, double& air, double& bir, double& cir, double& el,
                     double& qsatl, double& qsatldT, double& taf, double& qaf, double& um, double& ur, double& obu,
                     double& zldis, double& delq, double& t_veg);

/*! Calculate Monin-Obukhov length and wind speed, call photosynthesis, calculate ET & SH flux
Iterates until convergence, up to 40 iterations, calling friction velocity functions, then
photosynthesis for both sun & shade.

\param[in]  Land                       [LandType] struct containing information about landtype
\param[in]  dtime                      [double] timestep size (sec)
\param[in]  snl                        [int] number of snow layers
\param[in]  frac_veg_nosno             [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  frac_sno                   [double] fraction of ground covered by snow (0 to 1)
\param[in]  forc_hgt_u_patch           [double] observational height of wind at pft level [m]
\param[in]  forc_hgt_t_patch           [double] observational height of temperature at pft level [m]
\param[in]  forc_hgt_q_patch           [double] observational height of specific humidity at pft level [m]
\param[in]  fwet                       [double] fraction of canopy that is wet (0 to 1)
\param[in]  fdry                       [double] fraction of foliage that is green and dry [-]
\param[in]  laisun                     [double] sunlit leaf area
\param[in]  laisha                     [double] shaded leaf area
\param[in]  forc_rho                   [double] air density (kg/m**3)
\param[in]  snow_depth                 [double] snow height (m)
\param[in]  soilbeta                   [double] soil wetness relative to field capacity
\param[in]  frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
\param[in]  t_h2osfc                   [double] surface water temperature
\param[in]  sabv                       [double] solar radiation absorbed by vegetation (W/m**2)
\param[in]  h2ocan                     [double] canopy water (mm H2O)
\param[in]  htop                       [double] canopy top(m)
\param[in]  t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
\param[in]  air                        [double] atmos. radiation temporay set
\param[in]  bir                        [double] atmos. radiation temporay set
\param[in]  cir                        [double] atmos. radiation temporay set
\param[in]  ur                         [double] wind speed at reference height [m/s]
\param[in]  zldis                      [double] reference height "minus" zero displacement height [m]
\param[in]  displa                     [double]  displacement height (m)
\param[in]  elai                       [double] one-sided leaf area index with burying by snow
\param[in]  esai                       [double] one-sided stem area index with burying by snow
\param[in]  t_grnd                     [double] ground temperature (Kelvin)
\param[in]  forc_pbot                  [double]  atmospheric pressure (Pa)
\param[in]  forc_q                     [double] atmospheric specific humidity (kg/kg)
\param[in]  forc_th                    [double] atmospheric potential temperature (Kelvin)
\param[in]  z0mg                       [double] roughness length over ground, momentum [m]
\param[in]  z0mv                       [double] roughness length over vegetation, momentum [m]
\param[in]  z0hv                       [double] roughness length over vegetation, sensible heat [m]
\param[in]  z0qv                       [double] roughness length over vegetation, latent heat [m]
\param[in]  thm                        [double]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
\param[in]  thv                        [double] virtual potential temperature (kelvin)
\param[in]  qg                         [double] ground specific humidity [kg/kg]
\param[in]  psnveg                     [PSNVegData] constant vegetation data for current pft
\param[in]  nrad                       [int]  number of canopy layers above snow for radiative transfer
\param[in]  t10                        [double] 10-day running mean of the 2 m temperature (K)
\param[in]  tlai_z[nlevcan]            [double] pft total leaf area index for canopy layer
\param[in]  vcmaxcintsha               [double] leaf to canopy scaling coefficient - shade
\param[in]  vcmaxcintsun               [double] leaf to canopy scaling coefficient - sun
\param[in]  parsha_z[nlevcan]          [double] par absorbed per unit lai for canopy layer (w/m**2) - shade
\param[in]  parsun_z[nlevcan]          [double] par absorbed per unit lai for canopy layer (w/m**2) - sun
\param[in]  laisha_z[nlevcan]          [double] leaf area index for canopy layer, sunlit or shaded - shade
\param[in]  laisun_z[nlevcan]          [double] leaf area index for canopy layer, sunlit or shaded - sun
\param[in]  forc_pco2                  [double]  partial pressure co2 (Pa)
\param[in]  forc_po2                   [double]  partial pressure o2 (Pa))
\param[in]  dayl_factor                [double] scalar (0-1) for daylength effect on Vcmax
\param[out] btran                      [double]  transpiration wetness factor (0 to 1)
\param[out] qflx_tran_veg              [double] vegetation transpiration (mm H2O/s) (+ = to atm)
\param[out] qflx_evap_veg              [double] vegetation evaporation (mm H2O/s) (+ = to atm)
\param[out] eflx_sh_veg                [double] sensible heat flux from leaves (W/m**2) [+ to atm]
\param[out] wtg                        [double] heat conductance for ground [m/s]
\param[out] wtl0                       [double] normalized heat conductance for leaf [-]
\param[out] wta0                       [double] normalized heat conductance for air [-]
\param[out] wtal                       [double] normalized heat conductance for air and leaf [-]
\param[out] el                         [double] vapor pressure on leaf surface [pa]
\param[out] qsatl                      [double] leaf specific humidity [kg/kg]
\param[out] qsatldT                    [double] derivative of "qsatl" on "t_veg"
\param[out] taf                        [double] air temperature within canopy space [K]
\param[out] qaf                        [double] humidity of canopy air [kg/kg]
\param[out] um                         [double] wind speed including the stablity effect [m/s]
\param[out] dth                        [double] diff of virtual temp. between ref. height and surface
\param[out] dqh                        [double] diff of humidity between ref. height and surface
\param[out] obu                        [double] Monin-Obukhov length (m)
\param[out] temp1                      [double] relation for potential temperature profile
\param[out] temp2                      [double] relation for specific humidity profile
\param[out] temp12m                    [double] relation for potential temperature profile applied at 2-m
\param[out] temp22m                    [double] relation for specific humidity profile applied at 2-m
\param[out] tlbef                      [double] leaf temperature from previous iteration [K]
\param[out] delq                       [double] temporary
\param[out] dt_veg                     [double] change in t_veg, last iteration (Kelvin)
\param[out] t_veg                      [double]  vegetation temperature (Kelvin)
\param[out] wtgq                       [double] latent heat conductance for ground [m/s]
\param[out] wtalq                      [double] normalized latent heat cond. for air and leaf [-]
\param[out] wtlq0                      [double] normalized latent heat conductance for leaf [-]
\param[out] wtaq0                      [double] normalized latent heat conductance for air [-]
*/
template <class ArrayD1>
ACCELERATED
void stability_iteration(
    const LandType& Land, const double& dtime, const int& snl, const int& frac_veg_nosno, const double& frac_sno,
    const double& forc_hgt_u_patch, const double& forc_hgt_t_patch, const double& forc_hgt_q_patch, const double& fwet,
    const double& fdry, const double& laisun, const double& laisha, const double& forc_rho, const double& snow_depth,
    const double& soilbeta, const double& frac_h2osfc, const double& t_h2osfc, const double& sabv, const double& h2ocan,
    const double& htop, const ArrayD1 t_soisno, const double& air, const double& bir, const double& cir,
    const double& ur, const double& zldis, const double& displa, const double& elai, const double& esai,
    const double& t_grnd, const double& forc_pbot, const double& forc_q, const double& forc_th, const double& z0mg,
    const double& z0mv, const double& z0hv, const double& z0qv, const double& thm, const double& thv, const double& qg,
    const PSNVegData& psnveg, const int& nrad, const double& t10, const ArrayD1 tlai_z, const double& vcmaxcintsha,
    const double& vcmaxcintsun, const ArrayD1 parsha_z, const ArrayD1 parsun_z, const ArrayD1 laisha_z,
    const ArrayD1 laisun_z, const double& forc_pco2, const double& forc_po2, const double& dayl_factor, double& btran,
    double& qflx_tran_veg, double& qflx_evap_veg, double& eflx_sh_veg, double& wtg, double& wtl0, double& wta0,
    double& wtal, double& el, double& qsatl, double& qsatldT, double& taf, double& qaf, double& um, double& dth,
    double& dqh, double& obu, double& temp1, double& temp2, double& temp12m, double& temp22m, double& tlbef,
    double& delq, double& dt_veg, double& t_veg, double& wtgq, double& wtalq, double& wtlq0, double& wtaq0);

/*! Calculate water and energy fluxes for vegetated surfaces.

\param[in]  Land                       [LandType] struct containing information about landtype
\param[in]  dtime                      [double] timestep size (sec)
\param[in]  snl                        [int] number of snow layers
\param[in]  frac_veg_nosno             [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  frac_sno                   [double] fraction of ground covered by snow (0 to 1)
\param[in]  t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
\param[in]  frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
\param[in]  t_h2osfc                   [double] surface water temperature
\param[in]  sabv                       [double] solar radiation absorbed by vegetation (W/m**2)
\param[in]  qg_snow                    [double] specific humidity at snow surface [kg/kg]
\param[in]  qg_soil                    [double] specific humidity at soil surface [kg/kg]
\param[in]  qg_h2osfc                  [double] specific humidity at h2osfc [kg/kg]
\param[in]  dqgdT                      [double] d(qg)/dT
\param[in]  htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]
\param[in]  wtg                        [double] heat conductance for ground [m/s]
\param[in]  wtl0                       [double]  normalized heat conductance for leaf [-]
\param[in]  wta0                       [double] normalized heat conductance for air [-]
\param[in]  wtal                       [double] normalized heat conductance for air and leaf [-]
\param[in]  air                        [double] atmos. radiation temporay set
\param[in]  bir                        [double] atmos. radiation temporay set
\param[in]  cir                        [double] atmos. radiation temporay set
\param[in]  qsatl                      [double] leaf specific humidity [kg/kg]
\param[in]  qsatldT                    [double] derivative of "qsatl" on "t_ve
\param[in]  dth                        [double] diff of virtual temp. between ref. height and surface
\param[in]  dqh                        [double] diff of humidity between ref. height and surface
\param[in]  temp1                      [double] relation for potential temperature profile
\param[in]  temp2                      [double] relation for specific humidity profile
\param[in]  temp12m                    [double] relation for potential temperature profile applied at 2-m
\param[in]  temp22m                    [double] relation for specific humidity profile applied at 2-m
\param[in]  tlbef                      [double] leaf temperature from previous iteration [K]
\param[in]  delq                       [double] temporary
\param[in]  dt_veg                     [double] change in t_veg, last iteration (Kelvin)
\param[in]  t_veg                      [double]  vegetation temperature (Kelvin)
\param[in]  t_grnd                     [double] ground temperature (Kelvin)
\param[in]  forc_pbot                  [double]  atmospheric pressure (Pa)
\param[in]  qflx_tran_veg              [double] vegetation transpiration (mm H2O/s) (+ = to atm)
\param[in]  qflx_evap_veg              [double] vegetation evaporation (mm H2O/s) (+ = to atm)
\param[in]  eflx_sh_veg                [double] sensible heat flux from leaves (W/m**2) [+ to atm]
\param[in]  forc_q                     [double] atmospheric specific humidity (kg/kg)
\param[in]  forc_rho                   [double] air density (kg/m**3)
\param[in]  thm                        [double]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
\param[in]  emv                        [double] vegetation emissivity
\param[in]  emg                        [double] ground emissivity
\param[in]  forc_lwrad                 [double]  downward infrared (longwave) radiation (W/m**2)
\param[in]  wtgq                       [double] latent heat conductance for ground [m/s]
\param[in]  wtalq                      [double] normalized latent heat cond. for air and leaf [-]
\param[in]  wtlq0                      [double] normalized latent heat conductance for leaf [-]
\param[in]  wtaq0                      [double] normalized latent heat conductance for air [-]
\param[out] h2ocan                     [double] canopy water (mm H2O)
\param[out] eflx_sh_grnd               [double] sensible heat flux from ground (W/m**2) [+ to atm]
\param[out] eflx_sh_snow               [double] sensible heat flux from snow (W/m**2) [+ to atm]
\param[out] eflx_sh_soil               [double] sensible heat flux from soil (W/m**2) [+ to atm]
\param[out] eflx_sh_h2osfc             [double] sensible heat flux from h2osfc (W/m**2) [+ to atm]
\param[out] qflx_evap_soi              [double] soil evaporation (mm H2O/s) (+ = to atm)
\param[out] qflx_ev_snow               [double] evaporation flux from snow (W/m**2) [+ to atm]
\param[out] qflx_ev_soil               [double] evaporation flux from soil (W/m**2) [+ to atm]
\param[out] qflx_ev_h2osfc             [double] evaporation flux from h2osfc (W/m**2) [+ to atm]
\param[out] dlrad                      [double] downward longwave radiation below the canopy [W/m2]
\param[out] ulrad                      [double] upward longwave radiation above the canopy [W/m2]
\param[out] cgrnds                     [double] deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
\param[out] cgrndl                     [double] deriv of soil latent heat flux wrt soil temp [w/m**2/k]
\param[out] cgrnd                      [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
\param[out] t_ref2m                    [double]  2 m height surface air temperature (Kelvin)
\param[out] t_ref2m_r                  [double]  Rural 2 m height surface air temperature (Kelvin)
\param[out] q_ref2m                    [double]  2 m height surface specific humidity (kg/kg)
\param[out] rh_ref2m_r                 [double]  Rural 2 m height surface relative humidity (%)
\param[out] rh_ref2m                   [double]  2 m height surface relative humidity (%)
*/
template <class ArrayD1>
ACCELERATED
void compute_flux(const LandType& Land, const double& dtime, const int& snl, const int& frac_veg_nosno,
                  const double& frac_sno, const ArrayD1 t_soisno, const double& frac_h2osfc, const double& t_h2osfc,
                  const double& sabv, const double& qg_snow, const double& qg_soil, const double& qg_h2osfc,
                  const double& dqgdT, const double& htvp, const double& wtg, const double& wtl0, const double& wta0,
                  const double& wtal, const double& air, const double& bir, const double& cir, const double& qsatl,
                  const double& qsatldT, const double& dth, const double& dqh, const double& temp1, const double& temp2,
                  const double& temp12m, const double& temp22m, const double& tlbef, const double& delq,
                  const double& dt_veg, const double& t_veg, const double& t_grnd, const double& forc_pbot,
                  const double& qflx_tran_veg, const double& qflx_evap_veg, const double& eflx_sh_veg,
                  const double& forc_q, const double& forc_rho, const double& thm, const double& emv, const double& emg,
                  const double& forc_lwrad, const double& wtgq, const double& wtalq, const double& wtlq0,
                  const double& wtaq0, double& h2ocan, double& eflx_sh_grnd, double& eflx_sh_snow, double& eflx_sh_soil,
                  double& eflx_sh_h2osfc, double& qflx_evap_soi, double& qflx_ev_snow, double& qflx_ev_soil,
                  double& qflx_ev_h2osfc, double& dlrad, double& ulrad, double& cgrnds, double& cgrndl, double& cgrnd,
                  double& t_ref2m, double& t_ref2m_r, double& q_ref2m, double& rh_ref2m, double& rh_ref2m_r);

} // namespace ELM::canopy_fluxes

#include "canopy_fluxes_impl.hh"
