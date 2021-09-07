/*! \file BareGroundFluxes.h
\brief Functions derived from BareGroundFluxesMod.F90

Compute sensible and latent fluxes and their derivatives with respect
to ground temperature using ground temperatures from previous time step.

Call sequence: InitializeFlux_BG() -> StabilityIteration_BG() -> ComputeFlux_BG()

Note:
the "displa" displacement height variable used here is local to these functions.
It is not the same as the displa calculated in CanopyTemperature
*/

#pragma once

#include "ELMConstants.h"
#include "LandType.h"
#include <algorithm>
#include <cmath>

namespace ELM {

/*! Initialize variables and call MoninObukIni() for bare-ground cells.

\param[in]  Land             [LandType] struct containing information about landtype
\param[in]  frac_veg_nosno   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  forc_u           [double] atmospheric wind speed in east direction (m/s)
\param[in]  forc_v           [double] atmospheric wind speed in north direction (m/s)
\param[in]  forc_q           [double] atmospheric specific humidity (kg/kg)
\param[in]  forc_th          [double] atmospheric potential temperature (Kelvin)
\param[in]  forc_hgt_u_patch [double] observational height of wind at pft level [m]
\param[in]  thm              [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
\param[in]  thv              [double] virtual potential temperature (kelvin)
\param[in]  t_grnd           [double] ground temperature (Kelvin)
\param[in]  qg               [double] ground specific humidity [kg/kg]
\param[in]  z0mg             [double] roughness length over ground, momentum [m]
\param[out] dlrad            [double] downward longwave radiation below the canopy [W/m2]
\param[out] ulrad            [double] upward longwave radiation above the canopy [W/m2]
\param[out] zldis            [double] reference height "minus" zero displacement height [m]
\param[out] displa           [double] displacement height [m]
\param[out] dth              [double] diff of virtual temp. between ref. height and surface
\param[out] dqh              [double] diff of humidity between ref. height and surface
\param[out] obu              [double] Monin-Obukhov length (m)
\param[out] ur               [double] wind speed at reference height [m/s]
\param[out] um               [double] wind speed including the stablity effect [m/s]
*/
void InitializeFlux_BG(const LandType &Land, const int &frac_veg_nosno, const double &forc_u, const double &forc_v,
                       const double &forc_q, const double &forc_th, const double &forc_hgt_u_patch, const double &thm,
                       const double &thv, const double &t_grnd, const double &qg, const double &z0mg, double &dlrad,
                       double &ulrad, double &zldis, double &displa, double &dth, double &dqh, double &obu, double &ur,
                       double &um);

/*! Perform fixed 3-stage stability iteration over bare-ground cells.
Determine friction velocity and potential temperature and humidity
profiles of the surface boundary layer.

\param[in]  Land             [LandType] struct containing information about landtype
\param[in]  frac_veg_nosno   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  forc_hgt_t_patch [double] observational height of temperature at pft level [m]
\param[in]  forc_hgt_u_patch [double] observational height of wind at pft level [m]
\param[in]  forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
\param[in]  z0mg             [double] roughness length over ground, momentum [m]
\param[in]  zldis            [double] reference height "minus" zero displacement height [m]
\param[in]  displa           [double] displacement height [m]
\param[in]  dth              [double] diff of virtual temp. between ref. height and surface
\param[in]  dqh              [double] diff of humidity between ref. height and surface
\param[in]  ur               [double] wind speed at reference height [m/s]
\param[out] z0hg             [double] roughness length over ground, sensible heat [m]
\param[out] z0qg             [double] roughness length over ground, latent heat [m]
\param[out] obu              [double] Monin-Obukhov length (m)
\param[out] um               [double] wind speed including the stablity effect [m/s]
\param[out] temp1            [double] relation for potential temperature profile
\param[out] temp2            [double] relation for specific humidity profile
\param[out] temp12m          [double] relation for potential temperature profile applied at 2-m
\param[out] temp22m          [double] relation for specific humidity profile applied at 2-m
\param[out] ustar            [double] friction velocity [m/s]
*/
void StabilityIteration_BG(const LandType &Land, const int &frac_veg_nosno, const double &forc_hgt_t_patch,
                           const double &forc_hgt_u_patch, const double &forc_hgt_q_patch, const double &z0mg,
                           const double &zldis, const double &displa,
                           const double &dth, const double &dqh, const double &ur, const double &forc_q,
                           const double &forc_th, const double &thv, double &z0hg, double &z0qg, double &obu,
                           double &um, double &temp1, double &temp2, double &temp12m, double &temp22m, double &ustar);

/*! Calculate bare-ground water and energy fluxes.

\param[in]  Land                       [LandType] struct containing information about landtype
\param[in]  frac_veg_nosno             [int] fraction of vegetation not covered by snow (0 OR 1) [-]
\param[in]  snl                        [int] number of snow layers
\param[in]  forc_rho                   [double] density (kg/m**3)
\param[in]  soilbeta                   [double] soil wetness relative to field capacity
\param[in]  dqgdT                      [double] temperature derivative of "qg"
\param[in]  htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]
\param[in]  t_h2osfc                   [double] surface water temperature
\param[in]  qg_snow                    [double] specific humidity at snow surface [kg/kg]
\param[in]  qg_soil                    [double] specific humidity at soil surface [kg/kg]
\param[in]  qg_h2osfc                  [double] specific humidity at h2osfc surface [kg/kg]
\param[in]  t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
\param[in]  forc_pbot                  [double] atmospheric pressure (Pa)
\param[in]  dth                        [double] diff of virtual temp. between ref. height and surface
\param[in]  dqh                        [double] diff of humidity between ref. height and surface
\param[in]  temp1                      [double] relation for potential temperature profile
\param[in]  temp2                      [double] relation for specific humidity profile
\param[in]  temp12m                    [double] relation for potential temperature profile applied at 2-m
\param[in]  temp22m                    [double] relation for specific humidity profile applied at 2-m
\param[in]  ustar                      [double] friction velocity [m/s]
\param[in]  forc_q                     [double] atmospheric specific humidity (kg/kg)
\param[in]  thm                        [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
\param[out] cgrnds                     [double] deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
\param[out] cgrndl                     [double] deriv of soil latent heat flux wrt soil temp [w/m**2/k]
\param[out] cgrnd                      [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
\param[out] eflx_sh_grnd               [double] sensible heat flux from ground (W/m**2) [+ to atm]
\param[out] eflx_sh_tot                [double] total sensible heat flux (W/m**2) [+ to atm]
\param[out] eflx_sh_snow               [double] sensible heat flux from snow (W/m**2) [+ to atm]
\param[out] eflx_sh_soil               [double] sensible heat flux from soil (W/m**2) [+ to atm]
\param[out] eflx_sh_h2osfc             [double] sensible heat flux from soil (W/m**2) [+ to atm]
\param[out] qflx_evap_soi              [double] soil evaporation (mm H2O/s) (+ = to atm)
\param[out] qflx_evap_tot              [double] qflx_evap_soi + qflx_evap_can + qflx_tran_veg
\param[out] qflx_ev_snow               [double] evaporation flux from snow (W/m**2) [+ to atm]
\param[out] qflx_ev_soil               [double] evaporation flux from soil (W/m**2) [+ to atm]
\param[out] qflx_ev_h2osfc             [double] evaporation flux from h2osfc (W/m**2) [+ to atm]
\param[out] t_ref2m                    [double] 2 m height surface air temperature (Kelvin)
\param[out] t_ref2m_r                  [double] Rural 2 m height surface air temperature (Kelvin)
\param[out] q_ref2m                    [double] 2 m height surface specific humidity (kg/kg)
\param[out] rh_ref2m_r                 [double] Rural 2 m height surface relative humidity (%)
\param[out] rh_ref2m                   [double] 2 m height surface relative humidity (%)
*/
template <class ArrayD1>
void ComputeFlux_BG(const LandType &Land, const int &frac_veg_nosno, const int &snl, const double &forc_rho,
                    const double &soilbeta, const double &dqgdT, const double &htvp, const double &t_h2osfc,
                    const double &qg_snow, const double &qg_soil, const double &qg_h2osfc, const ArrayD1 t_soisno,
                    const double &forc_pbot, const double &dth, const double &dqh, const double &temp1,
                    const double &temp2, const double &temp12m, const double &temp22m, const double &ustar,
                    const double &forc_q, const double &thm, double &cgrnds, double &cgrndl, double &cgrnd,
                    double &eflx_sh_grnd, double &eflx_sh_tot, double &eflx_sh_snow, double &eflx_sh_soil,
                    double &eflx_sh_h2osfc, double &qflx_evap_soi, double &qflx_evap_tot, double &qflx_ev_snow,
                    double &qflx_ev_soil, double &qflx_ev_h2osfc, double &t_ref2m, double &t_ref2m_r, double &q_ref2m,
                    double &rh_ref2m, double &rh_ref2m_r);

} // namespace ELM

#include "BareGroundFluxes_impl.hh"
