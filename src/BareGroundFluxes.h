/* functions derived from BareGroundFluxesMod.F90
Compute sensible and latent fluxes and their derivatives with respect
to ground temperature using ground temperatures from previous time step.

Note:
the "displa" displacement height variable used here is local to these functions.
It is not the same as the displa calculated in CanopyTemperature 
*/

#pragma once

#include <algorithm>
#include <cmath>
#include "clm_constants.h"
#include "landtype.h"

namespace ELM {

class BareGroundFluxes {

private:
  double zldis;   // reference height "minus" zero displacement height [m]
  double displa;  // displacement height [m]
  double dth;     // diff of virtual temp. between ref. height and surface
  double dqh;     // diff of humidity between ref. height and surface
  double obu;     // Monin-Obukhov length (m)
  double ur;      // wind speed at reference height [m/s]
  double um;      // wind speed including the stablity effect [m/s]
  double temp1;   // relation for potential temperature profile  
  double temp2;   // relation for specific humidity profile
  double temp12m; // relation for potential temperature profile applied at 2-m
  double temp22m; // relation for specific humidity profile applied at 2-m
  double ustar;   // friction velocity [m/s]

  // local vars
  double frac_veg_nosno;
  double forc_q;
  double forc_th;
  double forc_hgt_u_patch;
  double thm;
  double thv;
  double z0mg;


public:

/* BareGroundFluxes::InitializeFlux
DESCRIPTION:
initialize fluxes and call MoninObukIni()

INPUTS:
Land                [LandType] struct containing information about landtype 
frac_veg_nosno_in   [int] fraction of vegetation not covered by snow (0 OR 1) [-]
forc_u              [double] atmospheric wind speed in east direction (m/s)
forc_v              [double] atmospheric wind speed in north direction (m/s)
forc_q_in           [double] atmospheric specific humidity (kg/kg)
forc_th_in          [double] atmospheric potential temperature (Kelvin)
forc_hgt_u_patch_in [double] observational height of wind at pft level [m]
thm_in              [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
thv_in              [double] virtual potential temperature (kelvin)
t_grnd              [double] ground temperature (Kelvin)
qg                  [double] ground specific humidity [kg/kg]         
z0mg_in             [double] roughness length over ground, momentum [m]

OUTPUTS:
dlrad            [double] downward longwave radiation below the canopy [W/m2]
ulrad            [double] upward longwave radiation above the canopy [W/m2]  
*/
  void InitializeFlux(
    const LandType& Land,
    const int& frac_veg_nosno_in,
    const double& forc_u,
    const double& forc_v,
    const double& forc_q_in,
    const double& forc_th_in,
    const double& forc_hgt_u_patch_in,
    const double& thm_in,
    const double& thv_in,
    const double& t_grnd,
    const double& qg,
    const double& z0mg_in,
    double& dlrad,
    double& ulrad);


/* BareGroundFluxes::StabilityIteration()
DESCRIPTION:
calculate Monin-Obukhov length and wind speed

INPUTS:
Land             [LandType] struct containing information about landtype
forc_hgt_t_patch [double] observational height of temperature at pft level [m]
forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
z0mg             [double] roughness length over ground, momentum [m]
zii              [double] convective boundary height [m]
beta             [double] coefficient of convective velocity [-]

OUTPUTS:
z0hg             [double] roughness length over ground, sensible heat [m]
z0qg             [double] roughness length over ground, latent heat [m]
*/
  void StabilityIteration(
    const LandType& Land,
    const double& forc_hgt_t_patch,
    const double& forc_hgt_q_patch,
    const double& z0mg,
    const double& zii,
    const double& beta,
    double& z0hg,
    double& z0qg);


/* BareGroundFluxes::ComputeFlux
DESCRIPTION:
calculated bare ground water and energy fluxes

INPUTS:
Land                       [LandType] struct containing information about landtype
snl                        [int] number of snow layers
forc_rho                   [double] density (kg/m**3)
soilbeta                   [double] soil wetness relative to field capacity
dqgdT                      [double] temperature derivative of "qg"
htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]
t_h2osfc                   [double] surface water temperature
qg_snow                    [double] specific humidity at snow surface [kg/kg]
qg_soil                    [double] specific humidity at soil surface [kg/kg]
qg_h2osfc                  [double] specific humidity at h2osfc surface [kg/kg]
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
forc_pbot                  [double]  atmospheric pressure (Pa)

OUTPUTS:
cgrnds                     [double] deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
cgrndl                     [double] deriv of soil latent heat flux wrt soil temp [w/m**2/k]
cgrnd                      [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
eflx_sh_grnd               [double] sensible heat flux from ground (W/m**2) [+ to atm]
eflx_sh_tot                [double] total sensible heat flux (W/m**2) [+ to atm]
eflx_sh_snow               [double] sensible heat flux from snow (W/m**2) [+ to atm]
eflx_sh_soil               [double] sensible heat flux from soil (W/m**2) [+ to atm]
eflx_sh_h2osfc             [double] sensible heat flux from soil (W/m**2) [+ to atm]
qflx_evap_soi              [double] soil evaporation (mm H2O/s) (+ = to atm)
qflx_evap_tot              [double] qflx_evap_soi + qflx_evap_can + qflx_tran_veg
qflx_ev_snow               [double] evaporation flux from snow (W/m**2) [+ to atm]  
qflx_ev_soil               [double] evaporation flux from soil (W/m**2) [+ to atm]  
qflx_ev_h2osfc             [double] evaporation flux from h2osfc (W/m**2) [+ to atm]
t_ref2m                    [double]  2 m height surface air temperature (Kelvin)
t_ref2m_r                  [double]  Rural 2 m height surface air temperature (Kelvin)
q_ref2m                    [double]  2 m height surface specific humidity (kg/kg)
rh_ref2m_r                 [double]  Rural 2 m height surface relative humidity (%)
rh_ref2m                   [double]  2 m height surface relative humidity (%)
*/
  template<class dArray_type>
  void ComputeFlux(
    const LandType& Land,
    const int& snl,
    const double& forc_rho,
    const double& soilbeta,
    const double& dqgdT,
    const double& htvp,
    const double& t_h2osfc,
    const double& qg_snow,
    const double& qg_soil,
    const double& qg_h2osfc,
    const dArray_type t_soisno,
    const double& forc_pbot,
    double& cgrnds,
    double& cgrndl,
    double& cgrnd,
    double& eflx_sh_grnd,
    double& eflx_sh_tot,
    double& eflx_sh_snow,
    double& eflx_sh_soil,
    double& eflx_sh_h2osfc,
    double& qflx_evap_soi,
    double& qflx_evap_tot,
    double& qflx_ev_snow,
    double& qflx_ev_soil,
    double& qflx_ev_h2osfc,
    double& t_ref2m,
    double& t_ref2m_r,
    double& q_ref2m,
    double& rh_ref2m,
    double& rh_ref2m_r);

};

} // namespace


#include "BareGroundFluxes_impl.h"

