/*! \file friction_velocity.h
\brief Internal functions derived from FrictionVelocityMod.F90

Split into 5 functions - wind, temperature, humidity, 2-m temp, 2-m humidity.
The friction velocity scheme is based on the work of Zeng et al. (1998):
Intercomparison of bulk aerodynamic algorithms for the computation
of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
Vol. 11, 2628-2644.
*/

#pragma once

namespace ELM::friction_velocity {

/*! Initialization of the Monin-Obukhov length. The scheme is based on the work of
Zeng et al. (1998): Intercomparison of bulk aerodynamic algorithms for the computation of sea
surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11, 2628-2644. (internal)

\param[in]  ur    [double]  wind speed at reference height [m/s]
\param[in]  thv   [double]  virtual potential temperature (kelvin)
\param[in]  dthv  [double]  diff of vir. poten. temp. between ref. height and surface
\param[in]  zldis [double]  reference height "minus" zero displacement height [m]
\param[in]  z0m   [double]  roughness length, momentum [m]
\param[out] um    [double]  wind speed including the stability effect [m/s]
\param[out] obu   [double]  monin-obukhov length (m)
*/
void monin_obukhov_length(const double& ur, const double& thv, const double& dthv, const double& zldis,
                          const double& z0m, double& um, double& obu);

/*! Calculation of the friction velocity of surface boundary layer. (internal)

forc_hgt_u_patch   [double] observational height of wind at pft level [m]
\param[in]  displa [double] displacement height (m)
\param[in]  um     [double] wind speed including the stability effect [m/s]
\param[in]  obu    [double] monin-obukhov length (m)
\param[in]  z0m    [double] roughness length over vegetation, momentum [m]
\param[out] ustar  [double] friction velocity [m/s]
*/
void friction_velocity_wind(const double& forc_hgt_u_patch, const double& displa, const double& um, const double& obu,
                            const double& z0m, double& ustar);

/*! Calculation of the relation for potential temperature of surface boundary layer. (internal)

\param[in]  forc_hgt_t_patch [double] observational height of temperature at pft level [m]
\param[in]  displa           [double] displacement height (m)
\param[in]  obu              [double] monin-obukhov length (m)
\param[in]  z0h              [double] roughness length over vegetation, sensible heat [m]
\param[out] temp1            [double] relation for potential temperature profile
*/
void friction_velocity_temp(const double& forc_hgt_t_patch, const double& displa, const double& obu, const double& z0h,
                            double& temp1);

/*! Calculation of the relation for potential humidity of surface boundary layer. (internal)
\param[in]  forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
\param[in]  forc_hgt_t_patch [double] observational height of temperature at pft level
\param[in]  displa           [double] displacement height (m)
\param[in]  obu              [double] monin-obukhov length (m)
\param[in]  z0h              [double] roughness length over vegetation, sensible heat [m]
\param[in]  z0q              [double] roughness length over vegetation, latent heat [m]
\param[in]  temp1            [double] relation for potential temperature profile
\param[our] temp2            [double] relation for specific humidity profile
*/
void friction_velocity_humidity(const double& forc_hgt_q_patch, const double& forc_hgt_t_patch, const double& displa,
                                const double& obu, const double& z0h, const double& z0q, const double& temp1,
                                double& temp2);

/*! Calculation of the relation for potential temperature at 2-m. (internal)
\param[in]  obu     [double] monin-obukhov length (m)
\param[in]  z0h     [double] roughness length over vegetation, sensible heat [m]
\param[out] temp12m [double] relation for potential temperature profile applied at 2-m
*/
void friction_velocity_temp2m(const double& obu, const double& z0h, double& temp12m);

/*! Calculation of the relation for potential humidity at 2-m. (internal)
INPUTS:
\param[in]  obu     [double] monin-obukhov length (m)
\param[in]  z0h     [double] roughness length over vegetation, sensible heat [m]
\param[in]  z0q     [double] roughness length over vegetation, latent heat [m]
\param[out] temp22m [double] relation for specific humidity profile applied at 2-m
*/
void friction_velocity_humidity2m(const double& obu, const double& z0h, const double& z0q, const double& temp12m,
                                  double& temp22m);

} // namespace ELM::friction_velocity
