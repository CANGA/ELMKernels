/*! \file canopy_hydrology.h
\brief Functions derived from CanopyHydrologyMod.F90

Compute interception, canopy fluxes to ground, fraction vegetation wet/dry, fraction of surface covered by H2O, FSCA,
initialize new snow layers.

Call sequence: interception() -> Irrigation() -> ground_flux() -> fraction_wet() -> snow_init() -> fraction_h2osfc()
*/
#pragma once

#include "elm_constants.h"
#include "land_data.h"
#include <algorithm>
#include <cmath>

#include "kokkos_includes.hh"

namespace ELM::canopy_hydrology {

/*! Calculate interception and partition incoming precipitation
into canopy storage, canopy runoff, rain and snow throughfall.
Also calculates fraction of precipitation that is rain/snow.

\param[in]     Land              [LandType] struct containing information about landtype
\param[in]     frac_veg_nosno    [int]      fraction of veg not covered by snow (0/1 now) [-]
\param[in]     forc_rain         [double]   rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
\param[in]     forc_snow         [double]   snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
\param[in]     dewmx             [double]   Maximum allowed dew [mm]
\param[in]     elai              [double]   one-sided leaf area index with burying by snow
\param[in]     esai              [double]   one-sided stem area index with burying by snow
\param[in]     dtime             [double]   time step length (sec)
\param[in,out] h2ocan            [double]   total canopy water (mm H2O)
\param[out]    qflx_candrip      [double]   rate of canopy runoff and snow falling off canopy [mm/s]
\param[out]    qflx_through_snow [double]   direct snow throughfall [mm/s]
\param[out]    qflx_through_rain [double]   direct rain throughfall [mm/s]
\param[out]    fracsnow          [double]   frac of precipitation that is snow [-]
\param[out]    fracrain          [double]   frac of precipitation that is rain [-]
*/
ACCELERATE
void interception(const LandType& Land, const int& frac_veg_nosno, const double& forc_rain, const double& forc_snow,
                  const double& dewmx, const double& elai, const double& esai, const double& dtime, double& h2ocan,
                  double& qflx_candrip, double& qflx_through_snow, double& qflx_through_rain, double& fracsnow,
                  double& fracrain);

/*! Determine whether we're irrigating here; set qflx_irrig appropriately.

\param[in]  Land               [LandType] struct containing information about landtype
\param[in]  irrig_rate         [double]   current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]
\param[out] n_irrig_steps_left [int] number of time steps for which we still need to irrigate today
\param[out] qflx_irrig         [double]   irrigation amount (mm/s)
*/
ACCELERATE
void Irrigation(const LandType& Land, const double& irrig_rate, int& n_irrig_steps_left, double& qflx_irrig);

/*! Add liquid and solid water inputs to ground surface after interception and
canopy storage losses.

\param[in]  Land              [LandType] struct containing information about landtype
\param[in]  do_capsnow        [bool]     true => do snow capping
\param[in]  frac_veg_nosno    [int]      fraction of veg not covered by snow (0/1 now) [-]
\param[in]  forc_rain         [double]   rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
\param[in]  forc_snow         [double]   snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
\param[in]  qflx_irrig        [double]   irrigation amount (mm/s)
\param[in]  qflx_candrip      [double]   rate of canopy runoff and snow falling off canopy [mm/s]
\param[in]  qflx_through_snow [double]   direct snow throughfall [mm/s]
\param[in]  qflx_through_rain [double]   direct rain throughfall [mm/s]
\param[in]  fracsnow          [double]   frac of precipitation that is snow [-]
\param[in]  fracrain          [double]   frac of precipitation that is rain [-]
\param[out] qflx_prec_grnd    [double]   water onto ground including canopy runoff [kg/(m2 s)]
\param[out] qflx_snwcp_liq    [double]   excess rainfall due to snow capping (mm H2O /s) [+]
\param[out] qflx_snwcp_ice    [double]   excess snowfall due to snow capping (mm H2O /s) [+]
\param[out] qflx_snow_grnd    [double]   snow on ground after interception (mm H2O/s) [+]
\param[out] qflx_rain_grnd    [double]   rain on ground after interception (mm H2O/s) [+]
*/
ACCELERATE
void ground_flux(const LandType& Land, const bool& do_capsnow, const int& frac_veg_nosno, const double& forc_rain,
                 const double& forc_snow, const double& qflx_irrig, const double& qflx_candrip,
                 const double& qflx_through_snow, const double& qflx_through_rain, const double& fracsnow,
                 const double& fracrain, double& qflx_prec_grnd, double& qflx_snwcp_liq, double& qflx_snwcp_ice,
                 double& qflx_snow_grnd, double& qflx_rain_grnd);

/*! Determine fraction of vegetated surfaces which are wet and
fraction of elai which is dry. The variable ``fwet'' is the
fraction of all vegetation surfaces which are wet including
stem area which contribute to evaporation. The variable ``fdry''
is the fraction of elai which is dry because only leaves
can transpire.  Adjusted for stem area which does not transpire.

\param[in]  Land           [LandType] struct containing information about landtype
\param[in]  frac_veg_nosno [int]      fraction of veg not covered by snow (0/1 now) [-]
\param[in]  dewmx          [double]   Maximum allowed dew [mm]
\param[in]  elai           [double]   one-sided leaf area index with burying by snow
\param[in]  esai           [double]   one-sided stem area index with burying by snow
\param[in]  h2ocan         [double]   total canopy water (mm H2O)
\param[out] fwet           [double]   fraction of canopy that is wet (0 to 1)
\param[out] fdry           [double]   fraction of foliage that is green and dry [-] (new)
*/
ACCELERATE
void fraction_wet(const LandType& Land, const int& frac_veg_nosno, const double& dewmx, const double& elai,
                  const double& esai, const double& h2ocan, double& fwet, double& fdry);

/*! Initialize new snow layer if the snow accumulation exceeds 10 mm, compute fractional SCA.

\param[in]     Land                          [LandType] struct containing information about landtype
\param[in]     dtime                         [double] time step length (sec)
\param[in]     do_capsnow                    [bool] true => do snow capping
\param[in]     oldfflag                      [int]  use old fsno parameterization
\param[in]     forc_t                        [double] atmospheric temperature (Kelvin)
\param[in]     t_grnd                        [double] ground temperature (Kelvin)
\param[in]     qflx_snow_grnd                [double] snow on ground after interception (mm H2O/s) [+]
\param[in]     qflx_snow_melt                [double] snow melt from previous time step
\param[in]     n_melt                        [double] SCA shape parameter [-]
\param[in,out] snow_depth                    [double] snow height (m)
\param[in,out] h2osno                        [double] snow water (mm H2O)
\param[in,out] int_snow                      [double] integrated snowfall [mm]
\param[out]    swe_old[nlevsno]              [double] snow water before update
\param[out]    h2osoi_liq[nlevgrnd+nlevsno]  [double] liquid water (kg/m2)
\param[out]    h2osoi_ice[nlevgrnd+nlevsno]  [double] ice lens (kg/m2)
\param[out]    t_soisno[nlevgrnd+nlevsno]    [double] soil temperature (Kelvin)
\param[out]    frac_iceold[nlevgrnd+nlevsno] [double] fraction of ice relative to the tot water
\param[out]    snl                           [int] number of snow layers
\param[out]    dz                            [double] layer thickness (m)
\param[out]    z                             [double] layer cell center elevation (m)
\param[out]    zi                            [double] layer interface elevation (m)
\param[out]    snw_rds[nlevsno]              [double] snow grain radius [m^-6, microns]
\param[out]    frac_sno_eff                  [double] fraction of ground covered by snow (0 to 1)
\param[out]    frac_sno                      [double] fraction of ground covered by snow (0 to 1)
*/
template <class ArrayD1>
ACCELERATE
void snow_init(const LandType& Land, const double& dtime, const bool& do_capsnow, const int& oldfflag,
               const double& forc_t, const double& t_grnd, const double& qflx_snow_grnd, const double& qflx_snow_melt,
               const double& n_melt,

               double& snow_depth, double& h2osno, double& int_snow, ArrayD1 swe_old, ArrayD1 h2osoi_liq,
               ArrayD1 h2osoi_ice, ArrayD1 t_soisno, ArrayD1 frac_iceold, int& snl, ArrayD1 dz, ArrayD1 z, ArrayD1 zi,
               ArrayD1 snw_rds, double& frac_sno_eff, double& frac_sno);

/*! Determine fraction of land surfaces which are submerged
based on surface microtopography and surface water storage.

\param[in]  Land                         [LandType] struct containing information about landtype
\param[in]  micro_sigma                  [double] microtopography pdf sigma (m)
\param[in]  h2osno                       [double] snow water (mm H2O)
\param[out] h2osfc                       [double] surface water (mm)
\param[out] h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
\param[out] frac_sno                     [double] fraction of ground covered by snow (0 to 1)
\param[out] frac_sno_eff                 [double] effective fraction of ground covered by snow (0 to 1)
\param[out] frac_h2osfc                  [double] fractional area with surface water greater than zero (0 to 1)
*/
template <class ArrayD1>
ACCELERATE
void fraction_h2osfc(const LandType& Land, const double& micro_sigma, const double& h2osno,

                     double& h2osfc, ArrayD1 h2osoi_liq, double& frac_sno, double& frac_sno_eff, double& frac_h2osfc);

} // namespace ELM::canopy_hydrology

#include "canopy_hydrology_impl.hh"
