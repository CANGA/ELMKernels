
#pragma once

#include "clm_constants.h"
#include "landtype.h"
#include "vegproperties.h"

namespace ELM {

/* InitializeFlux_Can()
INPUTS:
Land                            [LandType] struct containing information about landtype
snl                             [int] number of snow layers
frac_veg_nosno                  [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
frac_sno                        [double] fraction of ground covered by snow (0 to 1)
forc_hgt_u_patch                [double] observational height of wind at pft level [m]
thm                             [double] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
thv                             [double] virtual potential temperature (kelvin)
max_dayl                        [double] maximum daylength for this grid cell (s)
dayl                            [double] daylength (seconds)
altmax_indx                     [int] index corresponding to maximum active layer depth from current year
altmax_lastyear_indx            [int] index corresponding to maximum active layer depth from prior year
t_soisno[nlevgrnd+nlevsno]      [double] col soil temperature (Kelvin)
h2osoi_ice[nlevgrnd+nlevsno]    [double] ice lens (kg/m2)
h2osoi_liq[nlevgrnd+nlevsno]    [double] liquid water (kg/m2)
dz[nlevgrnd+nlevsno]            [double] layer thickness (m)
rootfr[nlevgrnd]                [double] fraction of roots in each soil layer
tc_stress                       [double] critical soil temperature for soil water stress (C)
sucsat[nlevgrnd]                [double] minimum soil suction (mm)
watsat[nlevgrnd]                [double] volumetric soil water at saturation (porosity)
bsw[nlevgrnd]                   [double] Clapp and Hornberger "b
smpso[numpft]                   [double] soil water potential at full stomatal opening (mm)
smpsc[numpft]                   [double] soil water potential at full stomatal closure (mm)
elai                            [double] one-sided leaf area index with burying by snow
esai                            [double] one-sided stem area index with burying by snow
emv                             [double] vegetation emissivity
emg                             [double] ground emissivity
qg                              [double] ground specific humidity [kg/kg]
t_grnd                          [double] ground temperature (Kelvin)
forc_t                          [double] atmospheric temperature (Kelvin)
forc_pbot                       [double] atmospheric pressure (Pa)
forc_lwrad                      [double] downward infrared (longwave) radiation (W/m**2)
forc_u                          [double] atmospheric wind speed in east direction (m/s)
forc_v                          [double] atmospheric wind speed in north direction (m/s)
forc_q                          [double] atmospheric specific humidity (kg/kg)
forc_th                         [double] atmospheric potential temperature (Kelvin)
z0mg                            [double] roughness length over ground, momentum [m]

OUTPUTS:
btran                           [double]  transpiration wetness factor (0 to 1)
z0mv                            [double] roughness length over vegetation, momentum [m]
z0hv                            [double] roughness length over vegetation, sensible heat [m]
z0qv                            [double] roughness length over vegetation, latent heat [m]
displa                          [double]  displacement height (m)
rootr[nlevgrnd]                 [double] effective fraction of roots in each soil layer
eff_porosity[nlevgrnd]          [double] effective soil porosity
h2osoi_liqvol[nlevgrnd+nlevsno] [double] liquid volumetric moisture
dayl_factor                     [double] scalar (0-1) for daylength effect on Vcmax
air                             [double] atmos. radiation temporay set
bir                             [double] atmos. radiation temporay set
cir                             [double] atmos. radiation temporay set
el                              [double] vapor pressure on leaf surface [pa]
qsatl                           [double] leaf specific humidity [kg/kg]
qsatldT                         [double] derivative of "qsatl" on "t_veg"
taf                             [double] air temperature within canopy space [K]
qaf                             [double] humidity of canopy air [kg/kg]
um                              [double] wind speed including the stablity effect [m/s]
ur                              [double] wind speed at reference height [m/s]
obu                             [double] Monin-Obukhov length (m)
zldis                           [double] reference height "minus" zero displacement height [m]
delq                            [double] temporary
t_veg                           [double]  vegetation temperature (Kelvin)
*/
template <class dArray_type>
void InitializeFlux_Can(const LandType &Land, const int &snl, const int &frac_veg_nosno, const double &frac_sno,
                        const double &forc_hgt_u_patch, const double &thm, const double &thv, const double &max_dayl,
                        const double &dayl, const int &altmax_indx, const int &altmax_lastyear_indx,
                        const dArray_type t_soisno, const dArray_type h2osoi_ice, const dArray_type h2osoi_liq,
                        const dArray_type dz, const dArray_type rootfr, const double &tc_stress,
                        const dArray_type sucsat, const dArray_type watsat, const dArray_type bsw,
                        const dArray_type smpso, const dArray_type smpsc, const double &elai, const double &esai,
                        const double &emv, const double &emg, const double &qg, const double &t_grnd,
                        const double &forc_t, const double &forc_pbot, const double &forc_lwrad, const double &forc_u,
                        const double &forc_v, const double &forc_q, const double &forc_th, const double &z0mg,
                        double &btran, double &displa, double &z0mv, double &z0hv, double &z0qv, dArray_type rootr,
                        dArray_type eff_porosity, dArray_type h2osoi_liqvol, double &dayl_factor, double &air,
                        double &bir, double &cir, double &el, double &qsatl, double &qsatldT, double &taf, double &qaf,
                        double &um, double &ur, double &obu, double &zldis, double &delq, double &t_veg);

/* StabilityIteration_Can()
DESCRIPTION:
calculate Monin-Obukhov length and wind speed, call photosynthesis, calculate ET & SH flux
Iterates until convergence, up to 40 iterations, calling friction velocity functions, then
photosynthesis for sun & shade.

INPUTS:
Land                       [LandType] struct containing information about landtype
dtime                      [double] timestep size (sec)
snl                        [int] number of snow layers
frac_veg_nosno             [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
frac_sno                   [double] fraction of ground covered by snow (0 to 1)
forc_hgt_u_patch           [double] observational height of wind at pft level [m]
forc_hgt_t_patch           [double] observational height of temperature at pft level [m]
forc_hgt_q_patch           [double] observational height of specific humidity at pft level [m]
dleaf[numpft]              [double] characteristic leaf dimension (m)
fwet                       [double] fraction of canopy that is wet (0 to 1)
fdry                       [double] fraction of foliage that is green and dry [-]
laisun                     [double] sunlit leaf area
laisha                     [double] shaded leaf area
forc_rho                   [double] air density (kg/m**3)
snow_depth                 [double] snow height (m)
soilbeta                   [double] soil wetness relative to field capacity
frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
t_h2osfc                   [double] surface water temperature
sabv                       [double] solar radiation absorbed by vegetation (W/m**2)
h2ocan                     [double] canopy water (mm H2O)
htop                       [double] canopy top(m)
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
air                        [double] atmos. radiation temporay set
bir                        [double] atmos. radiation temporay set
cir                        [double] atmos. radiation temporay set
ur                         [double] wind speed at reference height [m/s]
zldis                      [double] reference height "minus" zero displacement height [m]
displa                     [double]  displacement height (m)
elai                       [double] one-sided leaf area index with burying by snow
esai                       [double] one-sided stem area index with burying by snow
t_grnd                     [double] ground temperature (Kelvin)
forc_pbot                  [double]  atmospheric pressure (Pa)
forc_q                     [double] atmospheric specific humidity (kg/kg)
forc_th                    [double] atmospheric potential temperature (Kelvin)
z0mg                       [double] roughness length over ground, momentum [m]
z0mv                       [double] roughness length over vegetation, momentum [m]
z0hv                       [double] roughness length over vegetation, sensible heat [m]
z0qv                       [double] roughness length over vegetation, latent heat [m]
thm                        [double]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
thv                        [double] virtual potential temperature (kelvin)
qg                         [double] ground specific humidity [kg/kg]
veg                        [VegProperties] struct containing vegetation constant parameters
nrad                       [int]  number of canopy layers above snow for radiative transfer
t10                        [double] 10-day running mean of the 2 m temperature (K)
tlai_z[nlevcan]            [double] pft total leaf area index for canopy layer
vcmaxcintsha               [double] leaf to canopy scaling coefficient - shade
vcmaxcintsun               [double] leaf to canopy scaling coefficient - sun
parsha_z[nlevcan]          [double] par absorbed per unit lai for canopy layer (w/m**2) - shade
parsun_z[nlevcan]          [double] par absorbed per unit lai for canopy layer (w/m**2) - sun
laisha_z[nlevcan]          [double] leaf area index for canopy layer, sunlit or shaded - shade
laisun_z[nlevcan]          [double] leaf area index for canopy layer, sunlit or shaded - sun
forc_pco2                  [double]  partial pressure co2 (Pa)
forc_po2                   [double]  partial pressure o2 (Pa))
dayl_factor                [double] scalar (0-1) for daylength effect on Vcmax

OUTPUTS:
btran                      [double]  transpiration wetness factor (0 to 1)
qflx_tran_veg              [double] vegetation transpiration (mm H2O/s) (+ = to atm)
qflx_evap_veg              [double] vegetation evaporation (mm H2O/s) (+ = to atm)
eflx_sh_veg                [double] sensible heat flux from leaves (W/m**2) [+ to atm]
wtg                        [double] heat conductance for ground [m/s]
wtl0                       [double] normalized heat conductance for leaf [-]
wta0                       [double] normalized heat conductance for air [-]
wtal                       [double] normalized heat conductance for air and leaf [-]
el                         [double] vapor pressure on leaf surface [pa]
qsatl                      [double] leaf specific humidity [kg/kg]
qsatldT                    [double] derivative of "qsatl" on "t_veg"
taf                        [double] air temperature within canopy space [K]
qaf                        [double] humidity of canopy air [kg/kg]
um                         [double] wind speed including the stablity effect [m/s]
dth                        [double] diff of virtual temp. between ref. height and surface
dqh                        [double] diff of humidity between ref. height and surface
obu                        [double] Monin-Obukhov length (m)
temp1                      [double] relation for potential temperature profile
temp2                      [double] relation for specific humidity profile
temp12m                    [double] relation for potential temperature profile applied at 2-m
temp22m                    [double] relation for specific humidity profile applied at 2-m
tlbef                      [double] leaf temperature from previous iteration [K]
delq                       [double] temporary
dt_veg                     [double] change in t_veg, last iteration (Kelvin)
t_veg                      [double]  vegetation temperature (Kelvin)
wtgq                       [double] latent heat conductance for ground [m/s]
wtalq                      [double] normalized latent heat cond. for air and leaf [-]
wtlq0                      [double] normalized latent heat conductance for leaf [-]
wtaq0                      [double] normalized latent heat conductance for air [-]
*/
template <class dArray_type>
void StabilityIteration_Can(const LandType &Land, const double &dtime, const int &snl, const int &frac_veg_nosno,
                            const double &frac_sno, const double &forc_hgt_u_patch, const double &forc_hgt_t_patch,
                            const double &forc_hgt_q_patch, const dArray_type dleaf, const double &fwet,
                            const double &fdry, const double &laisun, const double &laisha, const double &forc_rho,
                            const double &snow_depth, const double &soilbeta, const double &frac_h2osfc,
                            const double &t_h2osfc, const double &sabv, const double &h2ocan, const double &htop,
                            const dArray_type t_soisno, const double &air, const double &bir, const double &cir,
                            const double &ur, const double &zldis, const double &displa, const double &elai,
                            const double &esai, const double &t_grnd, const double &forc_pbot, const double &forc_q,
                            const double &forc_th, const double &z0mg, const double &z0mv, const double &z0hv,
                            const double &z0qv, const double &thm, const double &thv, const double &qg,
                            const VegProperties &veg, const int &nrad, const double &t10, const dArray_type tlai_z,
                            const double &vcmaxcintsha, const double &vcmaxcintsun, const dArray_type parsha_z,
                            const dArray_type parsun_z, const dArray_type laisha_z, const dArray_type laisun_z,
                            const double &forc_pco2, const double &forc_po2, const double &dayl_factor, double &btran,
                            double &qflx_tran_veg, double &qflx_evap_veg, double &eflx_sh_veg, double &wtg,
                            double &wtl0, double &wta0, double &wtal, double &el, double &qsatl, double &qsatldT,
                            double &taf, double &qaf, double &um, double &dth, double &dqh, double &obu, double &temp1,
                            double &temp2, double &temp12m, double &temp22m, double &tlbef, double &delq,
                            double &dt_veg, double &t_veg, double &wtgq, double &wtalq, double &wtlq0, double &wtaq0);

/* ComputeFlux_Can()
INPUTS:
Land                       [LandType] struct containing information about landtype
dtime                      [double] timestep size (sec)
snl                        [int] number of snow layers
frac_veg_nosno             [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
frac_sno                   [double] fraction of ground covered by snow (0 to 1)
t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
t_h2osfc                   [double] surface water temperature
sabv                       [double] solar radiation absorbed by vegetation (W/m**2)
qg_snow                    [double] specific humidity at snow surface [kg/kg]
qg_soil                    [double] specific humidity at soil surface [kg/kg]
qg_h2osfc                  [double] specific humidity at h2osfc [kg/kg]
dqgdT                      [double] d(qg)/dT
htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]
wtg                        [double] heat conductance for ground [m/s]
wtl0                       [double]  normalized heat conductance for leaf [-]
wta0                       [double] normalized heat conductance for air [-]
wtal                       [double] normalized heat conductance for air and leaf [-]
air                        [double] atmos. radiation temporay set
bir                        [double] atmos. radiation temporay set
cir                        [double] atmos. radiation temporay set
qsatl                      [double] leaf specific humidity [kg/kg]
qsatldT                    [double] derivative of "qsatl" on "t_ve
dth                        [double] diff of virtual temp. between ref. height and surface
dqh                        [double] diff of humidity between ref. height and surface
temp1                      [double] relation for potential temperature profile
temp2                      [double] relation for specific humidity profile
temp12m                    [double] relation for potential temperature profile applied at 2-m
temp22m                    [double] relation for specific humidity profile applied at 2-m
tlbef                      [double] leaf temperature from previous iteration [K]
delq                       [double] temporary
dt_veg                     [double] change in t_veg, last iteration (Kelvin)
t_veg                      [double]  vegetation temperature (Kelvin)
t_grnd                     [double] ground temperature (Kelvin)
forc_pbot                  [double]  atmospheric pressure (Pa)
qflx_tran_veg              [double] vegetation transpiration (mm H2O/s) (+ = to atm)
qflx_evap_veg              [double] vegetation evaporation (mm H2O/s) (+ = to atm)
eflx_sh_veg                [double] sensible heat flux from leaves (W/m**2) [+ to atm]
forc_q                     [double] atmospheric specific humidity (kg/kg)
forc_rho                   [double] air density (kg/m**3)
thm                        [double]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
emv                        [double] vegetation emissivity
emg                        [double] ground emissivity
forc_lwrad                 [double]  downward infrared (longwave) radiation (W/m**2)
wtgq                       [double] latent heat conductance for ground [m/s]
wtalq                      [double] normalized latent heat cond. for air and leaf [-]
wtlq0                      [double] normalized latent heat conductance for leaf [-]
wtaq0                      [double] normalized latent heat conductance for air [-]

OUTPUTS:
h2ocan                     [double] canopy water (mm H2O)
eflx_sh_grnd               [double] sensible heat flux from ground (W/m**2) [+ to atm]
eflx_sh_snow               [double] sensible heat flux from snow (W/m**2) [+ to atm]
eflx_sh_soil               [double] sensible heat flux from soil (W/m**2) [+ to atm]
eflx_sh_h2osfc             [double] sensible heat flux from h2osfc (W/m**2) [+ to atm]
qflx_evap_soi              [double] soil evaporation (mm H2O/s) (+ = to atm)
qflx_ev_snow               [double] evaporation flux from snow (W/m**2) [+ to atm]
qflx_ev_soil               [double] evaporation flux from soil (W/m**2) [+ to atm]
qflx_ev_h2osfc             [double] evaporation flux from h2osfc (W/m**2) [+ to atm]
dlrad                      [double] downward longwave radiation below the canopy [W/m2]
ulrad                      [double] upward longwave radiation above the canopy [W/m2]
cgrnds                     [double] deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
cgrndl                     [double] deriv of soil latent heat flux wrt soil temp [w/m**2/k]
cgrnd                      [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
t_ref2m                    [double]  2 m height surface air temperature (Kelvin)
t_ref2m_r                  [double]  Rural 2 m height surface air temperature (Kelvin)
q_ref2m                    [double]  2 m height surface specific humidity (kg/kg)
rh_ref2m_r                 [double]  Rural 2 m height surface relative humidity (%)
rh_ref2m                   [double]  2 m height surface relative humidity (%)
*/
template <class dArray_type>
void ComputeFlux_Can(const LandType &Land, const double &dtime, const int &snl, const int &frac_veg_nosno,
                     const double &frac_sno, const dArray_type t_soisno, const double &frac_h2osfc,
                     const double &t_h2osfc, const double &sabv, const double &qg_snow, const double &qg_soil,
                     const double &qg_h2osfc, const double &dqgdT, const double &htvp, const double &wtg,
                     const double &wtl0, const double &wta0, const double &wtal, const double &air, const double &bir,
                     const double &cir, const double &qsatl, const double &qsatldT, const double &dth,
                     const double &dqh, const double &temp1, const double &temp2, const double &temp12m,
                     const double &temp22m, const double &tlbef, const double &delq, const double &dt_veg,
                     const double &t_veg, const double &t_grnd, const double &forc_pbot, const double &qflx_tran_veg,
                     const double &qflx_evap_veg, const double &eflx_sh_veg, const double &forc_q,
                     const double &forc_rho, const double &thm, const double &emv, const double &emg,
                     const double &forc_lwrad, const double &wtgq, const double &wtalq, const double &wtlq0,
                     const double &wtaq0, double &h2ocan, double &eflx_sh_grnd, double &eflx_sh_snow,
                     double &eflx_sh_soil, double &eflx_sh_h2osfc, double &qflx_evap_soi, double &qflx_ev_snow,
                     double &qflx_ev_soil, double &qflx_ev_h2osfc, double &dlrad, double &ulrad, double &cgrnds,
                     double &cgrndl, double &cgrnd, double &t_ref2m, double &t_ref2m_r, double &q_ref2m,
                     double &rh_ref2m, double &rh_ref2m_r);

} // namespace ELM

#include "CanopyFluxes_impl.hh"