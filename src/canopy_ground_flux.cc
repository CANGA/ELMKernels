/*
DESCRIPTION:
Add liquid and solid water inputs to ground surface after interception and 
canopy storage losses.

INPUTS:
ctype             [int]    urban column type
do_capsnow        [bool]   true => do snow capping
frac_veg_nosno    [int]    fraction of veg not covered by snow (0/1 now) [-]
forc_rain         [double] rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
forc_snow         [double] snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
qflx_candrip      [double] rate of canopy runoff and snow falling off canopy [mm/s]
qflx_through_snow [double] direct rain throughfall [mm/s]
qflx_through_rain [double] direct snow throughfall [mm/s]
fracsnow          [double] frac of precipitation that is snow [-]
fracrain          [double] frac of precipitation that is rain [-]
qflx_irrig        [double] irrigation amount (mm/s)

OUTPUTS:
qflx_prec_grnd    [double] water onto ground including canopy runoff [kg/(m2 s)]
qflx_snwcp_liq    [double] excess rainfall due to snow capping (mm H2O /s) [+]
qflx_snwcp_ice    [double] excess snowfall due to snow capping (mm H2O /s) [+]
qflx_snow_grnd    [double] snow on ground after interception (mm H2O/s) [+]
qflx_rain_grnd    [double] rain on ground after interception (mm H2O/s) [+]
*/

#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

void CanopyGroundFlux(const int& ctype,
  const bool& do_capsnow,
  const int& frac_veg_nosno,
  const double& forc_rain,
  const double& forc_snow,
  const double& qflx_candrip,
  const double& qflx_through_snow,
  const double& qflx_through_rain,
  const double& fracsnow,
  const double& fracrain,
  const double& qflx_irrig,


  double& qflx_prec_grnd,
  double& qflx_snwcp_liq,
  double& qflx_snwcp_ice,
  double& qflx_snow_grnd,
  double& qflx_rain_grnd)
{

  double qflx_prec_grnd_snow, qflx_prec_grnd_rain;
  double qflx_dirct_rain, qflx_leafdrip; // these are not currently needed, but should be passed out if SedYield is desired 

// Precipitation onto ground (kg/(m2 s))
  if ((ctype != icol_sunwall) && (ctype != icol_shadewall)) {
    if (frac_veg_nosno == 0) {
      qflx_prec_grnd_snow = forc_snow;
      qflx_prec_grnd_rain = forc_rain;
      qflx_dirct_rain = forc_rain;
      qflx_leafdrip = 0.0;
    } else {
      qflx_prec_grnd_snow = qflx_through_snow + (qflx_candrip * fracsnow);
      qflx_prec_grnd_rain = qflx_through_rain + (qflx_candrip * fracrain);
      qflx_dirct_rain = qflx_through_rain;
      qflx_leafdrip = qflx_candrip * fracrain;
    }
  } else {
  // Urban sunwall and shadewall have no intercepted precipitation
    qflx_prec_grnd_snow = 0.0;
    qflx_prec_grnd_rain = 0.0;
    qflx_dirct_rain = 0.0;
    qflx_leafdrip = 0.0;
  }

  // Add irrigation water directly onto ground (bypassing canopy interception)
  qflx_prec_grnd_rain = qflx_prec_grnd_rain + qflx_irrig;
  // Total water onto ground
  qflx_prec_grnd = qflx_prec_grnd_snow + qflx_prec_grnd_rain;

  if (do_capsnow) {
    qflx_snwcp_liq = qflx_prec_grnd_rain;
    qflx_snwcp_ice = qflx_prec_grnd_snow;

    qflx_snow_grnd = 0.0;
    qflx_rain_grnd = 0.0;
  } else {
    qflx_snwcp_liq = 0.0;
    qflx_snwcp_ice = 0.0;
    qflx_snow_grnd = qflx_prec_grnd_snow;          // ice onto ground (mm/s)
    qflx_rain_grnd = qflx_prec_grnd_rain;          // liquid water onto ground (mm/s)
  }
}