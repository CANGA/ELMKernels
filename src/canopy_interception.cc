/*
DESCRIPTION:
Calculate interception and partition incoming precipitation
into canopy storage, canopy runoff, rain and snow throughfall.
Also calculates fraction of precipitation that is rain/snow.

INPUTS:
ltype             [int]    landunit type
urbpoi            [bool]   true => landunit is an urban point
ctype             [int]    urban column type
frac_veg_nosno    [int]    fraction of veg not covered by snow (0/1 now) [-]
forc_rain         [double] rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
forc_snow         [double] snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
dewmx             [double] Maximum allowed dew [mm]
elai              [double] one-sided leaf area index with burying by snow
esai              [double] one-sided stem area index with burying by snow
dtime             [double] time step length (sec)

OUTPUTS:
h2ocan            [double] total canopy water (mm H2O) 
qflx_candrip      [double] rate of canopy runoff and snow falling off canopy [mm/s]
qflx_through_snow [double] direct rain throughfall [mm/s]
qflx_through_rain [double] direct snow throughfall [mm/s]
qflx_prec_intr    [double] interception of precipitation [mm/s] 
fracsnow          [double] frac of precipitation that is snow [-]
fracrain          [double] frac of precipitation that is rain [-]
*/

#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

void CanopyInterception(const int& ltype,
    const bool& urbpoi,
    const int& ctype,
    const int& frac_veg_nosno,
    const double& forc_rain,
    const double& forc_snow,
    const double& dewmx,
    const double& elai,
    const double& esai,
    const double& dtime,

    
    double& h2ocan,
    double& qflx_candrip,
    double& qflx_through_snow,
    double& qflx_through_rain,
    double& qflx_prec_intr,
    double& fracsnow,
    double& fracrain)


{
  //Canopy interception/storage and throughfall
  //Add precipitation to leaf water
  if (ltype == istsoil || ltype == istwet || urbpoi || ltype == istcrop) {
    qflx_candrip = 0.0;      // rate of canopy runoff
    qflx_through_snow = 0.0; // rain precipitation direct through canopy
    qflx_through_rain = 0.0; // snow precipitation direct through canopy
    qflx_prec_intr = 0.0;    // total intercepted precipitation
    fracsnow = 0.0;          // fraction of input precip that is snow
    fracrain = 0.0;          // fraction of input precip that is rain

    if (ctype != icol_sunwall && ctype != icol_shadewall) {

      if (frac_veg_nosno == 1 && (forc_rain + forc_snow) > 0.0) {

        // determine fraction of input precipitation that is snow and rain
        fracsnow = forc_snow/(forc_snow + forc_rain);
        fracrain = forc_rain/(forc_snow + forc_rain);

        // The leaf water capacities for solid and liquid are different,
        // generally double for snow, but these are of somewhat less
        // significance for the water budget because of lower evap. rate at
        // lower temperature.  Hence, it is reasonable to assume that
        // vegetation storage of solid water is the same as liquid water.
        double h2ocanmx = dewmx * (elai + esai);

        // Coefficient of interception
        // set fraction of potential interception to max 0.25
        double fpi = 0.25 * (1.0 - exp(-0.5*(elai + esai)));

        // Direct throughfall
        qflx_through_snow = forc_snow * (1.0-fpi);
        qflx_through_rain = forc_rain * (1.0-fpi);

        // Intercepted precipitation [mm/s]
        qflx_prec_intr = (forc_snow + forc_rain) * fpi;

        // Water storage of intercepted precipitation and dew
        h2ocan = std::max(0.0, (h2ocan + dtime*qflx_prec_intr));

        // Initialize rate of canopy runoff and snow falling off canopy
        qflx_candrip = 0.0;

        // Excess water that exceeds the leaf capacity
        double xrun = (h2ocan - h2ocanmx)/dtime;

        // Test on maximum dew on leaf
        // Note if xrun > 0 then h2ocan must be at least h2ocanmx
        if (xrun > 0.0) {
           qflx_candrip = xrun;
           h2ocan = h2ocanmx;
         }
       }
     }
   } else if (ltype == istice||ltype == istice_mec) {
    h2ocan            = 0.0;
    qflx_candrip      = 0.0;
    qflx_through_snow = 0.0;
    qflx_through_rain = 0.0;
    qflx_prec_intr    = 0.0;
    fracsnow          = 0.0;
    fracrain          = 0.0;
  }
}
