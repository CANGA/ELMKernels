/*
DESCRIPTION:
functions derived from CanopyHydrologyMod.F90

*/

#include "CanopyHydrology.h"
#include "clm_constants.h"
#include "landtype.h"
#include <algorithm>
#include <cmath>

namespace ELM {

/* CanopyHydrology::Interception()
DESCRIPTION:
Calculate interception and partition incoming precipitation
into canopy storage, canopy runoff, rain and snow throughfall.
Also calculates fraction of precipitation that is rain/snow.

INPUTS:
Land              [LandType] struct containing information about landtype
frac_veg_nosno    [int]      fraction of veg not covered by snow (0/1 now) [-]
forc_rain         [double]   rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
forc_snow         [double]   snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
dewmx             [double]   Maximum allowed dew [mm]
elai              [double]   one-sided leaf area index with burying by snow
esai              [double]   one-sided stem area index with burying by snow
dtime             [double]   time step length (sec)

OUTPUTS:
h2ocan            [double]   total canopy water (mm H2O)
qflx_candrip      [double]   rate of canopy runoff and snow falling off canopy [mm/s]
qflx_through_snow [double]   direct snow throughfall [mm/s]
qflx_through_rain [double]   direct rain throughfall [mm/s]
fracsnow          [double]   frac of precipitation that is snow [-]
fracrain          [double]   frac of precipitation that is rain [-]
*/
void Interception(const LandType &Land, const int &frac_veg_nosno, const double &forc_rain, const double &forc_snow,
                  const double &dewmx, const double &elai, const double &esai, const double &dtime, double &h2ocan,
                  double &qflx_candrip, double &qflx_through_snow, double &qflx_through_rain, double &fracsnow,
                  double &fracrain) {

  if (!Land.lakpoi) {
    // Canopy interception/storage and throughfall
    // Add precipitation to leaf water
    if (Land.ltype == istsoil || Land.ltype == istwet || Land.urbpoi || Land.ltype == istcrop) {
      qflx_candrip = 0.0;          // rate of canopy runoff
      qflx_through_snow = 0.0;     // rain precipitation direct through canopy
      qflx_through_rain = 0.0;     // snow precipitation direct through canopy
      double qflx_prec_intr = 0.0; // total intercepted precipitation
      fracsnow = 0.0;              // fraction of input precip that is snow
      fracrain = 0.0;              // fraction of input precip that is rain

      if (Land.ctype != icol_sunwall && Land.ctype != icol_shadewall) {
        if (frac_veg_nosno == 1 && (forc_rain + forc_snow) > 0.0) {
          // determine fraction of input precipitation that is snow and rain
          fracsnow = forc_snow / (forc_snow + forc_rain);
          fracrain = forc_rain / (forc_snow + forc_rain);
          // The leaf water capacities for solid and liquid are different,
          // generally double for snow, but these are of somewhat less
          // significance for the water budget because of lower evap. rate at
          // lower temperature.  Hence, it is reasonable to assume that
          // vegetation storage of solid water is the same as liquid water.
          double h2ocanmx = dewmx * (elai + esai);
          // Coefficient of interception
          // set fraction of potential interception to max 0.25
          double fpi = 0.25 * (1.0 - exp(-0.5 * (elai + esai)));
          // Direct throughfall
          qflx_through_snow = forc_snow * (1.0 - fpi);
          qflx_through_rain = forc_rain * (1.0 - fpi);
          // Intercepted precipitation [mm/s]
          qflx_prec_intr = (forc_snow + forc_rain) * fpi;
          // Water storage of intercepted precipitation and dew
          h2ocan = std::max(0.0, (h2ocan + dtime * qflx_prec_intr));
          // Initialize rate of canopy runoff and snow falling off canopy
          qflx_candrip = 0.0;
          // Excess water that exceeds the leaf capacity
          double xrun = (h2ocan - h2ocanmx) / dtime;
          // Test on maximum dew on leaf
          // Note if xrun > 0 then h2ocan must be at least h2ocanmx
          if (xrun > 0.0) {
            qflx_candrip = xrun;
            h2ocan = h2ocanmx;
          }
        }
      }
    } else if (Land.ltype == istice || Land.ltype == istice_mec) {
      h2ocan = 0.0;
      qflx_candrip = 0.0;
      qflx_through_snow = 0.0;
      qflx_through_rain = 0.0;
      double qflx_prec_intr = 0.0;
      fracsnow = 0.0;
      fracrain = 0.0;
    }
  }
} // Interception

/* CanopyHydrology::Irrigation()
DESCRIPTION:
Determine whether we're irrigating here; set qflx_irrig appropriately.
note - will likely need to figure out different method for updating n_irrig_steps_left

INPUTS:
Land               [LandType] struct containing information about landtype
irrig_rate         [double]   current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]

OUTPUTS:
n_irrig_steps_left [int] number of time steps for which we still need to irrigate today
qflx_irrig         [double]   irrigation amount (mm/s)
*/
void Irrigation(const LandType &Land, const double &irrig_rate, int &n_irrig_steps_left, double &qflx_irrig) {
  if (!Land.lakpoi) {
    if (n_irrig_steps_left > 0) {
      qflx_irrig = irrig_rate;
      n_irrig_steps_left -= 1;
    } else {
      qflx_irrig = 0.0;
    }
  }
} // Irrigation

/* CanopyHydrology::GroundFlux()
DESCRIPTION:
Add liquid and solid water inputs to ground surface after interception and
canopy storage losses.

INPUTS:
Land              [LandType] struct containing information about landtype
do_capsnow        [bool]     true => do snow capping
frac_veg_nosno    [int]      fraction of veg not covered by snow (0/1 now) [-]
forc_rain         [double]   rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
forc_snow         [double]   snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
qflx_irrig        [double]   irrigation amount (mm/s)
qflx_candrip      [double]   rate of canopy runoff and snow falling off canopy [mm/s]
qflx_through_snow [double]   direct snow throughfall [mm/s]
qflx_through_rain [double]   direct rain throughfall [mm/s]
fracsnow          [double]   frac of precipitation that is snow [-]
fracrain          [double]   frac of precipitation that is rain [-]

OUTPUTS:
qflx_prec_grnd    [double]   water onto ground including canopy runoff [kg/(m2 s)]
qflx_snwcp_liq    [double]   excess rainfall due to snow capping (mm H2O /s) [+]
qflx_snwcp_ice    [double]   excess snowfall due to snow capping (mm H2O /s) [+]
qflx_snow_grnd    [double]   snow on ground after interception (mm H2O/s) [+]
qflx_rain_grnd    [double]   rain on ground after interception (mm H2O/s) [+]
*/
void GroundFlux(const LandType &Land, const bool &do_capsnow, const int &frac_veg_nosno, const double &forc_rain,
                const double &forc_snow, const double &qflx_irrig, const double &qflx_candrip,
                const double &qflx_through_snow, const double &qflx_through_rain, const double &fracsnow,
                const double &fracrain, double &qflx_prec_grnd, double &qflx_snwcp_liq, double &qflx_snwcp_ice,
                double &qflx_snow_grnd, double &qflx_rain_grnd) {

  if (!Land.lakpoi) {
    double qflx_prec_grnd_snow, qflx_prec_grnd_rain;
    // Precipitation onto ground (kg/(m2 s))
    if ((Land.ctype != icol_sunwall) && (Land.ctype != icol_shadewall)) {
      if (frac_veg_nosno == 0) {
        qflx_prec_grnd_snow = forc_snow;
        qflx_prec_grnd_rain = forc_rain;
      } else {
        qflx_prec_grnd_snow = qflx_through_snow + (qflx_candrip * fracsnow);
        qflx_prec_grnd_rain = qflx_through_rain + (qflx_candrip * fracrain);
      }
    } else {
      // Urban sunwall and shadewall have no intercepted precipitation
      qflx_prec_grnd_snow = 0.0;
      qflx_prec_grnd_rain = 0.0;
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
      qflx_snow_grnd = qflx_prec_grnd_snow; // ice onto ground (mm/s)
      qflx_rain_grnd = qflx_prec_grnd_rain; // liquid water onto ground (mm/s)
    }
  }
} // GroundFlux

/* CanopyHydrology::FracWet()
DESCRIPTION:
Determine fraction of vegetated surfaces which are wet and
fraction of elai which is dry. The variable ``fwet'' is the
fraction of all vegetation surfaces which are wet including
stem area which contribute to evaporation. The variable ``fdry''
is the fraction of elai which is dry because only leaves
can transpire.  Adjusted for stem area which does not transpire.

INPUTS:
Land           [LandType] struct containing information about landtype
frac_veg_nosno [int]      fraction of veg not covered by snow (0/1 now) [-]
dewmx          [double]   Maximum allowed dew [mm]
elai           [double]   one-sided leaf area index with burying by snow
esai           [double]   one-sided stem area index with burying by snow
h2ocan         [double]   total canopy water (mm H2O)

OUTPUTS:
fwet           [double]   fraction of canopy that is wet (0 to 1)
fdry           [double]   fraction of foliage that is green and dry [-] (new)
*/
void FracWet(const LandType &Land, const int &frac_veg_nosno, const double &dewmx, const double &elai,
             const double &esai, const double &h2ocan, double &fwet, double &fdry) {

  if (!Land.lakpoi) {
    if (frac_veg_nosno == 1) {
      if (h2ocan > 0.0) {
        double vegt = frac_veg_nosno * (elai + esai);
        double dewmxi = 1.0 / dewmx;
        fwet = pow(((dewmxi / vegt) * h2ocan), 2.0 / 3.0);
        fwet = std::min(fwet, 1.0); // Check for maximum limit of fwet
      } else {
        fwet = 0.0;
      }
      fdry = (1.0 - fwet) * elai / (elai + esai);
    } else {
      fwet = 0.0;
      fdry = 0.0;
    }
  }
} // FracWet

} // namespace ELM
