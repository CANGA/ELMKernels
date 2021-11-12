
// h2ocan update - eq 7.9 in CLM 4.5 tech note - 
//h2ocan(n+1) = h2ocan(n) + intercepted_precip * dt - qflx_candrip * dt - evap * dt
// - evap * dt portion gets calculated in CanopyFluxes


#include "ELMConstants.h"
#include "LandType.h"
#include <algorithm>
#include <cmath>

namespace ELM {

// determine fraction of input precipitation that is snow and rain
double phase_precip_fraction(const double &total_precip, const double &phase_precip) {
  return total_precip > 0.0 ? phase_precip / total_precip : 0.0;
}

// set fraction of potential interception to max 0.25
double fraction_of_potential_interception(const double &elai, const double &esai) {
  return 0.25 * (1.0 - exp(-0.5 * (elai + esai)));
}

// Direct throughfall
double phase_throughfall(const double &phase_precip, const double &fpi) {
  return phase_precip > 0.0 ? phase_precip * (1.0 - fpi) : 0.0;
}


// maximum allowed water on canopy [mm]
double max_canopy_water(const double &dewmx, const double &elai, const double &esai) {
  return dewmx * (elai + esai);
}

// Intercepted precipitation [mm/s]
double intercepted_precip(const double total_precip, const double &fpi) {
  return total_precip * fpi;
}

// Excess water that exceeds the leaf capacity
double excess_canopy_water(const double &h2ocan, const double &h2ocanmx, const double &dtime) {
  return (h2ocan - h2ocanmx) / dtime;
}


// Canopy storage of interception and canopy drainage
void canopy_storage_and_drainage(const LandType &Land, const int &frac_veg_nosno, const double &total_precip, 
  const double &fpi, const double &dtime, const double &dewmx, const double &elai, const double &esai, double &h2ocan, 
  double &qflx_candrip) {
  if (!Land.lakpoi) {
    if (Land.ltype != istice && Land.ltype != istice_mec) {
      if (Land.ctype != icol_sunwall && Land.ctype != icol_shadewall) {
        if (frac_veg_nosno == 1 && total_precip > 0.0) {
          h2ocan += std::max(0.0, (dtime * intercepted_precip(total_precip, fpi) ));
          qflx_candrip = 0.0;
          double h2ocanmx = max_canopy_water(dewmx, elai, esai);
          double xrun = excess_canopy_water(h2ocan, h2ocanmx, dtime);
          if (xrun > 0.0) {
            qflx_candrip = xrun;
            h2ocan = h2ocanmx;
          }
        }
      }
    } else {
      h2ocan = 0.0;
      qflx_candrip = 0.0;
    }
  }
}

double phase_ground_precip(const LandType &Land, const int &frac_veg_nosno, const double &total_precip, 
  const double &phase_precip, const double &fpi, const double &qflx_candrip) {
  double phase_qflx_prec_grnd = 0.0;
  if ((Land.ctype != icol_sunwall) && (Land.ctype != icol_shadewall)) {
    if (frac_veg_nosno == 0) {
      phase_qflx_prec_grnd = phase_precip;
    } else {
      phase_qflx_prec_grnd = phase_throughfall(phase_precip, fpi) + (qflx_candrip * 
        phase_precip_fraction(total_precip, phase_precip));
    }
  } 
  return phase_qflx_prec_grnd;
}


void ground_precip_fluxes(const LandType &Land, const bool &do_capsnow, const int &frac_veg_nosno, 
  const double &total_precip, const double &forc_rain, const double &forc_snow, const double &fpi, 
  const double &qflx_candrip, const double &qflx_irrig,double &qflx_snwcp_liq, double &qflx_snwcp_ice, 
  double &qflx_snow_grnd, double &qflx_rain_grnd, double &qflx_prec_grnd) {
  if (!Land.lakpoi) {
    if (do_capsnow) {
      // Add irrigation water directly onto ground (bypassing canopy interception)
      qflx_snwcp_liq = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_rain, fpi, qflx_candrip) + qflx_irrig;
      qflx_snwcp_ice = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_snow, fpi, qflx_candrip);
      qflx_snow_grnd = 0.0;
      qflx_rain_grnd = 0.0;
      // Total water onto ground
      qflx_prec_grnd = qflx_snwcp_ice + qflx_snwcp_liq;
    } else {
      qflx_snwcp_liq = 0.0;
      qflx_snwcp_ice = 0.0;
      // Add irrigation water directly onto ground (bypassing canopy interception)
      qflx_rain_grnd = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_rain, fpi, qflx_candrip) + qflx_irrig;;
      qflx_snow_grnd = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_snow, fpi, qflx_candrip);
      // Total water onto ground
      qflx_prec_grnd = qflx_snow_grnd + qflx_rain_grnd;
    }
  }
}




void interception_physics(const LandType &Land, const bool &do_capsnow, const int &frac_veg_nosno, 
  const double &forc_rain, const double &forc_snow, const double &dtime, const double &dewmx, const double &elai, 
  const double &esai, double &h2ocan, const double &qflx_irrig, double &qflx_snwcp_liq, 
  double &qflx_snwcp_ice, double &qflx_snow_grnd, double &qflx_rain_grnd, double &qflx_prec_grnd) {

  double qflx_candrip;
  double total_precip = forc_snow + forc_rain;
  
  double fpi = fraction_of_potential_interception(elai, esai);
  
  canopy_storage_and_drainage(Land, frac_veg_nosno, total_precip, 
    fpi, dtime, dewmx, elai, esai, h2ocan, qflx_candrip);
  
  ground_precip_fluxes(Land, do_capsnow, frac_veg_nosno, total_precip, forc_rain, forc_snow, fpi, qflx_candrip, 
    qflx_irrig, qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd, qflx_rain_grnd, qflx_prec_grnd);

}














} // namespace ELM

