
#include "incident_shortwave.h"
#include "elm_constants.h"
#include <cmath>

namespace ns = ELM::incident_shortwave;

const double TWO_PI = ELM::ELM_PI * 2.0;
const double TWO_OVER_PI = 2.0 / ELM::ELM_PI;
const double PI_OVER_TWO = ELM::ELM_PI / 2.0;

// declination angle calc from ats/landlab
double ns::declination_angle(int doy) { return 23.45 * ELM_PI / 180.0 * cos(TWO_PI / 365.0 * (172.0 - doy)); }

// declination angle calc from ELM lnd_import szenith()/shr_orb_cosz()
double ns::declination_angle2(int doy) { return 23.45 * ELM_PI / 180.0 * sin(TWO_PI * (284.0 + doy) / 365.0); }

// cosine of the solar zenith angle
double ns::coszen(double latrad, double lonrad, double jday) {
  double decrad = declination_angle2(floor(jday));
  double cosz = sin(latrad) * sin(decrad) - cos(latrad) * cos(decrad) * cos((jday - floor(jday)) * TWO_PI + lonrad);
  return cosz > 0.001 ? cosz : 0.001;
}

// average coszen functions
// derived from shr_orb_avg_cosz() in shr_orb_mod.F90

// adjust variable (latitude or declination angle) so that its tangent will be defined
double ns::ensure_tan_defined(double var) {
  return (var == PI_OVER_TWO) ? var - 1.0e-05 : (var == -PI_OVER_TWO) ? var + 1.0e-05 : var;
}

// convert model dt from seconds to radians wrt daylength
double ns::dt_radians(double dt) { return dt * TWO_PI / 86400.0; }

// define dt start time of day on the period -pi to pi
double ns::dt_start_rad(double jday, double lonrad) {
  // adjust t to be between -2pi and 2pi
  double t_start = (jday - floor(jday)) * TWO_PI + lonrad - ELM_PI;
  return (t_start >= ELM_PI) ? t_start - TWO_PI : (t_start < -ELM_PI) ? t_start + TWO_PI : t_start;
}

// define time of day at end of dt
double ns::dt_end_rad(double t_start, double dtrad) { return t_start + dtrad; }

// define the cosine of the half-day length [0 to pi]
// adjust for cases of all daylight or all night
double ns::coshalfday(double latrad, double declin) {
  double cos_h = -tan(ensure_tan_defined(latrad)) * tan(ensure_tan_defined(declin));
  return (cos_h <= -1.0) ? ELM_PI : (cos_h >= 1.0) ? 0.0 : acos(cos_h);
}

// define the hour angle
// force it to be between -cos_h and cos_h
// consider the situation when the night period is too short
void ns::avg_hourangle(double t_start, double t_end, double dtrad, double cos_h, double *hour_angle) {
  if (t_end >= ELM_PI && t_start <= ELM_PI && ELM_PI - cos_h <= dtrad) {
    hour_angle[0] = std::min(std::max(t_start, -cos_h), cos_h);
    hour_angle[1] = cos_h;
    hour_angle[2] = TWO_PI - cos_h;
    hour_angle[3] = std::min(std::max(t_end, TWO_PI - cos_h), TWO_PI + cos_h);
  } else if (t_end >= -ELM_PI && t_start <= -ELM_PI && ELM_PI - cos_h <= dtrad) {
    hour_angle[0] = std::min(std::max(t_start, -TWO_PI - cos_h), -TWO_PI + cos_h);
    hour_angle[1] = -TWO_PI + cos_h;
    hour_angle[2] = -cos_h;
    hour_angle[3] = std::min(std::max(t_end, -cos_h), cos_h);
  } else {
    if (t_start > ELM_PI) {
      hour_angle[0] = std::min(std::max(t_start - TWO_PI, -cos_h), cos_h);
    } else if (t_start < -ELM_PI) {
      hour_angle[0] = std::min(std::max(t_start + TWO_PI, -cos_h), cos_h);
    } else {
      hour_angle[0] = std::min(std::max(t_start, -cos_h), cos_h);
    }
    if (t_end > ELM_PI) {
      hour_angle[1] = std::min(std::max(t_end - TWO_PI, -cos_h), cos_h);
    } else if (t_end < -ELM_PI) {
      hour_angle[1] = std::min(std::max(t_end + TWO_PI, -cos_h), cos_h);
    } else {
      hour_angle[1] = std::min(std::max(t_end, -cos_h), cos_h);
    }
    hour_angle[2] = 0.0;
    hour_angle[3] = 0.0;
  }
}

// perform a time integration to obtain cosz if desired
// output is valid over the period from t to t + dt
double ns::integrate_cosz(double t_start, double t_end, double dtrad, double cos_h, double latrad, double declin) {
  // define terms needed in the cosine zenith angle equation
  double aa = sin(latrad) * sin(declin);
  double bb = cos(latrad) * cos(declin);
  double ha[4];
  avg_hourangle(t_start, t_end, dtrad, cos_h, ha);
  return (ha[1] > ha[0] || ha[3] > ha[2]) ? (aa * (ha[1] - ha[0]) + bb * (sin(ha[1]) - sin(ha[0]))) / dtrad +
                                                (aa * (ha[3] - ha[2]) + bb * (sin(ha[3]) - sin(ha[2]))) / dtrad
                                          : 0.0;
}

// evaluate average cosine(zenith) for a given dt
double ns::average_cosz(double latrad, double lonrad, double declin, double dt, double jday) {
  double dtrad = dt_radians(dt);
  double t_start = dt_start_rad(jday, lonrad);
  double t_end = dt_end_rad(t_start, dtrad);
  double cos_h = coshalfday(latrad, declin);
  return integrate_cosz(t_start, t_end, dtrad, cos_h, latrad, declin);
}
