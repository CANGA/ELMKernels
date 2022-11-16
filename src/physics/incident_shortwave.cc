
#include "incident_shortwave.h"
#include "elm_constants.h"
#include <cmath>

namespace ns = ELM::incident_shortwave;
using ELM::ELMconst::ELM_PI;

namespace {
static constexpr double TWO_PI{ELM_PI * 2.0};
static constexpr double PI_OVER_TWO{ELM_PI / 2.0};
}

// declination angle calc from ats/landlab
double ns::declination_angle_cos(const int& doy) { return 23.45 * ELM_PI / 180.0 * cos(TWO_PI / 365.0 * (172.0 - doy)); }

// declination angle calc from ELM lnd_import szenith()/shr_orb_cosz()
double ns::declination_angle_sin(const int& doy) { return 23.45 * ELM_PI / 180.0 * sin(TWO_PI * (284.0 + doy) / 365.0); }

// cosine of the solar zenith angle
double ns::coszen(const double& latrad, const double& lonrad, const double& jday) {
  const double decrad{declination_angle_sin(floor(jday))};
  double cosz = sin(latrad) * sin(decrad) - cos(latrad) * cos(decrad) * cos((jday - floor(jday)) * TWO_PI + lonrad);
  return cosz > 0.001 ? cosz : 0.001;
}

// average coszen functions
// derived from shr_orb_avg_cosz() in shr_orb_mod.F90

// adjust variable (latitude or declination angle) so that its tangent will be defined
double ns::ensure_tan_defined(const double& var) {
  return (var == PI_OVER_TWO) ? var - 1.0e-05 : (var == -PI_OVER_TWO) ? var + 1.0e-05 : var;
}

// convert model dt from seconds to radians wrt daylength
double ns::dt_radians(const double& dt) { return dt * TWO_PI / 86400.0; }

// define dt start time of day on the period -pi to pi
double ns::dt_start_rad(const double& jday, const double& lonrad)
{
  // adjust t to be between -2pi and 2pi
  double t_start = (jday - floor(jday)) * TWO_PI + lonrad - ELM_PI;
  return (t_start >= ELM_PI) ? t_start - TWO_PI : (t_start < -ELM_PI) ? t_start + TWO_PI : t_start;
}

// define time of day at end of dt
double ns::dt_end_rad(const double& t_start, const double& dtrad)
{ 
  return t_start + dtrad;
}

// define the cosine of the half-day length [0 to pi]
// adjust for cases of all daylight or all night
double ns::coshalfday(const double& latrad, const double& declin)
{
  double cos_h = -tan(ensure_tan_defined(latrad)) * tan(ensure_tan_defined(declin));
  return (cos_h <= -1.0) ? ELM_PI : (cos_h >= 1.0) ? 0.0 : acos(cos_h);
}

// define the hour angle
// force it to be between -cos_h and cos_h
// consider the situation when the night period is too short
void ns::avg_hourangle(const double& t_start, const double& t_end, const double& dtrad,
                       const double& cos_h, double hour_angle[4])
{
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
double ns::integrate_cosz(const double& t_start, const double& t_end, const double& dtrad,
                          const double& cos_h, const double& latrad, const double& declin)
{
  // define terms needed in the cosine zenith angle equation
  const double aa{sin(latrad) * sin(declin)};
  const double bb{cos(latrad) * cos(declin)};
  double ha[4];
  avg_hourangle(t_start, t_end, dtrad, cos_h, ha);
  return (ha[1] > ha[0] || ha[3] > ha[2]) ? 
    (aa * (ha[1] - ha[0]) + bb * (sin(ha[1]) - sin(ha[0]))) / dtrad +
    (aa * (ha[3] - ha[2]) + bb * (sin(ha[3]) - sin(ha[2]))) / dtrad
    : 0.0;
}

// evaluate average cosine(zenith) for a given dt
double ns::average_cosz(const double& latrad, const double&  lonrad, const double&  dt, const double& jday)
{
  const double dtrad{dt_radians(dt)};
  const double t_start{dt_start_rad(jday, lonrad)};
  const double t_end{dt_end_rad(t_start, dtrad)};
  const double declin{declination_angle_sin(static_cast<int>(jday))};
  const double cos_h{coshalfday(latrad, declin)};
  return integrate_cosz(t_start, t_end, dtrad, cos_h, latrad, declin);
}
