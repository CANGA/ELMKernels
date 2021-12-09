/*! \file daylength.h
\brief Functions derived from DaylengthMod.F90
*/

#include "daylength.h"
#include "elm_constants.h"
#include <cmath>
#include <limits>
#include <assert.h>

namespace ELM {

// Computes daylength (in seconds)
// Latitude and solar declination angle should both be specified in radians. decl must
// be strictly less than pi/2; lat must be less than pi/2 within a small tolerance.
double daylength(const double& lat, const double& decl) {
  // number of seconds per radian of hour-angle
  const double secs_per_radian = 13750.9871;
  // epsilon for defining latitudes "near" the pole
  const double lat_epsilon = 10.0 * std::numeric_limits<double>::epsilon();
  // Define an offset pole as slightly less than pi/2 to avoid problems with cos(lat) being negative
  const double pole = ELM_PI/2.0;
  const double offset_pole = pole - lat_epsilon;
  
  assert((std::abs(lat) >= (pole + lat_epsilon)) && 
    "lat must be less than pi/2 within a small tolerance");
  assert((std::abs(decl) >= pole) && "decl must be strictly less than pi/2");

  // Ensure that latitude isn't too close to pole, to avoid problems with cos(lat) being negative
  double my_lat = std::min(offset_pole, std::max(1.0 * offset_pole, lat));
  double temp = -(sin(my_lat)*sin(decl)) / (cos(my_lat) * cos(decl));
  temp = std::min(1.0, std::max(-1.0, temp));
  
  return 2.0 * secs_per_radian * acos(temp);
}


// Initialize maximum daylength, based on latitude and maximum declination
// maximum declination hardwired for present-day orbital parameters, 
// +/- 23.4667 degrees = +/- 0.409571 radians, use negative value for S. Hem
double max_daylength(const double& lat) {
 return (lat < 0.0) ? daylength(lat, -0.409571) : daylength(lat, 0.409571);
}

} // namespace ELM
