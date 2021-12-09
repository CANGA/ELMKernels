/*! \file daylength.h
\brief Functions derived from DaylengthMod.F90
*/
#pragma once

namespace ELM {

// Computes daylength [seconds]
// Latitude and solar declination angle should both be specified in radians. decl must
// be strictly less than pi/2; lat must be less than pi/2 within a small tolerance.
// lat  [double]    latitude [radians]
// decl [double]    solar declination angle [radians]
double daylength(const double& lat, const double& decl);


// compute maximum daylength [seconds]
// based on latitude and maximum declination
// maximum declination hardwired for present-day orbital parameters, 
// +/- 23.4667 degrees = +/- 0.409571 radians, use negative value for S. Hem
// lat  [double]    latitude [radians]
double max_daylength(const double& lat);

} // namespace ELM
