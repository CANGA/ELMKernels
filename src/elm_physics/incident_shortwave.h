
#pragma once

namespace ELM::incident_shortwave {
  
// declination angle calc from ats/landlab
// doy [int]      day of year
// returns delination angle [radians]
double declination_angle(int doy);

// declination angle calc from ELM lnd_import szenith()/shr_orb_cosz()
// doy [int]      integer day of year
// returns delination angle [radians]
double declination_angle2(int doy);

// cosine of the solar zenith angle
// latrad [double]    latitude [radians]
// lonrad [double]    longitude [radians]
// jday   [double]    Julian day of year
double coszen(double latrad, double lonrad, double jday);




// average coszen functions
// derived from shr_orb_avg_cosz() in shr_orb_mod.F90

// adjust variable (latitude or declination angle) so that its tangent will be defined
// var [double]    variable to be bounded
double ensure_tan_defined(double var);

// convert model dt from seconds to radians wrt daylength
// dt [double]    timestep [seconds]
// returns model dt [radians]
double dt_radians (double dt);

// define dt start time of day on the period -2*p1 to 2*pi
// jday   [double]    Julian day of year
// lonrad [double]    longitude [radians]
// returns dt start time of day [-2pi to 2pi]
double dt_start_rad (double jday, double lonrad);

// define time of day at end of dt 
// t_start [double]    start time of day [-2pi to 2pi]
// dtrad   [double]    timestep [radians]
// returns dt end time of day [-2pi to 2pi]
double dt_end_rad (double t_start, double dtrad);

// define the cosine of the half-day length [0 to pi]
// adjust for cases of all daylight or all night
// latrad [double]    latitude [radians]
// declin [double]    solar declination angle [radians]
double coshalfday (double latrad, double declin);

// define the hour angle
// force it to be between -cos_h and cos_h
// consider the situation when the night period is too short
// input:
// t_start        [double]    start time of day [-2pi to 2pi]
// t_end          [double]    end time of day [-2pi to 2pi]
// dtrad          [double]    timestep [radians]
// cos_h          [double]    cosine of the half-day length
// output:
// hour_angle[4]  [double]    hour angle parameters
void avg_hourangle(double t_start, double t_end, double dtrad, double cos_h, double *hour_angle);

// perform a time integration to obtain cosz if desired
// output is valid over the period from t to t + dt
// t_start  [double]    start time of day [-2pi to 2pi]
// t_end    [double]    end time of day [-2pi to 2pi]
// dtrad    [double]    timestep [radians]
// cos_h    [double]    cosine of the half-day length
// latrad   [double]    latitude [radians]
// declin   [double]    solar declination angle [radians]
// returns time-integrated cosine of the solar zenith angle
double integrate_cosz(double t_start, double t_end, double dtrad, double cos_h, 
  double latrad, double declin);


/* evaluate average cosine(zenith) for a given dt

    A New Algorithm for Calculation of Cosine Solar Zenith Angle
    Author: Linjiong Zhou
    E-mail: linjiongzhou@hotmail.com
    Date  : 2015.02.22
    Ref.  : Zhou et al., GRL, 2015

in:
latrad  [double] latitude [radians]
lonrad  [double] longititude [radians]
declin  [double] solar declination angle [radians]
dt      [double] time increment over which to average cosz [seconds]
jday    [double] Julian date at the beginning of dt

out:
cosine of the solar zenith angle averaged over dt
*/
double average_cosz(double latrad, double lonrad, double declin, double dt, double jday);

} // namespace ELM::incident_shortwave

