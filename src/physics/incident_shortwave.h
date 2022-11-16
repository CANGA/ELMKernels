
#pragma once

// there are three methods to calculate zenith angle
// they all produce similar results for the lat/lon tested here.

// for now a single value for coszen is appropriate
// but a slope based factor would necessitate per-cell values

// first method - average cosz for dt_start to dt_end
//    auto decday = ELM::Utils::decimal_doy(current) + 1.0;
//    assign(S->coszen, ELM::incident_shortwave::average_cosz(S->lat_r, S->lon_r, dtime, decday));

// second method - point cosz at dt_start + dt/2
//    auto thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, decday + dtime / 86400.0 /2.0);

// third method - calc avg cosz over forcing dt (larger than model dt)
// then calc point dt at start + dt/2
// and use to calculate cosz_factor
//    ELM::Utils::Date forc_dt_start{forc_FSDS.get_data_start_time()};
//    forc_dt_start.increment_seconds(round(forc_FSDS.forc_t_idx(time_plus_half_dt, forc_FSDS.get_data_start_time()) * forc_FSDS.get_forc_dt_secs()));
//    double cosz_forc_decday = ELM::Utils::decimal_doy(forc_dt_start) + 1.0;
//    auto cosz_forcdt_avg = ELM::incident_shortwave::average_cosz(lat_r, lon_r, forc_FSDS.get_forc_dt_secs(), cosz_forc_decday);
//    auto thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, decday + dtime / 86400.0 /2.0);
//    cosz_factor = (thiscosz > 0.001) ? std::min(thiscosz/cosz_forcdt_avg, 10.0) : 0.0;

namespace ELM::incident_shortwave {

// declination angle calc from ats/landlab
// doy [int]      day of year
// returns delination angle [radians]
double declination_angle_cos(const int& doy);

// declination angle calc from ELM lnd_import szenith()/shr_orb_cosz()
// doy [int]      integer day of year
// returns delination angle [radians]
double declination_angle_sin(const int& doy);

// cosine of the solar zenith angle
// latrad [double]    latitude [radians]
// lonrad [double]    longitude [radians]
// jday   [double]    Julian day of year
double coszen(const double& latrad, const double& lonrad, const double& jday);

// average coszen functions
// derived from shr_orb_avg_cosz() in shr_orb_mod.F90

// adjust variable (latitude or declination angle) so that its tangent will be defined
// var [double]    variable to be bounded
double ensure_tan_defined(const double& var);

// convert model dt from seconds to radians wrt daylength
// dt [double]    timestep [seconds]
// returns model dt [radians]
double dt_radians(const double& dt);

// define dt start time of day on the period -2*p1 to 2*pi
// jday   [double]    Julian day of year
// lonrad [double]    longitude [radians]
// returns dt start time of day [-2pi to 2pi]
double dt_start_rad(const double& jday, const double& lonrad);

// define time of day at end of dt
// t_start [double]    start time of day [-2pi to 2pi]
// dtrad   [double]    timestep [radians]
// returns dt end time of day [-2pi to 2pi]
double dt_end_rad(const double& t_start, const double& dtrad);

// define the cosine of the half-day length [0 to pi]
// adjust for cases of all daylight or all night
// latrad [double]    latitude [radians]
// declin [double]    solar declination angle [radians]
double coshalfday(const double& latrad, const double& declin);

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
void avg_hourangle(const double& t_start, const double& t_end, const double& dtrad,
                   const double& cos_h, double hour_angle[4]);

// perform a time integration to obtain cosz if desired
// output is valid over the period from t to t + dt
// t_start  [double]    start time of day [-2pi to 2pi]
// t_end    [double]    end time of day [-2pi to 2pi]
// dtrad    [double]    timestep [radians]
// cos_h    [double]    cosine of the half-day length
// latrad   [double]    latitude [radians]
// declin   [double]    solar declination angle [radians]
// returns time-integrated cosine of the solar zenith angle
double integrate_cosz(const double& t_start, const double& t_end, const double& dtrad,
                      const double& cos_h, const double& latrad, const double& declin);

/* evaluate average cosine(zenith) for a given dt

    A New Algorithm for Calculation of Cosine Solar Zenith Angle
    Author: Linjiong Zhou
    E-mail: linjiongzhou@hotmail.com
    Date  : 2015.02.22
    Ref.  : Zhou et al., GRL, 2015

in:
latrad  [double] latitude [radians]
lonrad  [double] longititude [radians]
dt      [double] time increment over which to average cosz [seconds]
jday    [double] Julian date at the beginning of dt

out:
cosine of the solar zenith angle averaged over dt
*/
double average_cosz(const double& latrad, const double&  lonrad, const double&  dt, const double& jday);

} // namespace ELM::incident_shortwave
