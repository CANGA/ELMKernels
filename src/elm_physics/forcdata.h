
#pragma once

#include <cmath>
#include <string>
#include <array>
#include <tuple>

#include "elm_constants.h"
#include "mpi_types.hh"
#include "utils.hh"
#include "array.hh"
#include "read_input.hh"


#ifdef HAVE_PNETCDF
#include "read_pnetcdf.hh"
#else
#include "read_netcdf.hh"
#endif

#include "forcing_physics.h"

/*
Class to read, parse, and operate on atmospheric forcing input



instantaneous measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
are assumed to be point measurements and are interpolated for use in this model

continuous flux data like FSDS and PREC are assumed to be averaged over their timesteps
and are used without interpoating


ELM uses two different time schemes for forcing data. Instantaneous values are defined at 
t=forc_dt/2, t=3forc_dt/2 etc. Flux values are defined at t=forc_dt, t=2forc_dt etc.

We choose to define all forcing inputs for this model at multiples of t=forc_dt.
   
POINT DATA
 t_idx                       0      1      2       3
 file time                   0     dt     2dt     3dt
 centered start time         0     dt     2dt     3dt   -- this model, also ELM flux data
 staggered start time    -dt/2   dt/2   3dt/2   5dt/2   -- this is how ELM structures point data

  staggered S[]; centered C[];

      -fdt/2     fdt/2    fdt3/2            - S forcing times
              t=0       t=fdt     t=2fdt    - C forcing times
          X----|----X----|----X----|      
        S[0]       S[1]      S[2]
              C[0]      C[1]      C[2]

 S        X---------X---------X       alignment of forcing data
 C             X---------X---------X
              /           \
             /             \ physics timestepping with 4 model_dt per forc_dt 
            /               \  
        t=0/                 \ t = forc_dt = 4 * model_dt
          |-xx-|-xx-|-xx-|-xx-|
          C0                  C1
                                      
   We interpolate the ingested forcing data at the midpoint of the model timestep, 
   model_step_interp_time = model_timestep_start + model_dt / 2

*/



namespace ELM::forc_data_manager {

//using forcDataType = forcing_physics::forcDataType;
using forcDataType = ELMconstants::forcDataType;
//enum class forcDataType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };

// return name associated with enum type
template<forcDataType type>
constexpr auto vname ();

// return reference to variable that maps to dim_idx
template<typename T, typename U> 
constexpr T& get_dim_ref(const U dim_idx, T& t, T& x, T& y);

// return reference to variable that maps to dim_idx
template<typename T, typename U> 
constexpr T& get_dim_ref(const U dim_idx, T& t, T& x);


template<typename ArrayD1, typename ArrayD2, forcDataType type>
class ForcData {

private:
  ArrayD2 data_;
  GO ntimes_{0}, ncells_{0};
  double forc_dt_{0.0};
  std::string fname_, varname_;
  Utils::Date data_start_time_;
  Utils::Date file_start_time_;

public:
  constexpr ForcData(const std::string& filename, const Utils::Date &file_start_time, const GO ntimes, const GO ncells);

  // interface to update forcing file info
  constexpr void update_file_info(const Utils::Date& new_file_start_time, const std::string& new_filename);

  // interface to update working data start time
  constexpr void update_data_start_time(const GO t_idx);

  // calculate t_idx at model_time and check bounds
  // assumes model_time is centered on the model_dt interval, ie  = model_step_start + model_dt/2
  constexpr GO forc_t_idx_checks(const double& model_dt, const Utils::Date& model_time, const Utils::Date& forc_record_start_time) const;

  // calculate t_idx at model_time relative to forc_record_start_time
  constexpr GO forc_t_idx(const Utils::Date& model_time, const Utils::Date& forc_record_start_time) const;

  // calculate linear interpolation of [t1,t2] interval at t = model_time
  // only used for instantaneous point measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
  // !! ASSUMES instantaneous forcing measured at time corresponding to t_idx
  // ie given a model_t bounded by the forcing data interval [lb, ub], 
  // lb_time = data_start_time_ + forc_dt * t_idx and ub_time = data_start_time_ + forc_dt * (t_idx + 1)
  // forc_data_times_of_measurement =  {0, forc_dt, ..., Nforc_dt}
  // the other option is to define the values staggered by +- forc_dt/2
  constexpr std::pair<double,double> forcing_time_weights(const GO t_idx, const Utils::Date& model_time) const;

  // return a tuple of references to the passed in parameters based on the ordering of the file array
  // requires data(ntimes, nlon * nlat) and file_data(*,*,*) in {ntimes,nlon,nlat}
  template<typename T>
  constexpr auto input_idx_order(const Comm_type& comm, T& t, T& x, T& y) const;

  // requires data(ntimes, ncells) and file_data(ntimes, ncells)or(ncells, ntimes)
  template<typename T>
  constexpr auto input_idx_order(const Comm_type& comm, T& t, T& x) const;

  // read forcing data from a file 
  constexpr void read_atm_forcing(const Utils::DomainDecomposition<2> &dd, const Utils::Date& model_time, const GO ntimes);

  // get forcing data for the current timestep
  // interpolate point values
  // process data
  template<typename... Args>
  constexpr void get_forcing(const double& model_dt, const Utils::Date& model_time, Args&&...args);

};

} // namespace ELM::forc_data_manager

#include "forcdata_impl.hh"
