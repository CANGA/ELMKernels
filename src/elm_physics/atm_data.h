
#pragma once

#include "array.hh"
#include "atm_physics.h"
#include "elm_constants.h"
#include "read_input.hh"
#include "utils.hh"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <map>

#include "kokkos_includes.hh"
#include "invoke_kernel.hh"


/*
Class to read, parse, and operate on atmospheric forcing input


instantaneous measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
are assumed to be point measurements and are interpolated for use in this model

continuous flux data like FSDS and PREC are assumed to be averaged over their timesteps
and are used without interpolating

ELM uses two different time schemes for forcing data. Instantaneous values are defined at
a forc_dt/2 offset from t=0 -- forc_T(n) (n=0) = -forc_dt/2, forc_dt/2, 3forc_dt/2, ..., (n)forc_dt - forc_dt/2

Continuous flux values in ELM are defined without any offset -- forc_T0 = 0, forc_T1 = forc_dt, forc_T2 = 2forc_dt.

We choose to define all forcing inputs for this model at multiples of t=forc_dt, the way ELM defines continuous data.


EXAMPLE TABLE OF FORCING DATA ORIENTATION
ELM time indexing for point data compared to the indexing method used in this model.
with model_dt = 1800s and forc_dt = 3600s:

                                      ELM                       THIS MODEL (equivalent to ELM for continuous data)
nstep t_start t_end ELM_forc_lb_t ELM_forc_lb_t forc_t_idices | NEW_forc_lb_t NEW_forc_lb_t forc_t_idices
1        0     1800         -1800          1800        (0,1)  |             0          3600         (0,1)
2     1800     3600          1800          5400        (1,2)  |             0          3600         (0,1)
3     3600     5400          1800          5400        (1,2)  |          3600          7200         (1,2)
4     5400     7200          5400          9000        (2,3)  |          3600          7200         (1,2)
5     7200     9000          5400          9000        (2,3)  |          7200         10800         (2,3)
|


VISUAL REPRESENTATION OF FORCING TIME ORIENTATION

  staggered S[]; -- ELM point data
  centered C[]; -- this model

      -fdt/2      fdt/2     fdt3/2          - S forcing times
              t=0      t=fdt    t=2fdt      - C forcing times
          X----|----X----|----X----|
        S[0]       S[1]      S[2]
              C[0]      C[1]      C[2]

 S        X---------X---------X       alignment of forcing data
 C             X---------X---------X
              /           \
             /             \ physics timestepping with 4 model_dt per forc_dt
            /               \
        t=0/                 \ t = forc_dt = 4 * model_dt
 C        |-xx-|-xx-|-xx-|-xx-| model timesteps within 1 forc_dt
          C0                  C1

   We interpolate the ingested forcing data at the midpoint of the model timestep (xx),
   model_step_interp_time = model_timestep_start + model_dt / 2

*/

namespace ELM::atm_utils {

// return name associated with enum type
template <AtmForcType ftype> constexpr auto get_varname();

// return reference to variable that maps to dim_idx for a file_array with 3 dimensions
template <typename T, typename U> constexpr T& get_dim_ref(const U& dim_idx, T& t, T& x, T& y);

// return reference to variable that maps to dim_idx for a file_array with 2 dimensions
template <typename T, typename U> constexpr T& get_dim_ref(const U& dim_idx, T& t, T& x);

} // namespace ELM::atm_utils


namespace ELM {
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
class AtmDataManager {

public:

  // public to provide access from driver
  ArrayD2 data; // 2D (ntimes, ncells) array-like object of forcing data -- device array

  constexpr AtmDataManager(const std::string& filename, const Utils::Date& file_start_time,
                           const size_t& ntimes, const size_t& ncells);
  ~AtmDataManager() = default;

  // interface to update forcing file info
  constexpr void update_file_info(const Utils::Date& new_file_start_time, const std::string& new_filename);

  // interface to update working data start time
  constexpr void update_data_start_time(const size_t& t_idx);

  // interface to return date of working data start time
  constexpr Utils::Date get_data_start_time();

  // interface to return forc_dt_ in days
  constexpr double get_forc_dt_days();

  // interface to return forc_dt_ in seconds
  constexpr double get_forc_dt_secs();

  // included for testing
  constexpr size_t forc_t_idx_aligned(const double& delta_to_model_time, const double& model_dt,
                                      const Utils::Date& model_time, const Utils::Date& forc_record_start_time) const;

  // calculate t_idx at model_time and check bounds
  // assumes model_time is centered on the model_dt interval, ie parameter model_time = model_step_start_time +
  // model_dt/2
  constexpr size_t forc_t_idx_check_bounds(const double& model_dt, const Utils::Date& model_time,
                                           const Utils::Date& forc_record_start_time) const;

  // calculate t_idx at model_time relative to forc_record_start_time
  constexpr size_t forc_t_idx(const Utils::Date& model_time, const Utils::Date& forc_record_start_time) const;

  // calculate linear interpolation of [t1,t2] interval at t = model_time
  // only used for instantaneous point measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
  // !! ASSUMES instantaneous forcing measured at time corresponding to t_idx
  // ie given a model_t bounded by the forcing data interval [lb, ub],
  // lb_time = data_start_time_ + forc_dt * t_idx and ub_time = data_start_time_ + forc_dt * (t_idx + 1)
  // forc_data_times_of_measurement =  {0, forc_dt, ..., Nforc_dt}
  // the other option is to define the values staggered by +- forc_dt/2
  constexpr std::pair<double, double> forcing_time_weights(const size_t& t_idx, const Utils::Date& model_time) const;

  // read forcing data from a file
  template <typename h_ArrayD2>
  constexpr void read_atm_forcing(h_ArrayD2 h_data, const Utils::DomainDecomposition<2>& dd, const Utils::Date& model_time,
                                  const size_t& ntimes);

  // read forcing data from a file - update file info and call main read_atm method
  template <typename h_ArrayD2>
  constexpr void read_atm_forcing(h_ArrayD2 h_data, const Utils::DomainDecomposition<2>& dd, const Utils::Date& model_time,
                                  const size_t& ntimes, const Utils::Date& new_file_start_time,
                                  const std::string& new_filename);

  // get forcing data for the current timestep
  // interpolate point values
  // process data
  template <typename... Args>
  constexpr void get_atm_forcing(const double& model_dt, const Utils::Date& model_time, Args&&...args);

private:
  // return reference to arg in Args that matches position of dimension dimname in file array
  template <typename... Args, size_t D>
  constexpr auto& get_ref_to_dim(const std::string& dimname, const Comm_type& comm, const std::array<int, D>& dimids,
                                 Args&&...args) const;

  // return a tuple of references to the passed in parameters based on the ordering of the file array
  // requires data(ntimes, nlon * nlat) and file_data(*,*,*) in {ntimes,nlon,nlat}
  template <typename T>
  constexpr auto order_inputs(const Comm_type& comm, T& t, T& x, T& y) const;

  // requires data(ntimes, ncells) and file_data(ntimes, ncells)or(ncells, ntimes)
  template <typename T>
  constexpr auto order_inputs(const Comm_type& comm, T& t, T& x) const;

  std::string varname_, fname_; // variable name and full file name including path - "src/xyz/file.nc"
  Utils::Date file_start_time_; // date object containing file dataset start time
  size_t ntimes_, ncells_;      // dimensions for data_
  Utils::Date data_start_time_; // start time of forcing read from file into data_
  double forc_dt_{0.0};        // forcing data timestep (days) - recalc at every read; assume constant between reads
  // the next two variables allow compatibility with ELM forcing data that is stored in a short int
  // format and then scaled and potentially added to
  double scale_factor_{1.0}; // factor for scaling input data - maybe needed when using some ELM input data
  double add_offset_{0.0};   // offset to add to input data - maybe needed when using some ELM input data
};

} // namespace ELM

#include "atm_data_impl.hh"
