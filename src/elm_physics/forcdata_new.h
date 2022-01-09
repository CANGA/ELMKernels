#include <type_traits>
#include <cmath>
#include <string>
#include <array>
#include <tuple>

#include <iostream>

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


CURRENT ASSUMPTION: 




*/

  // calculate linear interpolation of [t1,t2] interval at t = model_time
  // only used for instantaneous point measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
  // !! ASSUMES instantaneous forcing measured at time corresponding to t_idx
  // ie given a model_t bounded by the forcing data interval [lb, ub], 
  // lb_time = data_start_time_ + forc_dt * t_idx and ub_time = data_start_time_ + forc_dt * (t_idx + 1)
  // forc_data_times_of_measurement =  {0, forc_dt, ..., Nforc_dt}

  // the other option is to define the values staggered by +- forc_dt/2
  // 
  // t_idx                       0      1      2       3
  // file time                   0     dt     2dt     3dt
  // centered start time         0     dt     2dt     3dt   -- this model
  // staggered start time    -dt/2   dt/2   3dt/2   5dt/2   -- this is how ELM does it
  //
  //  staggered S[]; centered C[];
  //
  //      -fdt/2     fdt/2    fdt3/2            - S forcing times
  //              t=0       t=fdt     t=2fdt    - C forcing times
  //          X----|----X----|----X----|      
  //        S[0]       S[1]      S[2]
  //              C[0]      C[1]      C[2]
  //
  // S        X---------X---------X       alignment of forcing data
  // C             X---------X---------X
  //              /           \
  //             /             \ physics timestepping with 4 model_dt per forc_dt 
  //            /               \  
  //        t=0/                 \ t = forc_dt = 4 * model_dt
  //          |-xx-|-xx-|-xx-|-xx-|
  //          C0                  C1
  //                                      
  //   We interpolate the ingested forcing data at the midpoint of the model timestep, 
  //   model_step_interp_time = model_timestep_start + model_dt / 2
  //



namespace ELM::forc_data_manager {

using forcDataType = forcing_physics::forcDataType;

//enum class forcDataType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };

template<forcDataType type>
constexpr auto vname () {
  std::string name;
  return [&name] {
    if constexpr(type == forcDataType::TBOT) name = "TBOT";
    if constexpr(type == forcDataType::PBOT) name = "PSRF";
    if constexpr(type == forcDataType::QBOT) name = "QBOT";
    if constexpr(type == forcDataType::RH) name = "RH";
    if constexpr(type == forcDataType::FLDS) name = "FLDS";
    if constexpr(type == forcDataType::FSDS) name = "FSDS"; 
    if constexpr(type == forcDataType::PREC) name = "PRECTmms";
    if constexpr(type == forcDataType::WIND) name = "WIND";
    if constexpr(type == forcDataType::ZBOT) name = "ZBOT";
    return name;
  }();
}


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
  ForcData(const std::string& filename, const Utils::Date &file_start_time, const GO ntimes, const GO ncells) 
    : fname_(filename), file_start_time_(file_start_time), ntimes_(ntimes),
      ncells_(ncells), varname_(vname<type>()), data_(vname<type>(), ntimes, ncells)  { }

  // interface to update forcing file info
  void update_file_info(const Utils::Date& new_file_start_time, const std::string& new_filename) {
    file_start_time_ = new_file_start_time;
    fname_ = new_filename;
  }

  // interface to update working data start time
  void update_data_start_time(const GO t_idx) {
    data_start_time_ = file_start_time_;
    data_start_time_.increment_seconds(static_cast<int>(round(86400 * forc_dt_)) * t_idx);
  }

  // calculate t_idx at model_time and check bounds
  // assumes model_time is centered on the model_dt interval, ie  = model_step_start + model_dt/2
  GO forc_t_idx_checks(const double& model_dt, const Utils::Date& model_time, const Utils::Date& forc_record_start_time) const {
    const double delta_to_model_time = Utils::days_since(model_time, forc_record_start_time);
    const double halfdt = model_dt * 0.5;
    const double model_dt_start = delta_to_model_time - halfdt;
    const double model_dt_end = delta_to_model_time + halfdt;
    const int idx_model_halfdt = static_cast<int>(delta_to_model_time / forc_dt_);
    const double forc_dt_start = idx_model_halfdt * forc_dt_;
    const double forc_dt_end = (idx_model_halfdt + 1) * forc_dt_;
    
    assert(forc_dt_ > 0.0 && "forc_dt is <= 0.0");
    constexpr double eps = 1.0e-8;
    assert((model_dt_start + eps * forc_dt_) >= 0.0 && "difference in dates is negative");
    
    assert((model_dt_start + eps * forc_dt_) >= forc_dt_start && 
      "model_dt start is less than forc_dt start, even with a 1e-8 * forc_dt error tolerance");
    assert((model_dt_end - eps * forc_dt_) <= forc_dt_end && 
      "model_dt end is greater than forc_dt end, even with a 1e-8 * forc_dt error tolerance");

    return idx_model_halfdt;
  }

  // calculate t_idx at model_time relative to forc_record_start_time
  GO forc_t_idx(const Utils::Date& model_time, const Utils::Date& forc_record_start_time) const {
    const double delta = Utils::days_since(model_time, forc_record_start_time);
    return static_cast<GO>( delta / forc_dt_);
  }

  // calculate linear interpolation of [t1,t2] interval at t = model_time
  // only used for instantaneous point measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
  // !! ASSUMES instantaneous forcing measured at time corresponding to t_idx
  // ie given a model_t bounded by the forcing data interval [lb, ub], 
  // lb_time = data_start_time_ + forc_dt * t_idx and ub_time = data_start_time_ + forc_dt * (t_idx + 1)
  // forc_data_times_of_measurement =  {0, forc_dt, ..., Nforc_dt}
  // the other option is to define the values staggered by +- forc_dt/2
  std::pair<double,double> forcing_time_weights(const GO t_idx, const Utils::Date& model_time) const {
    Utils::Date forc_start(data_start_time_);
    forc_start.increment_seconds(static_cast<int>(round(86400 * forc_dt_) * t_idx));
    const double elapsed_dt = Utils::days_since(model_time, forc_start) / forc_dt_;
    assert(elapsed_dt <= 1.0 && "time weights sum to > 1; model_time is greater than forc_t_start + forc_dt");
    assert(elapsed_dt >= 0.0 && "time weights sum to < 0");
    return std::make_pair(1.0 - elapsed_dt, elapsed_dt);
  }

  // return reference to variable that maps to dim_idx
  template<typename T, typename U> 
  T& get_dim_ref(const U dim_idx, T& t, T& x, T& y) { 
    if (dim_idx == 2) { return y; }
    if (dim_idx == 1) { return x; }
    if (dim_idx == 0) { return t; }
    throw std::runtime_error("ELM ERROR: NetCDF variable dimension index not in {0,1,2}");
  }
  
  // return reference to variable that maps to dim_idx
  template<typename T, typename U> 
  T& get_dim_ref(const U dim_idx, T& t, T& x) { 
    if (dim_idx == 1) { return x; }
    if (dim_idx == 0) { return t; }
    throw std::runtime_error("ELM ERROR: NetCDF variable dimension index not in {0,1}");
  }

  // return a tuple of references to the passed in parameters based on the ordering of the file array
  // requires data(ntimes, nlon * nlat) and file_data(*,*,*) in {ntimes,nlon,nlat}
  template<typename T>
  auto input_idx_order(const Comm_type& comm, T& t, T& x, T& y) {
    const auto dimids = 
        ELM::IO::get_vardimids<3>(comm, fname_, varname_);
  
    auto dim_id_idx = [=, &dimids] (const std::string& dimname) { 
        const auto dim_id = ELM::IO::get_dimid(comm, fname_, dimname);
        const auto itr = std::find(dimids.begin(), dimids.end(), dim_id);
        return std::distance(dimids.begin(), itr); };
  
    const auto i = dim_id_idx("DTIME");
    const auto j = dim_id_idx("lon");
    const auto k = dim_id_idx("lat");
    T& ii = get_dim_ref(i, t, x, y);
    T& jj = get_dim_ref(j, t, x, y);
    T& kk = get_dim_ref(k, t, x, y);
    return std::forward_as_tuple(ii, jj, kk);
  }

  // requires data(ntimes, ncells) and file_data(ntimes, ncells)or(ncells, ntimes)
  template<typename T>
  auto input_idx_order(const Comm_type& comm, T& t, T& x) {
    const auto dimids = 
        ELM::IO::get_vardimids<2>(comm, fname_, varname_);
  
    const auto dim_id_idx = [=, &dimids] (const std::string& dimname) { 
        const auto dim_id = ELM::IO::get_dimid(comm, fname_, dimname);
        const auto itr = std::find(dimids.begin(), dimids.end(), dim_id);
        return std::distance(dimids.begin(), itr); };
  
    const auto i = dim_id_idx("DTIME");
    const auto j = dim_id_idx("ncell"); 
    T& ii = get_dim_ref(i, t, x);
    T& jj = get_dim_ref(j, t, x);
    return std::forward_as_tuple(ii, jj);
  }

  // read forcing data from a file 
  void read_atm_forcing(const Utils::DomainDecomposition<2> &dd, const double& model_dt, const Utils::Date& model_time, const GO ntimes) {
    // resize if ntimes has changed - assume ncells_ doesn't change
    if (ntimes != data_.extent(0)) { ntimes_ = ntimes; data_.resize(ntimes, ncells_); }

    { // get forc_dt_ by differencing the first and second timestep
      // assume forc_dt doesn't change until next read
      const std::array<GO, 1> start = {0};
      const std::array<GO, 1> count = {2};
      ELM::Array<double, 1> arr_for_dt_measurement(2);
      IO::read(dd.comm, this->fname_, "DTIME", start, count, arr_for_dt_measurement.data());
      forc_dt_ = arr_for_dt_measurement(1) - arr_for_dt_measurement(0);
    }

    // get forcing time series time index (from file start time) immediately prior to model_time
    const auto file_t_idx = forc_t_idx_checks(model_dt, model_time, file_start_time_);
    update_data_start_time(file_t_idx);
    assert(data_.extent(0) == ntimes);
    assert(data_.extent(1) == dd.n_local[0] * dd.n_local[1]);

    // maps data_(ntimes, ncells) = arr_for_read(ii, jj, kk)
    // where (ii, jj, kk) are references to some arbitrary permutation of {ntimes, nlongitude, nlatitude}
    // get references to file array start indices
    const auto [si, sj, sk] = input_idx_order(dd.comm, file_t_idx, dd.start[0], dd.start[1]);
    std::array<GO, 3> start = {si, sj, sk};

    // get references to file array size
    auto [ci, cj, ck] = input_idx_order(dd.comm, ntimes, dd.n_local[0], dd.n_local[0]);
    std::array<GO, 3> count = {ci, cj, ck};
    ELM::Array<double, 3> arr_for_read(ci, cj, ck);

    // read data from file
    IO::read(dd.comm, this->fname_, this->varname_, start, count, arr_for_read.data());

    // get references to loop indices 
    GO i, j, k;
    auto [ii, jj, kk] = input_idx_order(dd.comm, i, j, k);
    // copy file data into model host array
    for (i = 0; i != ntimes; ++i) {
      for (j = 0; j != dd.n_local[0]; ++j) {
        for (k = 0; k != dd.n_local[1]; ++k) {
          data_(i, j * dd.n_local[1] + k) = arr_for_read(ii, jj, kk);
        }
      }
    }
  }

  // get forcing data for the current timestep
  // interpolate point values
  // process data
  template<typename... Args>
  constexpr void get_forcing(const double& model_dt, const Utils::Date& model_time, Args&&...args) {
  
    const GO t_idx = forc_t_idx_checks(model_dt, model_time, data_start_time_);
    const auto [wt1, wt2] = forcing_time_weights(t_idx, model_time);

    auto thisDataObject = [&] {
      if constexpr(type == forcDataType::TBOT) {
        return forcing_physics::ProcessTBOT(t_idx, wt1, wt2, data_, std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::PBOT) {
        return forcing_physics::ProcessPBOT(t_idx, wt1, wt2, data_, std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::QBOT || type == forcDataType::RH) {
        return forcing_physics::ProcessQBOT<ArrayD1, ArrayD2, type>(t_idx, wt1, wt2, data_, std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::FLDS) {
        return forcing_physics::ProcessFLDS(t_idx, wt1, wt2, data_, std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::FSDS) {
        return forcing_physics::ProcessFSDS(data_[t_idx], std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::PREC) {
        return forcing_physics::ProcessPREC(data_[t_idx], std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::WIND) {
        return forcing_physics::ProcessWIND(t_idx, wt1, wt2, data_, std::forward<Args>(args)...);
      } else if constexpr (type == forcDataType::ZBOT) {
        return forcing_physics::ProcessZBOT(std::forward<Args>(args)...);
      }
    }();

    for (int i = 0; i != ncells_; ++i) {
      std::invoke(thisDataObject, i);
    }
  }

};

} // namespace ELM::forc_data_manager


