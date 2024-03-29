
#pragma once

namespace ELM::atm_utils {

// return name associated with enum type
template <AtmForcType ftype>
constexpr auto get_varname() {
  return [] {
    std::string name;
    if constexpr (ftype == AtmForcType::TBOT) {
      name = {"TBOT"};
    }
    if constexpr (ftype == AtmForcType::PBOT) {
      name = {"PSRF"};
    }
    if constexpr (ftype == AtmForcType::QBOT) {
      name = {"QBOT"};
    }
    if constexpr (ftype == AtmForcType::RH) {
      name = {"RH"};
    }
    if constexpr (ftype == AtmForcType::FLDS) {
      name = {"FLDS"};
    }
    if constexpr (ftype == AtmForcType::FSDS) {
      name = {"FSDS"};
    }
    if constexpr (ftype == AtmForcType::PREC) {
      name = {"PRECTmms"};
    }
    if constexpr (ftype == AtmForcType::WIND) {
      name = {"WIND"};
    }
    if constexpr (ftype == AtmForcType::ZBOT) {
      name = {"ZBOT"};
    }
    return name;
  }();
}

// return reference to variable that maps to dim_idx for a file_array with 3 dimensions
template <typename T, typename U>
constexpr T& get_dim_ref(const U& dim_idx, T& t, T& x, T& y) {
  return (dim_idx == 2) ? y : (dim_idx == 1) ? x : (dim_idx == 0) ? t :
  throw std::runtime_error("ELM ERROR: NetCDF variable dimension index not in {0,1,2}");
}

// return reference to variable that maps to dim_idx for a file_array with 2 dimensions
template <typename T, typename U>
constexpr T& get_dim_ref(const U& dim_idx, T& t, T& x) {
  return (dim_idx == 1) ? x : (dim_idx == 0) ? t :
  throw std::runtime_error("ELM ERROR: NetCDF variable dimension index not in {0,1}");
}

} // namespace ELM::atm_utils

namespace ELM {
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr AtmDataManager<ArrayD1, ArrayD2, ftype>::
AtmDataManager(const std::string& filename,
               const Utils::Date& file_start_time,
               const size_t& ntimes, const size_t& ncells)
    : data(atm_utils::get_varname<ftype>(), ntimes, ncells),
      varname_{atm_utils::get_varname<ftype>()},
      fname_{filename},
      file_start_time_{file_start_time},
      ntimes_{ntimes},
      ncells_{ncells},
      data_start_time_{}
    {}

// interface to update forcing file info
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr void AtmDataManager<ArrayD1, ArrayD2, ftype>::
update_file_info(const Utils::Date& new_file_start_time,
                 const std::string& new_filename)
{
  file_start_time_ = new_file_start_time;
  fname_ = new_filename;
}

// interface to update working data start time
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr void AtmDataManager<ArrayD1, ArrayD2, ftype>::
update_data_start_time(const size_t& t_idx)
{
  data_start_time_ = file_start_time_;
  data_start_time_.increment_seconds(static_cast<int>(round(86400.0 * forc_dt_ * t_idx)));
}

// interface to return date of working data start time
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr Utils::Date AtmDataManager<ArrayD1, ArrayD2, ftype>::
get_data_start_time()
{
  return data_start_time_;
}

// interface to return forc_dt_ in days
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr double AtmDataManager<ArrayD1, ArrayD2, ftype>::
get_forc_dt_days()
{
  return forc_dt_;
}

// interface to return forc_dt_ in seconds
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr double AtmDataManager<ArrayD1, ArrayD2, ftype>::
get_forc_dt_secs()
{
  return 86400.0 * forc_dt_;
}

// calc t index when model_time == a forcing time interval
// assumes model_time is start of model_dt
// included for testing
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr size_t
AtmDataManager<ArrayD1, ArrayD2, ftype>::
forc_t_idx_aligned(const double& delta_to_model_time, const double& model_dt,
                   const Utils::Date& model_time,
                   const Utils::Date& forc_record_start_time) const
{
  const double model_dt_start = delta_to_model_time;
  const double model_dt_end = model_dt_start + model_dt;
  const int forc_t_idx = static_cast<int>(round(delta_to_model_time / forc_dt_));
  const double forc_dt_start = forc_t_idx * forc_dt_;
  const double forc_dt_end = (forc_t_idx + 1) * forc_dt_;
  constexpr double eps = 1.0e-8;
  assert((model_dt_start + eps * forc_dt_) >= forc_dt_start &&
         "model_dt start is less than forc_dt start, even with a 1e-8 * forc_dt error tolerance");
  assert((model_dt_end - eps * forc_dt_) <= forc_dt_end &&
         "model_dt end is greater than forc_dt end, even with a 1e-8 * forc_dt error tolerance");
  return forc_t_idx;
}

// calculate t_idx at model_time and check bounds
// model timesteps must fall within a single forcing time interval
// assumes model_time is centered on the model_dt interval, ie  = model_step_start + model_dt/2
// includes a (probably temporary) functionality to correctly calc the index in situations where
// model_time == a forcing time interval -- this is included for testing
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr size_t
AtmDataManager<ArrayD1, ArrayD2, ftype>::
forc_t_idx_check_bounds(const double& model_dt, const Utils::Date& model_time,
                        const Utils::Date& forc_record_start_time) const
{
  const double delta_to_model_time = Utils::days_since(model_time, forc_record_start_time);
  const bool aligned = (abs(remainder(delta_to_model_time, forc_dt_)) < 1.e-8) ? true : false;
  if (aligned) {
    return forc_t_idx_aligned(delta_to_model_time, model_dt, model_time, forc_record_start_time);
  }
  const double halfdt = model_dt * 0.5;
  const double model_dt_start = delta_to_model_time - halfdt;
  const double model_dt_end = delta_to_model_time + halfdt;
  const int forc_t_idx = static_cast<int>(delta_to_model_time / forc_dt_);
  const double forc_dt_start = forc_t_idx * forc_dt_;
  const double forc_dt_end = (forc_t_idx + 1) * forc_dt_;
  constexpr double eps = 1.0e-8;
  assert(forc_dt_ > 0.0 && "forc_dt is <= 0.0");
  assert((model_dt_start + eps * forc_dt_) >= 0.0 && "difference in dates is negative");
  assert((model_dt_start + eps * forc_dt_) >= forc_dt_start &&
         "model_dt start is less than forc_dt start, even with a 1e-8 * forc_dt error tolerance");
  assert((model_dt_end - eps * forc_dt_) <= forc_dt_end &&
         "model_dt end is greater than forc_dt end, even with a 1e-8 * forc_dt error tolerance");
  return forc_t_idx;
}

// calculate t_idx at model_time relative to forc_record_start_time
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr size_t AtmDataManager<ArrayD1, ArrayD2, ftype>::
forc_t_idx(const Utils::Date& model_time,
           const Utils::Date& forc_record_start_time) const
{
  const double delta = Utils::days_since(model_time, forc_record_start_time);
  return static_cast<size_t>(delta / forc_dt_);
}

// calculate linear interpolation of [t1,t2] interval at t = model_time
// only used for instantaneous point measurement data like TBOT, QBOT, PBOT, RH, FLDS, WIND
// !! ASSUMES instantaneous forcing measured at time corresponding to t_idx
// ie given a model_t bounded by the forcing data interval [lb, ub],
// lb_time = data_start_time_ + forc_dt * t_idx and ub_time = data_start_time_ + forc_dt * (t_idx + 1)
// forc_data_times_of_measurement =  {0, forc_dt, ..., Nforc_dt}
// the other option is to define the values staggered by +- forc_dt/2
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
constexpr std::pair<double, double>
AtmDataManager<ArrayD1, ArrayD2, ftype>::
forcing_time_weights(const size_t& t_idx, const Utils::Date& model_time) const
{
  Utils::Date forc_start(data_start_time_);
  forc_start.increment_seconds(static_cast<int>(round(86400 * forc_dt_) * t_idx));
  const double elapsed_dt = Utils::days_since(model_time, forc_start) / forc_dt_;
  assert(elapsed_dt <= 1.0 && "time weights sum to > 1; model_time is greater than forc_t_start + forc_dt");
  assert(elapsed_dt >= 0.0 && "time weights sum to < 0");
  return std::make_pair(1.0 - elapsed_dt, elapsed_dt);
}

// return reference to arg in Args that matches position of dimension dimname in file array
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
template <typename... Args, size_t D>
constexpr auto&
AtmDataManager<ArrayD1, ArrayD2, ftype>::
get_ref_to_dim(const std::string& dimname, const Comm_type& comm,
               const std::array<int, D>& dimids, Args&&...args) const
{
  auto dimname_idx = [this, comm, &dimids](const std::string& dimname) {
    const auto dim_id = ELM::IO::get_dimid(comm, fname_, dimname);
    const auto itr = std::find(dimids.begin(), dimids.end(), dim_id);
    return std::distance(dimids.begin(), itr);
  };

  const auto dim_idx = dimname_idx(dimname);
  return ELM::atm_utils::get_dim_ref(dim_idx, std::forward<Args>(args)...);
}

// return a tuple of references to the passed in parameters based on the ordering of the file array
// requires data(ntimes, nlon * nlat) and file_data(*,*,*) in {ntimes,nlon,nlat}
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
template <typename T>
constexpr auto AtmDataManager<ArrayD1, ArrayD2, ftype>::
order_inputs(const Comm_type& comm, T& t, T& x, T& y) const
{
  const auto dimids = ELM::IO::get_var_dimids<3>(comm, fname_, varname_);

  T& ii = get_ref_to_dim("DTIME", comm, dimids, t, x, y);
  T& jj = get_ref_to_dim("lon", comm, dimids, t, x, y);
  T& kk = get_ref_to_dim("lat", comm, dimids, t, x, y);
  return std::forward_as_tuple(ii, jj, kk);
}

// requires data(ntimes, ncells) and file_data(ntimes, ncells)or(ncells, ntimes)
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
template <typename T>
constexpr auto AtmDataManager<ArrayD1, ArrayD2, ftype>::
order_inputs(const Comm_type& comm, T& t, T& x) const
{
  const auto dimids = ELM::IO::get_var_dimids<2>(comm, fname_, varname_);

  T& ii = get_ref_to_dim(comm, "DTIME", dimids, t, x);
  T& jj = get_ref_to_dim(comm, "ncell", dimids, t, x);
  return std::forward_as_tuple(ii, jj);
}

// read forcing data from a file
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
template <typename h_ArrayD2>
constexpr void AtmDataManager<ArrayD1, ArrayD2, ftype>::
read_atm_forcing(h_ArrayD2 h_data, 
                 const Utils::DomainDecomposition<2>& dd,
                 const Utils::Date& model_time,
                 const size_t& ntimes)
{
  if (need_new_data_) {
    std::cout << "reading new "+atm_utils::get_varname<ftype>()+" data\n";
    // resize if ntimes has changed - assume ncells_ doesn't change
    if (ntimes != static_cast<size_t>(h_data.extent(0))) {
      ntimes_ = ntimes;
      NS::resize(h_data, ntimes, ncells_);
    }

    { // get forc_dt_ by differencing the first and second timestep
      // need to do this to init/update forc_dt_ before forc_t_idx() is called
      // assume forc_dt doesn't change until next read
      const std::array<size_t, 1> start = {0};
      const std::array<size_t, 1> count = {2};
      ELM::Array<double, 1> arr_for_dt_measurement(2);
      IO::read_netcdf(dd.comm, fname_, "DTIME", start, count, arr_for_dt_measurement.data());
      forc_dt_ = arr_for_dt_measurement(1) - arr_for_dt_measurement(0);
    }

    { // get scale factor and offset if available
      int err = IO::get_attribute(dd.comm, fname_, varname_, "scale_factor", scale_factor_);
      if (err < 0) {
        scale_factor_ = 1.0;
      }
      err = IO::get_attribute(dd.comm, fname_, varname_, "add_offset", add_offset_);
      if (err < 0) {
        add_offset_ = 0.0;
      }
    }

    // get forcing time series time index (from file start time) immediately prior to model_time
    const auto file_t_idx = forc_t_idx(model_time, file_start_time_);
    update_data_start_time(file_t_idx);
    // check data extents
    assert(static_cast<size_t>(h_data.extent(0)) == ntimes);
    assert(static_cast<size_t>(h_data.extent(1)) == dd.n_local[0] * dd.n_local[1]);

    // maps h_data(ntimes, ncells) = arr_for_read(ii, jj, kk)
    // where (ii, jj, kk) are references to some arbitrary permutation of {ntimes, nlongitude, nlatitude}
    // get references to file array start indices
    const auto [si, sj, sk] = order_inputs(dd.comm, file_t_idx, dd.start[0], dd.start[1]);
    std::array<size_t, 3> start = {si, sj, sk};

    // get references to file array size
    const auto [ci, cj, ck] = order_inputs(dd.comm, ntimes, dd.n_local[0], dd.n_local[1]);
    std::array<size_t, 3> count = {ci, cj, ck};

    // read data from file
    ELM::Array<double, 3> arr_for_read(ci, cj, ck);
    IO::read_netcdf(dd.comm, fname_, varname_, start, count, arr_for_read.data());

    // get references to loop indices
    size_t i, j, k;
    const auto [ii, jj, kk] = order_inputs(dd.comm, i, j, k);
    // copy file data into model host array
    for (i = 0; i != ntimes; ++i) {
      for (j = 0; j != dd.n_local[0]; ++j) {
        for (k = 0; k != dd.n_local[1]; ++k) {
          h_data(i, j * dd.n_local[1] + k) = arr_for_read(ii, jj, kk) * scale_factor_ + add_offset_;
        }
      }
    }
    need_new_data_ = false;
  } // if (need_new_data_)
}

// read forcing data from a file - update file info and call main read_atm method
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
template <typename h_ArrayD2>
constexpr void AtmDataManager<ArrayD1, ArrayD2, ftype>::
read_atm_forcing(h_ArrayD2 h_data, 
                 const Utils::DomainDecomposition<2>& dd,
                 const Utils::Date& model_time,
                 const size_t& ntimes,
                 const Utils::Date& new_file_start_time,
                 const std::string& new_filename)
{
  update_file_info(new_file_start_time, new_filename);
  read_atm_forcing(h_data, dd, model_time, ntimes);
}

// get forcing data for the current timestep
// interpolate point values
// process data
// assumes parameter model_time is model_start + model_dt/2
template <typename ArrayD1, typename ArrayD2, AtmForcType ftype>
template <typename... Args>
constexpr void AtmDataManager<ArrayD1, ArrayD2, ftype>::
get_atm_forcing(const double& model_dt,
                const Utils::Date& model_time,
                Args&&...args)
{
  const auto physics_object = [=, &args...] (const size_t& t_idx, const Utils::Date& model_time) {
    const auto [wt1, wt2] = forcing_time_weights(t_idx, model_time);
    if constexpr (ftype == AtmForcType::TBOT) {
      return atm_forcing_physics::
        ProcessTBOT(t_idx, wt1, wt2, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::PBOT) {
      return atm_forcing_physics::
        ProcessPBOT(t_idx, wt1, wt2, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::QBOT || ftype == AtmForcType::RH) {
      return atm_forcing_physics::
        ProcessQBOT<ArrayD1, ArrayD2, ftype>(t_idx, wt1, wt2, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::FLDS) {
      return atm_forcing_physics::
        ProcessFLDS(t_idx, wt1, wt2, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::FSDS) {
      return atm_forcing_physics::
        ProcessFSDS(t_idx, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::PREC) {
      return atm_forcing_physics::
        ProcessPREC(t_idx, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::WIND) {
      return atm_forcing_physics::
        ProcessWIND(t_idx, wt1, wt2, data, std::forward<Args>(args)...);
    } else if constexpr (ftype == AtmForcType::ZBOT) {
      return atm_forcing_physics::
        ProcessZBOT(std::forward<Args>(args)...);
    }
  };

  const size_t t_idx = forc_t_idx_check_bounds(model_dt, model_time, data_start_time_);
  const auto& forcing = physics_object(t_idx, model_time);
  apply_parallel_for(forcing, "ComputeAtmForcing_"+atm_utils::get_varname<ftype>(), static_cast<int>(ncells_));

  if (t_idx == ntimes_ - 2) // last possible calculation
    need_new_data_ = true;
}

template <typename ArrayD1, typename ArrayD2>
ELM::AtmForcObjects<ArrayD1, ArrayD2>::AtmForcObjects(
  const std::string& filename,
  const ELM::Utils::Date &file_start_time,
  const int atm_nsteps, const int ncells) :
    forc_TBOT(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::TBOT>>
      (filename, file_start_time, atm_nsteps, ncells)),
    forc_PBOT(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::PBOT>>
      (filename, file_start_time, atm_nsteps, ncells)),
    //forc_QBOT(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::QBOT>>
    //  (filename, file_start_time, atm_nsteps, ncells)),
    forc_QBOT(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::RH>>
      (filename, file_start_time, atm_nsteps, ncells)),
    forc_FLDS(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::FLDS>>
      (filename, file_start_time, atm_nsteps, ncells)),
    forc_FSDS(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::FSDS>>
      (filename, file_start_time, atm_nsteps, ncells)),
    forc_PREC(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::PREC>>
      (filename, file_start_time, atm_nsteps, ncells)),
    forc_WIND(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::WIND>>
      (filename, file_start_time, atm_nsteps, ncells)),
    forc_ZBOT(std::make_shared<ELM::AtmDataManager<ArrayD1, ArrayD2, AtmForcType::ZBOT>>
      (filename, file_start_time, atm_nsteps, ncells))
  {}

} // namespace ELM
