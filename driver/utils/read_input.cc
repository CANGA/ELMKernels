//! A set of utilities for testing ELM kernels in C++
#ifndef ELM_KERNEL_TEST_NETCDF_HH_
#define ELM_KERNEL_TEST_NETCDF_HH_

#include "read_input.hh"
#include <sstream>

namespace ELM {
namespace IO {

//
// Readers for forcing data.
// -----------------------------------------------------------------------------

//
// Returns shape in a forcing file, as { N_TIMES, N_LAT_GLOBAL, N_LON_GLOBAL }
//
std::array<GO, 3> get_forcing_dimensions(const MPI_Comm &comm, const std::string &dir, const std::string &basename,
                                         const std::string &varname, const Utils::Date &time_start, int n_months) {
  std::array<GO, 3> total_dims = {0, 0, 0};

  // files organized by month
  Utils::Date current(time_start);
  for (int mm = 0; mm != n_months; ++mm) {
    auto date = current.date();
    int month = std::get<1>(date);
    int year = std::get<0>(date);

    std::stringstream fname_full;
    fname_full << dir << "/" << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";
    std::string fname = fname_full.str();
    auto dims = get_dimensions<3>(comm, fname, varname);

    if (mm == 0) {
      total_dims[1] = dims[1];
      total_dims[2] = dims[2];
    } else {
      assert(total_dims[1] == dims[1]);
      assert(total_dims[2] == dims[2]);
    }
    total_dims[0] += dims[0];
    current.increment_month();
  }
  return total_dims;
}

//
// Read a forcing file
//
// Requires shape(arr) == { N_TIMES, N_LAT_LOCAL, N_LON_LOCAL }
//
void read_forcing(const std::string &dir, const std::string &basename, const std::string &varname,
                  const Utils::Date &time_start, int n_months, const Utils::DomainDecomposition<2> &dd,
                  Array<double, 3> &arr) {
  Utils::Date current(time_start);
  int i_times = 0;

  // my slice in space
  std::array<GO, 3> start = {0, dd.start[0], dd.start[1]};
  std::array<GO, 3> count = {0, dd.n_local[0], dd.n_local[1]};

  for (int mm = 0; mm != n_months; ++mm) {
    auto date = current.date();
    int month = std::get<1>(date);
    int year = std::get<0>(date);

    std::stringstream fname_full;
    fname_full << dir << "/" << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";

    // std::cout << "reading: " << fname_full.str() << std::endl;

    // this file's slice in time is the full thing
    auto dims = get_dimensions<3>(dd.comm, fname_full.str(), varname);
    count[0] = dims[0];

    // read the slice, into a specific location
    // std::cout << "  at (" << start[0] << "," << start[1] << "," << start[2] << ");"
    //           << " (" << count[0] << "," << count[1] << "," << count[2] << ")"
    //           << std::endl;
    // std::cout << "  into: " << i_times << std::endl;
    read(dd.comm, fname_full.str(), varname, start, count, arr[i_times].data());

    // increment the position and month
    i_times += dims[0];
    current.increment_month();
  }
}

//
// Readers for phenology data.
// -----------------------------------------------------------------------------

//
// Returns shape in a phenology file, as { N_TIMES, N_PFTS, N_LAT_GLOBAL, N_LON_GLOBAL }
//
std::array<GO, 4> get_phenology_dimensions(const Comm_type &comm, const std::string &dir, const std::string &basename,
                                           const std::string &varname, const Utils::Date &time_start, int n_months) {
  Utils::Date current(time_start);
  current.increment_month(n_months);
  current.increment_day(-1); // back up one day to get last day of previous month
  assert(time_start.year == current.year && "No current support for crossing years in phenology data?");

  std::stringstream fname_full;
  fname_full << dir << "/" << basename;
  auto phen_dims_one = get_dimensions<4>(comm, fname_full.str(), varname);
  phen_dims_one[0] = n_months;
  return phen_dims_one;
}

//
// Readers for pft constants data.
// -----------------------------------------------------------------------------

//
// Returns maxpfts in a pft_constants file
//
std::array<GO, 1> get_maxpfts(const Comm_type &comm, const std::string &dir, const std::string &basename,
                              const std::string &varname) {
  std::stringstream fname_full;
  fname_full << dir << "/" << basename;
  auto maxpfts = get_dimensions<1>(comm, fname_full.str(), varname);
  return maxpfts;
}

//
// Read a phenology file.
//
// Requires shape(arr) == { N_TIMES, N_PFTS, N_LAT_LOCAL, N_LON_LOCAL }
//
void read_phenology(const std::string &dir, const std::string &basename, const std::string &varname,
                    const Utils::Date &time_start, int n_months, const Utils::DomainDecomposition<2> &dd,
                    Array<double, 4> &arr) {
  // my slice in space
  assert(n_months > 0);
  std::array<GO, 4> start = {(GO)(std::get<1>(time_start.date()) - 1), 0, dd.start[0], dd.start[1]};

  assert(arr.extent(1) > 0);
  std::array<GO, 4> count = {(GO)n_months, (GO)arr.extent(1), dd.n_local[0], dd.n_local[1]};

  // std::cout << "Reading: start = " << start[0] << "," << start[1] << "," << start[2] << "," << start[3] << std::endl
  //           << "         count = " << count[0] << "," << count[1] << "," << count[2] << "," << count[3] << std::endl;

  std::stringstream fname_full;
  fname_full << dir << "/" << basename;

  read(dd.comm, fname_full.str(), varname, start, count, arr.data());
  return;
}

} // namespace IO
} // namespace ELM

#endif
