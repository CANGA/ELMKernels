//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_UTILS_HH_
#define ELM_UTILS_HH_

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "array.hh"

namespace ELM {
namespace Utils {

//
// Domain decomposition in x, y
//
std::pair<int, int> get_domain_decomposition_2D(int nprocs, int argc, char **argv) {
  int numprocs_x = 1;
  int numprocs_y = 1;
  int n = nprocs;
  int i = 2;
  while (n != 1) {
    if (n % i == 0) {
      if (numprocs_x - numprocs_y >= 0) {
        numprocs_y *= i;
      } else {
        numprocs_x *= i;
      }
      n = n / i;
    } else {
      i++;
    }
  }
  return std::make_pair(numprocs_x, numprocs_y);
}

//
// { min, max, sum } of arrays
//
#ifdef HAVE_MPI
// hide this from bad developers!
namespace Impl {
#endif

template <typename Array_type> std::array<double, 3> min_max_sum(const Array_type &arr) {
  auto min_max = std::minmax_element(arr.begin(), arr.end());
  double min = *min_max.first;
  double max = *min_max.second;
  double sum = std::accumulate(arr.begin(), arr.end(), 0.);
  return std::array<double, 3>{min, max, sum};
}

#ifdef HAVE_MPI
} // namespace Impl

template <typename Array_type> std::array<double, 3> min_max_sum(const MPI_Comm &comm, const Array_type &arr) {
  auto lmms = Impl::min_max_sum(arr);
  std::array<double, 3> gmms;

  MPI_Reduce(&std::get<0>(lmms), &std::get<0>(gmms), 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&std::get<1>(lmms), &std::get<1>(gmms), 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&std::get<2>(lmms), &std::get<2>(gmms), 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  return gmms;
}

#endif

//
// Uses air temp to convert from total precip to rain and snow.
//
// Initially all precip is in rain, afterwards it is in both.
//
template <class Array_t> void convert_precip_to_rain_snow(Array_t &rain, Array_t &snow, const Array_t &temp) {
  int nt = rain.extent(0);
  int ng = rain.extent(1);
  for (int k = 0; k < nt; k++) {
    for (int j = 0; j < ng; j++) {
      if (temp(k, j) < 273.15) {
        snow(k, j) = rain(k, j);
        rain(k, j) = 0.0;
      } else {
        snow(k, j) = 0.0;
      }
    }
  }
}

//
// Performance metric Clock
//
namespace Clock {
using time_point_type = std::chrono::high_resolution_clock::time_point;
using duration_type = std::chrono::duration<double>;

time_point_type time() { return std::chrono::high_resolution_clock::now(); }

#ifdef HAVE_MPI

std::array<double, 3> min_max_mean(const MPI_Comm &comm, duration_type duration) {
  double duration_d(duration.count());
  double min, max, mean;
  MPI_Reduce(&duration_d, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&duration_d, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&duration_d, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  mean /= numprocs;

  return std::array<double, 3>{min, max, mean};
}

#endif

} // namespace Clock

} // namespace Utils
} // namespace ELM

#endif
