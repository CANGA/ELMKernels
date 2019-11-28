//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>
#include <chrono>

#include "../utils/array.hh"

namespace ELM {
namespace Utils {

//
// Determine month_of_year from day_of_year
//
const static std::array<std::size_t,12> days_per_month = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
//
std::size_t
month_from_day(std::size_t doy) {
  doy = doy % 365;
  std::size_t i_month = 0;
  std::size_t total = 0;
  while (i_month < 12) {
    total += days_per_month[i_month];
    if (doy < total) {
      break;
    } else {
      i_month++;
    }
  }
  assert(i_month < 12 && "Broken month_from_day.");
  return i_month;
}


//
// Domain decomposition in x, y
//
std::pair<std::size_t, std::size_t>
get_domain_decomposition(int nprocs, int argc, char** argv) {
  assert(nprocs == 6 && "NPROCS must be 6 for now!");
  return std::make_pair(3,2);
}


//
// min/max/sum an array
//
template<size_t D>
std::array<double, 3>
min_max_sum(const MPI_Comm& comm, const Array<double,D>& arr)
{
  auto min_max = std::minmax_element(arr.begin(), arr.end());
  double min = *min_max.first;
  double max = *min_max.second;
  double sum = std::accumulate(arr.begin(), arr.end(), 0.);

  double gmin, gmax, gsum;
  MPI_Reduce(&min, &gmin, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  return std::array<double,3>{gmin, gmax, gsum};
}
  

namespace Clock {
using time_point_type = std::chrono::high_resolution_clock::time_point;
using duration_type = std::chrono::duration<double>;

time_point_type time() {
  return std::chrono::high_resolution_clock::now();
}

std::array<double,3> min_max_mean(const MPI_Comm& comm, duration_type duration) {
  double duration_d(duration.count());
  double min, max, mean;
  MPI_Reduce(&duration_d, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&duration_d, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&duration_d, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  mean /= numprocs;

  return std::array<double,3>{min, max, mean};
}

} // namespace Clock




} // namespace Utils
} // namespace ELM


#endif
