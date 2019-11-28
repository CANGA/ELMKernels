//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>

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
    




} // namespace Utils
} // namespace ELM


#endif
