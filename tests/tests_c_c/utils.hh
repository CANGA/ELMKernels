//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>

namespace ELM {
namespace Utils {

const static std::array<12, std::size_t> days_per_month = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

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


} // namespace Utils
} // namespace ELM


#endif
