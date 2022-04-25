
#pragma once

#ifndef ENABLE_KOKKOS
#define ACCELERATED
constexpr bool KOKKOS_ENABLED = false;
#else
#define ACCELERATED KOKKOS_INLINE_FUNCTION
constexpr bool KOKKOS_ENABLED = true;
#include "Kokkos_Core.hpp"
#endif


//?
//template<class object>
//void invoke_kernel() {
//  if constexpr (KOKKOS_ENABLED) {
//    kokkos::parallel_for();
//  } else {
//    std::invoke();
//  }
//}
