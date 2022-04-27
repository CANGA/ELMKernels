
#pragma once

#ifndef ENABLE_KOKKOS
#define ACCELERATED
constexpr bool KOKKOS_ENABLED = false;
#else
#define ACCELERATED KOKKOS_INLINE_FUNCTION
constexpr bool KOKKOS_ENABLED = true;
#include "Kokkos_Core.hpp"
#endif

#include "kokkos_types.hh"
