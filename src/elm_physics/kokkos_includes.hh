
#pragma once

#ifndef ENABLE_KOKKOS
#define ACCELERATE
#define ELM_LAMBDA [=]
constexpr bool KOKKOS_ENABLED = false;
namespace NS = ELM;
#else
#define ACCELERATE KOKKOS_INLINE_FUNCTION
#define ELM_LAMBDA KOKKOS_LAMBDA
constexpr bool KOKKOS_ENABLED = true;
#include "Kokkos_Core.hpp"
namespace NS = Kokkos;
#endif

#include "kokkos_types.hh"
