
#pragma once

//#ifndef ENABLE_KOKKOS
//
//#define ACCELERATE
//#define ELM_LAMBDA [=]
//namespace NS = ELM;
//
//#else

#define ACCELERATE KOKKOS_INLINE_FUNCTION
#define ELM_LAMBDA KOKKOS_LAMBDA
#include "Kokkos_Core.hpp"
namespace NS = Kokkos;

//#endif
