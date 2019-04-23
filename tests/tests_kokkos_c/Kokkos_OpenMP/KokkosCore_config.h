/* ---------------------------------------------
Makefile constructed configuration:
Thu Apr 18 13:30:50 EDT 2019
----------------------------------------------*/
#if !defined(KOKKOS_MACROS_HPP) || defined(KOKKOS_CORE_CONFIG_H)
#error "Do not include KokkosCore_config.h directly; include Kokkos_Macros.hpp instead."
#else
#define KOKKOS_CORE_CONFIG_H
#endif
/* Execution Spaces */
#define KOKKOS_ENABLE_OPENMP
#ifndef __CUDA_ARCH__
#define KOKKOS_USE_ISA_X86_64
#endif
/* General Settings */
#define KOKKOS_ENABLE_CXX11
#define KOKKOS_ENABLE_PROFILING
#define KOKKOS_ENABLE_DEPRECATED_CODE
/* Optimization Settings */
/* Cuda Settings */
#define KOKKOS_ARCH_AVX2
