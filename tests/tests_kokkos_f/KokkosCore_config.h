/* ---------------------------------------------
Makefile constructed configuration:
Wed Mar 27 16:19:36 EDT 2019
----------------------------------------------*/
#if !defined(KOKKOS_MACROS_HPP) || defined(KOKKOS_CORE_CONFIG_H)
#error "Do not include KokkosCore_config.h directly; include Kokkos_Macros.hpp instead."
#else
#define KOKKOS_CORE_CONFIG_H
#endif
/* Execution Spaces */
#define KOKKOS_ENABLE_CUDA
#define KOKKOS_COMPILER_CUDA_VERSION 91
#define KOKKOS_ENABLE_SERIAL
#ifndef __CUDA_ARCH__
#define KOKKOS_USE_ISA_X86_64
#endif
/* General Settings */
#define KOKKOS_ENABLE_CXX11
#define KOKKOS_ENABLE_PROFILING
#define KOKKOS_ENABLE_DEPRECATED_CODE
/* Optimization Settings */
/* Cuda Settings */
#define KOKKOS_ENABLE_CUDA_UVM
#define KOKKOS_ENABLE_CUDA_LAMBDA
#define KOKKOS_ARCH_AVX2
#define KOKKOS_ARCH_PASCAL
#define KOKKOS_ARCH_PASCAL60