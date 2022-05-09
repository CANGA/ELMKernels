
#pragma once

#ifndef ENABLE_KOKKOS

namespace ELM { template <typename T, size_t D> class Array; }
typedef ELM::Array<bool, 1> ArrayB1;
typedef ELM::Array<int, 1> ArrayI1;
typedef ELM::Array<int, 2> ArrayI2;
typedef ELM::Array<std::string, 1> ArrayS1;
typedef ELM::Array<double, 1> ArrayD1;
typedef ELM::Array<double, 2> ArrayD2;
typedef ELM::Array<double, 3> ArrayD3;


typedef ArrayD1 h_ArrayD1;
typedef ArrayD2 h_ArrayD2;
typedef ArrayD3 h_ArrayD3;

#else

typedef Kokkos::View<bool *> ArrayB1;
typedef Kokkos::View<int *> ArrayI1;
typedef Kokkos::View<int **> ArrayI2;
typedef Kokkos::View<std::string *> ArrayS1;
typedef Kokkos::View<double *> ArrayD1;
typedef Kokkos::View<double **> ArrayD2;
typedef Kokkos::View<double ***> ArrayD3;


typedef ArrayD1::HostMirror h_ArrayD1;
typedef ArrayD2::HostMirror h_ArrayD2;
typedef ArrayD3::HostMirror h_ArrayD3;

#endif
