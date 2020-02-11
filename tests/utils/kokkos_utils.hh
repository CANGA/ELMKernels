//! A set of utilities for testing ELM kernels in C++

#ifndef UTILS_HH_
#define UTILS_HH_

#include <cmath>
#include "mpi.h"
#include "Kokkos_Core.hpp"


namespace ELM {
namespace ELMKokkos {

template<class View_type>
std::array<double, 3>
min_max_sum1(const MPI_Comm& comm,
             const View_type& v)
{
  using T = typename View_type::value_type;
  typename View_type::const_type vc = v;

  std::array<double,3> results;
  {
    double result;
    Kokkos::parallel_reduce(
        "min", vc.extent(0),
        KOKKOS_LAMBDA(const int& i, double& min_val) {
          if (vc(i) < min_val) {
            min_val = vc(i);
          }
        }, Kokkos::Min<double>(result));
    
    double result_g;
    MPI_Reduce(&result, &result_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    results[0] = result_g;
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "max", vc.extent(0),
        KOKKOS_LAMBDA(const int& i, double& max_val) {
          if (vc(i) > max_val) {
            max_val = vc(i);
          }
        }, Kokkos::Max<double>(result));
    
    double result_g;
    MPI_Reduce(&result, &result_g, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    results[1] = result_g;
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "sum", vc.extent(0),
        KOKKOS_LAMBDA(const int& i, double& sum_val) {
          sum_val += vc(i);
        }, result);
    
    double result_g;
    MPI_Reduce(&result, &result_g, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    results[2] = result_g;
  }             
  return results;  
}


template<class View_type>
std::array<double,3>
min_max_sum2(const MPI_Comm& comm, const View_type& v)
{
  using T = typename View_type::value_type;
  typename View_type::const_type vc = v;
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> range({0,0},{vc.extent(0), vc.extent(1)});

  std::array<double,3> results;
  {
    double result;
    Kokkos::parallel_reduce(
        "min", range,
        KOKKOS_LAMBDA(const int& i, const int& j, double& min_val) {
          if (vc(i,j) < min_val) min_val = vc(i,j);
        }, Kokkos::Min<double>(result));
    
    double result_g;
    MPI_Reduce(&result, &result_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    results[0] = result_g;
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "max", range,
        KOKKOS_LAMBDA(const int& i, const int& j, double& max_val) {
          if (vc(i,j) > max_val) max_val = vc(i,j);
        }, Kokkos::Max<double>(result));
    
    double result_g;
    MPI_Reduce(&result, &result_g, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    results[1] = result_g;
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "sum", range,
        KOKKOS_LAMBDA(const int& i, const int& j, double& sum_val) {
          sum_val += vc(i,j);
        }, result);
    
    double result_g;
    MPI_Reduce(&result, &result_g, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    results[2] = result_g;
  }             
  return results;  
}


} // namespace Kokkos
} // namespace ELM


#endif
