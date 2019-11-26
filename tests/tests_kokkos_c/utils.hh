//! A set of utilities for testing ELM kernels in C++

#ifndef UTILS_HH_
#define UTILS_HH_

#include "mpi.h"
#include "Kokkos_Core.hpp"

template<typename T, template View_type>
std::array<3,T> min_max_sum1(const Teuchos::Comm<>& comm,
                             const View_type& v)
{
  View_type::const_view_type vc = v;

  std::array<3,T> results;
  {
    double result;
    Kokkos::parallel_reduce(
        "min", vc.extent(0),
        KOKKOS_LAMBDA(const int& i, T& min_val) {
          min_val = min(min_val, vc(i));
        }, result);
    
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &result, &results[0]);
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "max", vc.extent(0),
        KOKKOS_LAMBDA(const int& i, T& max_val) {
          max_val = max(max_val, vc(i));
        }, result);
    
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &result, &results[1]);
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "sum", vc.extent(0),
        KOKKOS_LAMBDA(const int& i, T& sum_val) {
          sum_val += vc(i);
        }, result);
    
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &result, &results[2]);
  }             
  return results;  
}


template<typename T, template View_type>
std::array<3,T> min_max_sum2(const View_type& v)
{
  View_type::const_view_type vc = v;
  Kokkos::MDRangePolicy<Kokkos::Range<2>> range({0,0},{vc.extent(0), vc.extent(1)});

  std::array<3,T> results;
  {
    double result;
    Kokkos::parallel_reduce(
        "min", range,
        KOKKOS_LAMBDA(const int& i, const int& j, T& min_val) {
          min_val = min(min_val, vc(i,j));
        }, result);
    
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &result, &results[0]);
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "max", range,
        KOKKOS_LAMBDA(const int& i, const int& j, T& max_val) {
          max_val = max(max_val, vc(i,j));
        }, result);
    
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &result, &results[1]);
  }             
  {
    double result;
    Kokkos::parallel_reduce(
        "sum", range,
        KOKKOS_LAMBDA(const int& i, const int& j, T& sum_val) {
          sum_val += vc(i,j);
        }, result);
    
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &result, &results[2]);
  }             
  return results;  
}





#endif
