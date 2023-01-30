
#pragma once

#include <functional>
#include "compile_options.hh"

// generic interface to parallel dispatch
// decouples programming model library syntax from core physics model
// currently supports Kokkos SIMD and serial C++ test code
// Kokkos parallel_reduce and parallel_scan planned for future

// variadic template interface adds flexibility
// functions should compile to zero-overhead invocations of incoming stateful lambda

namespace ELM {

  // implementation
  namespace impl {
#ifdef ENABLE_KOKKOS // Kokkos implementation

    // simple call to Kokkos::parallel_for
    // args can be (name, policy) or just (policy)
    // policy can be any valid Kokkos execution policy
    template <typename F, typename... Args>
    constexpr decltype(auto) apply_parallel_for_impl(F&& kernel, Args&&... args)
    {
      return Kokkos::parallel_for(std::forward<Args>(args)..., std::forward<F>(kernel));
    }

    // deprecated
    // call Kokkos::parallel_for using tuple for args
    template <typename F, typename T, std::size_t... I>
    constexpr decltype(auto) apply_parallel_for_tuple_impl(F&& kernel, T&& args_tuple, std::index_sequence<I...>)
    {
      return Kokkos::parallel_for(std::get<I>(std::forward<T>(args_tuple))..., std::forward<F>(kernel));
    }

#else

    // interface for serial test implementation
    template <typename F>
    constexpr decltype(auto) apply_parallel_for_impl(F&& kernel, int N)
    {
      return [] (F&& obj, int n) {
        for (int i = 0; i < n; ++i) std::invoke(obj, i);
      }(std::forward<F>(kernel), N);
    }

#endif
  } // namespace impl


  // generic interfaces

  // interface for SIMD kernel dispatch
  template <typename F, typename... Args>
  constexpr decltype(auto) apply_parallel_for(F&& kernel, Args&&... args)
  {
    return impl::apply_parallel_for_impl(
           std::forward<F>(kernel),
           std::forward<Args>(args)...);
  }

  // deprecated
  // interface for SIMD kernel dispatch
  // using tuple for args
  template <typename F, typename T>
  constexpr decltype(auto) apply_parallel_for_tuple(F&& kernel, T&& args_tuple)
  {
    return impl::apply_parallel_for_tuple_impl(
           std::forward<F>(kernel),
           std::forward<T>(args_tuple),
           std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<T>>>{});
  }

} // namespace ELM
