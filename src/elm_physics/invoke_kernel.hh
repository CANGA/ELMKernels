
#pragma once


#include <functional>
#include <tuple>
#include <utility>

#include "kokkos_includes.hh"

// is there a better way than ifdefs to do conditional compilation?
// wrapping in constexpr if still required Kokkos to be available 

#ifdef ENABLE_KOKKOS
namespace impl {
template <class F, class T, std::size_t... I>
constexpr decltype(auto) apply_kokkos_parallel_for(F&& obj, T&& args, const std::string& name, std::index_sequence<I...>)
{
  auto run_kernel = [] (F&& obj, T&& args, const std::string& name) {
    Kokkos::parallel_for(name, std::get<I>(std::forward<T>(args))..., std::forward<F>(obj));
  };

  return run_kernel(std::forward<F>(obj), std::forward<T>(args), name);
}
}  // namespace impl

template <class F, typename T>
constexpr decltype(auto) invoke_kernel(F&& obj, T&& args, const std::string& name = "") {
  return impl::apply_kokkos_parallel_for(
         std::forward<F>(obj), std::forward<T>(args), name,
         std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<T>>>{});
}
#else

template <class F, typename T>
constexpr decltype(auto) invoke_kernel(F&& obj, T&& args, const std::string& name = "") {
  // args will be a scalar int (index/num func calls)

  auto run_kernel = [] (F&& obj, int N) {
    for (int i = 0; i < N; ++i) {
      std::invoke(std::forward<F>(obj), i);
    }
  };

  return run_kernel(std::forward<F>(obj), std::get<0>(args));
}

#endif

