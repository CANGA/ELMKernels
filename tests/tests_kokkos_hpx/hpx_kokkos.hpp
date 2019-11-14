#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <hpx/include/async.hpp>
#include <hpx/parallel/executors/service_executors.hpp>
#include <iostream>
#include <type_traits>

namespace hpx {
namespace Kokkos {

template <typename... Args> inline void launch(Args &&... args) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  hpx::parallel::execution::sync_execute(main_thread_exec,
                                         std::forward<Args>(args)...);
}

template <typename... Args> inline hpx::future<void> launch_async(Args &&... args) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  return hpx::parallel::execution::async_execute(main_thread_exec,
                                                 std::forward<Args>(args)...);
}

template <typename ExecutionPolicy, typename F>
inline void parallel_for(ExecutionPolicy &&e, F &&f) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  hpx::parallel::execution::sync_execute(main_thread_exec, [&e, &f]() {
    ::Kokkos::parallel_for(std::forward<ExecutionPolicy>(e),
                           std::forward<F>(f));
  });
}

template <typename ExecutionPolicy, typename F>
inline hpx::future<void> parallel_for_async(ExecutionPolicy &&e, F &&f) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  return hpx::parallel::execution::async_execute(main_thread_exec, [&e, &f]() {
    ::Kokkos::parallel_for(std::forward<ExecutionPolicy>(e),
                           std::forward<F>(f));
  });
}

template <typename ExecutionPolicy, typename F, typename Init>
inline void parallel_reduce(ExecutionPolicy &&e, F &&f, Init &&i) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  hpx::parallel::execution::sync_execute(main_thread_exec, [&e, &f, &i]() {
    ::Kokkos::parallel_reduce(std::forward<ExecutionPolicy>(e),
                              std::forward<F>(f), std::forward<Init>(i));
  });
}

template <typename ExecutionPolicy, typename F, typename Init>
inline hpx::future<void> parallel_reduce_async(ExecutionPolicy &&e, F &&f,
                                               Init &&i) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  return hpx::parallel::execution::async_execute(
      main_thread_exec, [&e, &f, &i]() {
        ::Kokkos::parallel_reduce(std::forward<ExecutionPolicy>(e),
                                  std::forward<F>(f), std::forward<Init>(i));
      });
}

template <typename ExecutionPolicy, typename F>
inline void parallel_scan(ExecutionPolicy &&e, F &&f) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  hpx::parallel::execution::sync_execute(main_thread_exec, [&e, &f]() {
    ::Kokkos::parallel_scan(std::forward<ExecutionPolicy>(e),
                            std::forward<F>(f));
  });
}

template <typename ExecutionPolicy, typename F>
inline hpx::future<void> parallel_scan(ExecutionPolicy &&e, F &&f) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  return hpx::parallel::execution::async_execute(main_thread_exec, [&e, &f]() {
    ::Kokkos::parallel_scan(std::forward<ExecutionPolicy>(e),
                            std::forward<F>(f));
  });
}

template <typename... Args>
inline ::Kokkos::View<Args...> make_view(std::string name, std::size_t count) {
  hpx::parallel::execution::service_executor main_thread_exec(
      hpx::threads::executors::service_executor_type::main_thread);

  return hpx::parallel::execution::sync_execute(
      main_thread_exec,
      [&name, count]() { return ::Kokkos::View<Args...>(name, count); });
}

using ::Kokkos::create_mirror_view;
using ::Kokkos::deep_copy;

} // namespace Kokkos
} // namespace hpx
