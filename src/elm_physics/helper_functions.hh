

#pragma once

#include "kokkos_includes.hh"

namespace ELM::Utils {
  template <class Array_t> Array_t create(const std::string &name, int D0)
  { return Array_t(name, D0); }
  template <class Array_t> Array_t create(const std::string &name, int D0, int D1)
  { return Array_t(name, D0, D1); }
  template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2)
  { return Array_t(name, D0, D1, D2); }
  template <class Array_t, typename T> void assign(Array_t &arr, T val)
  { NS::deep_copy(arr, val); }
} // namespace ELM::Utils
