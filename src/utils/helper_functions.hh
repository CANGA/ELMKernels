

#pragma once

#include "compile_options.hh"

//namespace ELM {
//class AtmDataManager;
//}

namespace ELM::Utils {
  template <class Array_t> Array_t create(const std::string &name, int D0)
  { return Array_t(name, D0); }
  template <class Array_t> Array_t create(const std::string &name, int D0, int D1)
  { return Array_t(name, D0, D1); }
  template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2)
  { return Array_t(name, D0, D1, D2); }
  template <class Array_t, typename T> void assign(Array_t &arr, T&& val)
  { NS::deep_copy(arr, val); }

// template <AtmForcType ftype>
// atm_forc_util<ftype> create_forc_util(const std::string& filename,
//                                       const ELM::Utils::Date &file_start_time,
//                                       const int ntimes, const int ncells)
// { return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }
} // namespace ELM::Utils
