

#pragma once

#include "elm_constants.h"

#include "kokkos_includes.hh"

namespace ELM::solver {

  template <typename ArrayI1, typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void PDMA(const int& c,
            const ArrayI1 snl,
            const ArrayD3 LHS,
            const ArrayD2 A,
            const ArrayD2 B,
            const ArrayD2 Z,
            ArrayD2 RHS);

} // namespace ELM::solver

#include "pentadiagonal_solver_impl.hh"
