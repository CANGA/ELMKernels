
#pragma once

#include "helper_functions.hh"

namespace ELM::solver {

  // basic implementation of a pentadiagonal linear solver
  // will perform well on CPU, possibly poorly on GPU
  // should implement cyclic reduction to maximize GPU performance
  // S. S. Askar, A. A. Karawia, "On Solving Pentadiagonal Linear Systems via Transformations",
  // Mathematical Problems in Engineering, vol. 2015, Article ID 232456, 9 pages, 2015.
  // https://doi.org/10.1155/2015/232456
  template <typename ArrayI1, typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void PDMA(const int& c,
            const ArrayI1 snl,
            const ArrayD3 LHS,
            const ArrayD2 A,
            const ArrayD2 B,
            const ArrayD2 Z,
            ArrayD2 RHS)
  {
    using ELMdims::nlevsno;
    using ELM::Utils::create;

    // maximum size of system
    const int N = RHS.extent(1);
    // top active layer
    const int top = nlevsno - snl(c);

    // form coefficients
    // forward sweep
    // top layer calculations
    double U1 = 1.0 / LHS(c, top, 2);
    A(c, top) = LHS(c, top, 1) * U1;
    B(c, top) = LHS(c, top, 0) * U1;
    Z(c, top) = RHS(c, top) * U1;

    // 2nd from top
    double Y1 = LHS(c, top + 1, 3);
    U1 = 1.0 / (LHS(c, top + 1, 2) - A(c, top) * Y1);
    A(c, top + 1) = (LHS(c, top + 1, 1) - B(c, top) * Y1) * U1;
    B(c, top + 1) = LHS(c, top + 1, 0) * U1;
    Z(c, top + 1) = (RHS(c, top + 1) - Z(c, top) * Y1) * U1;

    // for cells (top + 2 : N - 3)
    for (int i = top + 2; i < N - 2; ++i) {
      Y1 = LHS(c, i, 3) - A(c, i - 2) * LHS(c, i, 4);
      U1 = 1.0 / (LHS(c, i, 2) - B(c, i - 2) * LHS(c, i, 4) - A(c, i - 1) * Y1);
      A(c, i) = (LHS(c, i, 1) - B(c, i - 1) * Y1) * U1;
      B(c, i) = LHS(c, i, 0) * U1;
      Z(c, i) = (RHS(c, i) - Z(c, i - 2) * LHS(c, i, 4) - Z(c, i - 1) * Y1) * U1;
    }

    // 2nd from bottom
    Y1 = LHS(c, N - 2, 3) - A(c, N - 4) * LHS(c, N - 2, 4);
    U1 = 1.0 / (LHS(c, N - 2, 2) - B(c, N - 4) * LHS(c, N - 2, 4) -
      A(c, N - 3) * Y1);
    A(c, N - 2) = (LHS(c, N - 2, 1) - B(c, N - 3) * Y1) * U1;

    // bottom cell
    const double Y2 = LHS(c, N - 1, 3) - A(c, N - 3) * LHS(c, N - 1, 4);
    const double  U2 = 1.0 / (LHS(c, N - 1, 2) - B(c, N - 3) * LHS(c, N - 1, 4) -
      A(c, N - 2) * Y2);
    Z(c, N - 2) = (RHS(c, N - 2) - Z(c, N - 3) * LHS(c, N - 2, 4) -
      Z(c, N - 3) * Y1) * U1;
    Z(c, N - 1) = (RHS(c, N - 1) - Z(c, N - 2) * LHS(c, N - 1, 4) -
      Z(c, N - 2) * Y2) * U2;

    // obtain solution via backward substitution
    RHS(c, N - 1) = Z(c, N - 1);
    RHS(c, N - 2) = Z(c, N - 2) - A(c, N - 2) * RHS(c, N - 1);
    for (int i = N - 3; i >= 0; --i)
      RHS(c, i) = Z(c, i) - A(c, i) * RHS(c, i + 1) - B(c, i) * RHS(c, i + 2);
  }

} // namespace ELM::solver
