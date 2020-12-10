/*
Description:
This does some sample physics

Inputs:
 in1
 in2
 in3

Outputs:
 out1
 out2

*/

#pragma once

#include "clm_constants.hh"

template<
  class ArrayDouble1D,    // a 1D array of non-const doubles, equivalent to double*
  class cArrayDouble1D,   // a 1D array of const doubles, equivalent to double const *
  class cArrayInt1D>      // a 1D array of const ints, equivalent to int const *
KOKKOS_INLINE_FUNCTION    // NATURE!
void Example(
  const cArrayDouble1D in1,
  const double& in2,
  const cArrayInt1D  in3,
  double& out1,
  const ArrayDouble1D out2) // note, this is equivalent to double * const.  We can
                            // change the contained values, but not the container
{
  // ... do stuff
  out2[0] = 1.0;
  out1 = in3[2] * in1[2];
  return;
}

//
// USE CASE 1
//
// simple driver: a single grid cell
int main1() {
  using ArrayDouble1D = Kokkos::View<double*>;
  using cArrayDouble1D = Kokkos::View<const double*>;
  using ArrayInt1D = Kokkos::View<int*>;

  // construct and initialize
  ArrayDouble1D in1(3); assign(in1, 1.0);
  ArrayDouble1D out2(3); assign(out2, 0.0);
  ArrayInt1D in3(3); assign(in3, 1);
  double out1;

  // call function
  Example((cArrayDouble1D) in1, 3.0, (cArrayInt1D) in3, out1, out2);
}


//
// USE CASE 2
//
// wrapped in a class for easier calling on multiple grid cells
//
template<
  class ArrayDouble2D,
  class ArrayInt2D,
  class ArrayDouble1D>
struct ExampleCaller {
  using cArrayDouble2D = ArrayDouble2D::const_view_type;
  using cArrayDouble1D = ArrayDouble1D::const_view_type;
  using cArrayInt2D = ArrayInt2D::const_view_type;

  ExampleCaller(cArrayDouble2D in1_,
                cArrayDouble1D in2_,
                cArrayInt2D in3_,
                ArrayDouble1D out1_,
                ArrayDouble2D out2_)
    : in1(in1_),
      in2(in2_),
      in3(in3_),
      out1(out1_),
      out2(out2_) {}

  void Compute() {
    Kokkos::parallel_for(in1.extent(0),
                         KOKKOS_LAMBDA(const int& g) {
                           Example(Kokkos::subview(in1, i, Kokkos::ALL),
                                   in2(i),
                                   Kokkos::subview(in3, i, Kokkos:ALL),
                                   out1(i),
                                   Kokkos::subview(out2, i, Kokkos::ALL));
                         });
  }

  cArrayDouble2D in1;
  cArrayDouble1D in2;
  cArrayInt2D in3;
  ArrayDouble1D out1;
  ArrayDouble2D out2;
};


int main2() {
  using ArrayDouble2D = Kokkos::View<double**>;
  using ArrayDouble1D = Kokkos::View<double*>;
  using ArrayInt2D = Kokkos::View<int**>;

  int n_grid_cells = 10;

  // construct and initialize
  ArrayDouble2D in1(n_grid_cells, 3); assign(in1, 1.0);
  ArrayDouble2D out2(n_grid_cells, 3); assign(out2, 0.0);
  ArrayInt2D in3(n_grid_cells, 3); assign(in3, 1);

  ArrayDouble1D in2(n_grid_cells); assign(in2, 11.0);
  ArrayDouble1D out1(n_grid_cells); assign(in2, 11.0);

  // call function
  ExampleCaller ex(in1, in2, in3, out1, out2);
  ex.Compute();
}
