// pseudocode prototype of proposed method to organize related functions into classes

// this class holds related functions derived from a single CLM source file
class RelatedPhysicsClass {

private:
  // variables shared between functions in this class
  // these variables are not needed by any code other than the functions in this class
  int p;
  double q[nlevgrnd];

public:
  template <class Array2d, class Array1d> void method1(Array2d x, Array2d y) {
    // do work
  }

  template <class Array2d, class Array1d> void method2(Array2d y, Array1d z) {
    // do work
  }

  template <class Array2d, class Array1d> void method3(Array1d z, Array2d x) {
    // do work
  }
};

// this struct provides a wrapper for calling across multiple grid cells
// it also provides separation between the physics library and Kokkos implementation
template <class Array2d, class Array1d> struct CallRelatedPhysics {

  // variable declarations go here

  // constructor - initialize incoming data to local variables
  // also pass in 1D Kokkos view of RelatedPhysicsClass
  CallRelatedPhysics(Array2d x_, Array2d y_, Array1d z_, Array1d RelatedPhys_)
      : x(x_) y(y_) z(z_) RelatedPhys(RelatedPhys_) {}

  void ComputePhysics() {
    Kokkos::parallel_for(
        num_cells, KOKKOS_LAMBDA(const int &i) {
          RelatedPhys(i).method1(Kokkos::subview(x, i, Kokkos::ALL), Kokkos::subview(y, i, Kokkos::ALL));
          RelatedPhys(i).method2(Kokkos::subview(y, i, Kokkos::ALL), z(i));
          RelatedPhys(i).method3(z(i), Kokkos::subview(x, i, Kokkos::ALL));
        });
  }
};

int main() {
  // set number of grid cells
  int num_cells = 10;
  // assign input and output variables as appropriate Kokkos view template types
  Array2d x;
  Array2d y;
  Array1d z;

  // create 1D Kokkos view from RelatedPhysicsClass
  // not sure if this is necessary?
  // seems like it might be, because data is now encapsulated in RelatedPhysicsClass
  // would use Array1d IRL, but view syntax used here for clarity
  Kokkos::View<RelatedPhysicsClass *> RelatedPhys("set of related physics functions", num_cells);

  // instantiate call wrapper
  CallRelatedPhysics caller(x, y, z, RelatedPhys);
  // call set of related physics functions
  caller.ComputePhysics();

  return 0;
}
