#include <iostream>
#include "utils.hh"

namespace ELM {
namespace Utils {


std::array<int,2>
square_numprocs(int nprocs) {
  int numprocs_x=1;
  int numprocs_y=1;
  int n=nprocs;
  int i=2;
  while(n!=1){
    if(n%i==0) {
      if(numprocs_x-numprocs_y>=0) {
        numprocs_y*=i;
      } else {
        numprocs_x*=i;
      }
      n=n/i;
    } else {
      i++;
    }
  }
  return { numprocs_x, numprocs_y };
}


DomainDecomposition<1>
create_domain_decomposition_1D(int nprocs, GO n_global, int proc_index)
{
  DomainDecomposition<1> d;
  d.n_procs[0] = nprocs;
  d.proc_index[0] = proc_index;
  d.n_global[0] = n_global;

  int n_local_small = n_global / nprocs;
  int n_local_big = n_local_small + 1;
  int n_procs_big = n_global % nprocs;

  d.n_local[0] = proc_index < n_procs_big ? n_local_big : n_local_small;
  if (proc_index < n_procs_big) {
    d.start[0] = (size_t) n_local_big * proc_index;
  } else {
    d.start[0] = (size_t) n_local_big * n_procs_big + n_local_small * (proc_index - n_procs_big);
  }
  return d;
}


DomainDecomposition<2>
create_domain_decomposition_2D(std::array<int,2> n_procs,
        std::array<GO,2> n_global,
        std::array<int,2> proc_index)
{
  // std::cout << "Creating with: n_procs = " << n_procs[0] << "," << n_procs[1] << std::endl
  //           << "              n_global = " << n_global[0] << "," << n_global[1] << std::endl
  //           << "              proc_ind = " << proc_index[0] << "," << proc_index[1] << std::endl;
    
  auto dd0 = create_domain_decomposition_1D(n_procs[0], n_global[0], proc_index[0]);
  auto dd1 = create_domain_decomposition_1D(n_procs[1], n_global[1], proc_index[1]);

  DomainDecomposition<2> d;
  d.n_procs[0] = dd0.n_procs[0];
  d.proc_index[0] = dd0.proc_index[0];
  d.n_global[0] = dd0.n_global[0];
  d.start[0] = dd0.start[0];
  d.n_local[0] = dd0.n_local[0];

  d.n_procs[1] = dd1.n_procs[0];
  d.proc_index[1] = dd1.proc_index[0];
  d.n_global[1] = dd1.n_global[0];
  d.start[1] = dd1.start[0];
  d.n_local[1] = dd1.n_local[0];
  return d;
}

#ifdef HAVE_MPI

namespace Clock {

std::array<double,3> min_max_mean(const MPI_Comm& comm, duration_type duration) {
  double duration_d(duration.count());
  double min, max, mean;
  MPI_Reduce(&duration_d, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&duration_d, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&duration_d, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  mean /= numprocs;

  return std::array<double,3>{min, max, mean};
}


} // namespace Clock

#endif


} // namespace Utils
} // namespace ELM
