#ifndef ELM_MPI_TYPES_HH_
#define ELM_MPI_TYPES_HH_

#include <cstddef>

#ifdef HAVE_MPI

#include "mpi.h"
using Comm_type = MPI_Comm;

#ifdef HAVE_PNETCDF
using GO = MPI_Offset;
#else
using GO = size_t;
#endif

#else

using MPI_Comm = int;
using Comm_type = MPI_Comm;
using GO = size_t;

#endif

#endif
