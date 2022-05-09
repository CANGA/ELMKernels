#ifndef ELM_UTILS_READ_INPUT_HH_
#define ELM_UTILS_READ_INPUT_HH_

#ifdef HAVE_LEGION
#include "legion.h"
#endif

#include "mpi_types.hh"
#include <stdio.h>
#include <string.h>

#ifdef HAVE_PNETCDF
#include "read_pnetcdf.hh"
#else
#include "read_netcdf.hh"
#endif

#include "array.hh"
#include "date_time.hh"

namespace ELM {
namespace IO {

//
// Readers for forcing data.
// -----------------------------------------------------------------------------

//
// Returns shape in a forcing file, as { N_TIMES, N_LAT_GLOBAL, N_LON_GLOBAL }
//
std::array<GO, 3> get_forcing_dimensions(const MPI_Comm &comm, const std::string &dir, const std::string &basename,
                                         const std::string &varname, const Utils::Date &time_start, int n_months);

//
// Read a forcing file
//
// Requires shape(arr) == { N_TIMES, N_LAT_LOCAL, N_LON_LOCAL }
//
void read_forcing(const std::string &dir, const std::string &basename, const std::string &varname,
                  const Utils::Date &time_start, int n_months, const Utils::DomainDecomposition<2> &dd,
                  Array<double, 3> &arr);

#ifndef HAVE_LEGION

//
// Read a variable from the forcing files.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL }
//
template <class Array_t>
inline void read_and_reshape_forcing(const std::string &dir, const std::string &basename, const std::string &varname,
                                     const Utils::Date &time_start, int n_months,
                                     const Utils::DomainDecomposition<2> &dd, Array_t &arr);

#else // HAVE_LEGION

//
// Read a variable from the forcing files into a Legion accessor.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL }
//
template <class Array_t>
inline void read_and_reshape_forcing(const std::string &dir, const std::string &basename, const std::string &varname,
                                     const Utils::Date &time_start, int n_months,
                                     const Utils::DomainDecomposition<2> &dd, Array_t &arr);

#endif // HAVE_LEGION

//
// Readers for pft constants data.
// -----------------------------------------------------------------------------

//
// Returns maxpfts in a pft_constants file
//
std::array<GO, 1> get_maxpfts(const Comm_type &comm, const std::string &dir, const std::string &basename,
                              const std::string &varname);

//
// Readers for phenology data.
// -----------------------------------------------------------------------------

//
// Returns shape in a phenology file, as { N_TIMES, N_PFTS, N_LAT_GLOBAL, N_LON_GLOBAL }
//
std::array<GO, 4> get_phenology_dimensions(const Comm_type &comm, const std::string &dir, const std::string &basename,
                                           const std::string &varname, const Utils::Date &time_start, int n_months);

//
// Read a phenology file.
//
// Requires shape(arr) == { N_TIMES, N_PFTS, N_LAT_LOCAL, N_LON_LOCAL }
//
void read_phenology(const std::string &dir, const std::string &basename, const std::string &varname,
                    const Utils::Date &time_start, int n_months, const Utils::DomainDecomposition<2> &dd,
                    Array<double, 4> &arr);

#ifndef HAVE_LEGION

//
// Read a variable from the forcing files.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL, N_PFTS }
//
template <class Array_t>
inline void read_and_reshape_phenology(const std::string &dir, const std::string &basename, const std::string &varname,
                                       const Utils::Date &time_start, int n_months,
                                       const Utils::DomainDecomposition<2> &dd, Array_t &arr);

#else // HAVE_LEGION

//
// Read a variable from the forcing files into a Legion accessor.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL, N_PFTS }
//
template <class Array_t>
inline void read_and_reshape_phenology(const std::string &dir, const std::string &basename, const std::string &varname,
                                       const Utils::Date &time_start, int n_months,
                                       const Utils::DomainDecomposition<2> &dd, int n_pfts, Array_t &arr);

#endif // HAVE_LEGION

//
// Writers by grid cell
// -----------------------------------------------------------------------------

//
// Assumes shape(arr) == { N_GRID_CELLS_LOCAL }
//
template <typename Array_t>
inline void reshape_and_write_grid_cell(const std::string &filename, const std::string &varname,
                                        const Utils::DomainDecomposition<2> &dd, const Array_t &arr);

//
// IMPLEMENTATION
//

//
// Readers for forcing data.
// -----------------------------------------------------------------------------

#ifndef HAVE_LEGION

//
// Read a variable from the forcing files.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL }
//
template <class Array_t>
inline void read_and_reshape_forcing(const std::string &dir, const std::string &basename, const std::string &varname,
                                     const Utils::Date &time_start, int n_months,
                                     const Utils::DomainDecomposition<2> &dd, Array_t &arr) {
  assert(static_cast<size_t>(arr.extent(1)) == dd.n_local[0] * dd.n_local[1]);
  Array<double, 3> arr_for_read(arr.extent(0), dd.n_local[0], dd.n_local[1]);
  read_forcing(dir, basename, varname, time_start, n_months, dd, arr_for_read);

  for (int i = 0; i != arr.extent(0); ++i) {
    for (int j = 0; j != static_cast<int>(dd.n_local[0]); ++j) {
      for (int k = 0; k != static_cast<int>(dd.n_local[1]); ++k) {
        arr(i, j * dd.n_local[1] + k) = arr_for_read(i, j, k);
      }
    }
  }
}

#else // HAVE_LEGION

//
// Read a variable from the forcing files into a Legion accessor.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL }
//
template <class Array_t>
inline void read_and_reshape_forcing(const std::string &dir, const std::string &basename, const std::string &varname,
                                     const Utils::Date &time_start, int n_months, int n_times,
                                     const Utils::DomainDecomposition<2> &dd, Array_t &arr) {
  Array<double, 3> arr_for_read(n_times, dd.n_local[0], dd.n_local[1]);
  read_forcing(dir, basename, varname, time_start, n_months, dd, arr_for_read);

  for (int i = 0; i != n_times; ++i) {
    for (int j = 0; j != dd.n_local[0]; ++j) {
      for (int k = 0; k != dd.n_local[1]; ++k) {
        arr(Legion::Point<2>{i, j * dd.n_local[1] + k}) = arr_for_read(i, j, k);
      }
    }
  }
}

#endif // HAVE_LEGION

#ifndef HAVE_LEGION

//
// Read a variable from the pft constant file.
//
// Requires shape(arr) == { MAXPFTS }
//
template <class Array_t>
inline void read_pft_var(const Comm_type &comm, const std::string& filename,
                         const std::string &varname, Array_t &arr)
{
  Array<double, 1> arr_for_read(arr.extent(0));
  std::array<GO, 1> start = {0};
  std::array<GO, 1> count = {(GO)arr.extent(0)};
  read(comm, filename, varname, start, count, arr_for_read.data());
  for (int i = 0; i != arr.extent(0); ++i) {
    arr(i) = arr_for_read(i);
  }
}

//
// Read entire string variable from the forcing files.
//
// Requires shape(arr) == { MAX_PFTS }
//
template <class Array_t>
inline void read_names(const Comm_type &comm, const std::string& filename,
                       const std::string &varname, const int strlen, Array_t &arr) {
  // ELM::Array of char
  Array<char, 1> arr_for_read(strlen * arr.extent(0) + 1);
  std::array<GO, 2> start = {0, 0};
  std::array<GO, 2> count = {(GO)arr.extent(0), (GO)strlen};
  read(comm, filename, varname, start, count, arr_for_read.data());
  // c string parsing
  char *token = strtok(arr_for_read.data(), " ");
  arr[0] = token;
  int i = 1;
  while (token != NULL && i != arr.extent(0)) {
    token = strtok(NULL, " ");
    arr[i] = token;
    i++;
  }
}


template <typename T, size_t D>
inline void read_netcdf(const Comm_type &comm, const std::string &filename, const std::string &varname,
                 const std::array<GO, D> &start, const std::array<GO, D> &count, T *arr) {
  std::array<GO, D> new_start;
  std::copy(begin(start), end(start), begin(new_start));
  std::array<GO, D> new_count;
  std::copy(begin(count), end(count), begin(new_count));
  read(comm, filename, varname, start, count, arr);
}


template <class Array_t>
inline void read_and_reshape_phenology(const std::string &dir, const std::string &basename, const std::string &varname,
                                       const Utils::Date &time_start, int n_months,
                                       const Utils::DomainDecomposition<2> &dd, Array_t &arr) {
  assert(static_cast<size_t>(arr.extent(1)) == dd.n_local[0] * dd.n_local[1]);
  Array<double, 4> arr_for_read(arr.extent(0), arr.extent(2), dd.n_local[0], dd.n_local[1]);
  read_phenology(dir, basename, varname, time_start, n_months, dd, arr_for_read);
  for (int i = 0; i != arr.extent(0); ++i) {
    for (int p = 0; p != arr.extent(2); ++p) {
      for (int j = 0; j != static_cast<int>(dd.n_local[0]); ++j) {
        for (int k = 0; k != static_cast<int>(dd.n_local[1]); ++k) {
          arr(i, j * dd.n_local[1] + k, p) = arr_for_read(i, p, j, k);
        }
      }
    }
  }
}

#else // HAVE_LEGION

//
// Read a variable from the forcing files into a Legion accessor.
//
// Requires shape(arr) == { N_TIMES, N_GRID_CELLS_LOCAL, N_PFTS }
//
template <class Array_t>
inline void read_and_reshape_phenology(const std::string &dir, const std::string &basename, const std::string &varname,
                                       const Utils::Date &time_start, int n_months,
                                       const Utils::DomainDecomposition<2> &dd, int n_pfts, Array_t &arr) {
  Array<double,4> arr_for_read(n_months, n_pfts], dd.n_local[0], dd.n_local[1]);
  read_phenology(dir, basename, varname, start_year, time_start, n_months, dd, arr_for_read);
  for (int i = 0; i != n_months; ++i) {
    for (int p = 0; p != n_pfts; ++p) {
      for (int j = 0; j != dd.n_local[0]; ++j) {
        for (int k = 0; k != dd.n_local[1]; ++k) {
          arr(Legion::Point<3>{i, j * dd.n_local[1] + k, p}) = arr_for_read(i, p, j, k);
        }
      }
    }
  }
}

#endif // HAVE_LEGION

//
// Writers by grid cell
// -----------------------------------------------------------------------------

//
// Assumes shape(arr) == { N_GRID_CELLS_LOCAL }
//
template <typename Array_t>
inline void reshape_and_write_grid_cell(const std::string &filename, const std::string &varname,
                                        const Utils::DomainDecomposition<2> &dd, const Array_t &arr) {
  Array<double, 2> arr_for_write(dd.n_local[0], dd.n_local[1]);
  for (int i = 0; i != dd.n_local[0]; ++i) {
    for (int j = 0; j != dd.n_local[1]; ++j) {
      arr_for_write(i, j) = arr[i * dd.n_local[1] + j];
    }
  }
  write<2>(filename, varname, dd, arr_for_write);
}

} // namespace IO
} // namespace ELM

#endif
