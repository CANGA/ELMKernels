#ifndef ELM_NETCDF_HH_
#define ELM_NETCDF_HH_

//
// Generic readers/writers using sequential NETCDF
//
// NOTE: systematically, an "const Comm_type& comm" argument appears in interfaces to make
// life easier for client codes.  This can be set to anything and is NEVER used
// in this file.

#include <array>
#include <iostream>
#include <string>

#include "array.hh"
#include "mpi_types.hh"
#include "netcdf.h"
#include "utils.hh"

#define NC_HANDLE_ERROR(status, what)                                                                                  \
  do {                                                                                                                 \
    if (status) {                                                                                                      \
      std::cout << __FILE__ << ':' << __LINE__ << ':' << what << " failed with rc = " << status << ':'                 \
                << nc_strerror(status) << '\n';                                                                        \
      abort();                                                                                                         \
    }                                                                                                                  \
  } while (0)

namespace ELM {
namespace IO {

//
// Generic readers/writers
// -----------------------------------------------------------------------------
inline void error(int status, const std::string &func, const std::string &file, const std::string &var = "") {
  std::string what = func + " \"" + file + ":" + var + "\"";
  NC_HANDLE_ERROR(status, what);
}

//
// Dimensions as they are in the file.
//
template <size_t D>
inline std::array<GO, D> get_dimensions(const Comm_type &comm, const std::string &filename,
                                        const std::string &varname) {
  int nc_id = -1;
  auto status = nc_open(filename.c_str(), NC_NOWRITE, &nc_id);
  error(status, "nc_open", filename);

  int var_id = -1;
  status = nc_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "nc_inq_varid", filename, varname);

  int n_dims;
  std::array<int, NC_MAX_VAR_DIMS> dim_ids;
  status = nc_inq_var(nc_id, var_id, nullptr, nullptr, &n_dims, dim_ids.data(), nullptr);
  error(status, "nc_inq_var", filename, varname);

  assert(n_dims == D);
  std::array<GO, D> dims;
  for (int i = 0; i != n_dims; ++i) {
    GO len_dim;
    status = nc_inq_dimlen(nc_id, dim_ids[i], &len_dim);
    error(status, "nc_inq_dimlen", filename, varname);
    dims[i] = len_dim;
  }

  status = nc_close(nc_id);
  error(status, "nc_close", filename);
  return dims;
}

//
// Read all of the dataset into all of the array.
//
template <size_t D>
inline void read(const Comm_type &comm, const std::string &filename, const std::string &varname,
                 Array<double, D> &arr) {
  std::array<size_t, D> start{};
  std::array<size_t, D> count;
  for (int i = 0; i != D; ++i)
    count[i] = (size_t)arr.extent(i);
  read(filename, varname, start, count, arr.data());
  return;
}

//
// Read some of the double dataset into some of the double array.
//
template <size_t D>
inline void read(const Comm_type &comm, const std::string &filename, const std::string &varname,
                 const std::array<size_t, D> &start, const std::array<size_t, D> &count, double *arr) {
  int nc_id = -1;
  auto status = nc_open(filename.c_str(), NC_NOWRITE, &nc_id);
  error(status, "nc_open", filename);

  int var_id = -1;
  status = nc_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "nc_inq_varid", filename, varname);

  status = nc_get_vara_double(nc_id, var_id, start.data(), count.data(), arr);
  error(status, "nc_get_vara_double", filename, varname);

  status = nc_close(nc_id);
  error(status, "nc_close", filename);
}

//
// Read some of the integer dataset into some of the int array.
//
template <size_t D>
inline void read(const Comm_type &comm, const std::string &filename, const std::string &varname,
                 const std::array<size_t, D> &start, const std::array<size_t, D> &count, int *arr) {
  int nc_id = -1;
  auto status = nc_open(filename.c_str(), NC_NOWRITE, &nc_id);
  error(status, "nc_open", filename);

  int var_id = -1;
  status = nc_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "nc_inq_varid", filename, varname);

  status = nc_get_vara_int(nc_id, var_id, start.data(), count.data(), arr);
  error(status, "nc_get_vara_int", filename, varname);

  status = nc_close(nc_id);
  error(status, "nc_close", filename);
}

//
// Read some of the string variable into some of the char array.
//
template <size_t D>
inline void read(const Comm_type &comm, const std::string &filename, const std::string &varname,
                           const std::array<size_t, D> &start, const std::array<size_t, D> &count, char data[]) {
  int nc_id = -1;
  auto status = nc_open(filename.c_str(), NC_NOWRITE, &nc_id);
  error(status, "nc_open", filename);

  int var_id = -1;
  status = nc_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "nc_inq_varid", filename, varname);

  status = nc_get_vara_text(nc_id, var_id, start.data(), count.data(), data);
  error(status, "nc_get_vara_text", filename, varname);

  status = nc_close(nc_id);
  error(status, "nc_close", filename);
}

//
// open for writing
//
inline void init_writing(const std::string &filename, const std::string &varname,
                         const Utils::DomainDecomposition<2> &dd) {
  // NOTE: this can only happen on one rank!
  int nc_id = -1;
  auto status = nc_create(filename.c_str(), NC_WRITE, &nc_id);
  error(status, "nc_create", filename);

  std::array<int, 2> dim_ids;
  status = nc_def_dim(nc_id, "lat", dd.n_global[0], &dim_ids[0]);
  error(status, "nc_def_dim", filename, "lat");
  status = nc_def_dim(nc_id, "lon", dd.n_global[1], &dim_ids[1]);
  error(status, "nc_def_dim", filename, "lon");

  int var_id = -1;
  status = nc_def_var(nc_id, varname.c_str(), NC_DOUBLE, 2, dim_ids.data(), &var_id);
  error(status, "nc_def_var_id", filename, varname);

  status = nc_enddef(nc_id);
  error(status, "nc_enddef", filename);

  status = nc_close(nc_id);
  error(status, "nc_close", filename);
}

//
// Write all
//
template <size_t D>
inline void write(const std::string &filename, const std::string &varname, const Utils::DomainDecomposition<2> &dd,
                  const Array<double, D> &arr) {
  int nc_id = -1;
  auto status = nc_open(filename.c_str(), NC_WRITE, &nc_id);
  error(status, "nc_open", filename);

  int var_id = -1;
  status = nc_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "nc_inq_varid", filename, varname);

  status = nc_put_vara_double(nc_id, var_id, dd.start.data(), dd.n_local.data(), (double *)arr.data());
  error(status, "nc_put_vara_double", filename, varname);
}

} // namespace IO
} // namespace ELM

#endif
