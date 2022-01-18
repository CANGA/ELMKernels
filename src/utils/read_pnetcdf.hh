#ifndef ELM_PNETCDF_HH_
#define ELM_PNETCDF_HH_

//
// Generic readers/writers using sequential NETCDF
//

#include <array>
#include <iostream>
#include <string>

#include "array.hh"
#include "mpi.h"
#include "pnetcdf.h"
#include "utils.hh"

#define NC_HANDLE_ERROR(status, what)                                                                                  \
  do {                                                                                                                 \
    if (status) {                                                                                                      \
      std::cout << __FILE__ << ':' << __LINE__ << ':' << what << " failed with rc = " << status << ':'                 \
                << ncmpi_strerror(status) << '\n';                                                                     \
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
inline std::array<GO, D> get_dimensions(const MPI_Comm &comm, const std::string &filename, const std::string &varname) {
  MPI_Info info;
  MPI_Info_create(&info);

  int nc_id = -1;
  auto status = ncmpi_open(comm, filename.c_str(), NC_NOWRITE, info, &nc_id);
  error(status, "nc_open", filename);

  int var_id = -1;
  status = ncmpi_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "ncmpi_inq_varid", filename, varname);

  int n_dims;
  std::array<int, 10> dim_ids; // NOTE: This nominally should be
                               // NC_MAX_VAR_DIMS (not 10), but that is NC_MAX_INT,
                               // which breaks this.
  status = ncmpi_inq_var(nc_id, var_id, nullptr, nullptr, &n_dims, dim_ids.data(), nullptr);
  error(status, "ncmpi_inq_var", filename, varname);

  assert(n_dims == D);
  std::array<GO, D> dims;
  for (int i = 0; i != n_dims; ++i) {
    MPI_Offset len_dim;
    status = ncmpi_inq_dimlen(nc_id, dim_ids[i], &len_dim);
    error(status, "ncmpi_inq_dimlen", filename, varname);
    dims[i] = (size_t)len_dim;
  }

  status = ncmpi_close(nc_id);
  error(status, "ncmpi_close", filename);

  MPI_Info_free(&info);
  return dims;
}

//
// get id of dimension variable named by varname
//
inline int get_dimid(const MPI_Comm &comm, const std::string &filename,
                                        const std::string &varname) {
  MPI_Info info;
  MPI_Info_create(&info);
  int nc_id = -1;
  auto status = ncmpi_open(comm, filename.c_str(), NC_NOWRITE, info, &nc_id);
  error(status, "ncmpi_open", filename);
  MPI_Info_free(&info);

  int dim_id = -1;
  status = ncmpi_inq_dimid (nc_id, varname.c_str(), &dim_id);
  error(status, "ncmpi_inq_varid", filename, varname);

  status = ncmpi_close(nc_id);
  error(status, "ncmpi_close", filename);
  return dim_id;
}

//
// return an array of dimension ids corresponding to the dimensions of variable named by varname
//
template <int D>
inline std::array<int, D> get_var_dimids(const MPI_Comm &comm, const std::string &filename,
                                        const std::string &varname) {
  MPI_Info info;
  MPI_Info_create(&info);
  int nc_id = -1;
  auto status = ncmpi_open(comm, filename.c_str(), NC_NOWRITE, info, &nc_id);
  error(status, "ncmpi_open", filename);
  MPI_Info_free(&info);

  int var_id = -1;
  status = ncmpi_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "ncmpi_inq_varid", filename, varname);

  std::array<int, D> dimids{0};
  status = ncmpi_inq_vardimid(nc_id, var_id, dimids.data());
  error(status, "ncmpi_inq_vardimid", filename, varname);

  status = ncmpi_close(nc_id);
  error(status, "ncmpi_close", filename);
  return dimids;
}

//
// Read some.
//
template <size_t D>
inline void read(const MPI_Comm &comm, const std::string &filename, const std::string &varname,
                 const std::array<MPI_Offset, D> &start, const std::array<MPI_Offset, D> &count, double *arr) {
  int nc_id = -1;
  MPI_Info info;
  MPI_Info_create(&info);
  auto status = ncmpi_open(comm, filename.c_str(), NC_NOWRITE, info, &nc_id);
  error(status, "ncmpi_open", filename);
  MPI_Info_free(&info);

  int var_id = -1;
  status = ncmpi_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "ncmpi_inq_varid", filename, varname);

  status = ncmpi_get_vara_double_all(nc_id, var_id, start.data(), count.data(), arr);
  error(status, "ncmpi_get_vara_double", filename, varname);

  status = ncmpi_close(nc_id);
  error(status, "ncmpi_close", filename);
}

//
// open for writing
//
inline void init_writing(const std::string &filename, const std::string &varname,
                         const Utils::DomainDecomposition<2> &dd) {
  // NOTE: this can only happen on one rank!
  int nc_id = -1;
  MPI_Info info;
  MPI_Info_create(&info);
  auto status = ncmpi_create(dd.comm, filename.c_str(), NC_WRITE, info, &nc_id);
  error(status, "ncmpi_create", filename);
  MPI_Info_free(&info);

  std::array<int, 2> dim_ids;
  status = ncmpi_def_dim(nc_id, "lat", dd.n_global[0], &dim_ids[0]);
  error(status, "ncmpi_def_dim", filename, "lat");
  status = ncmpi_def_dim(nc_id, "lon", dd.n_global[1], &dim_ids[1]);
  error(status, "ncmpi_def_dim", filename, "lon");

  int var_id = -1;
  status = ncmpi_def_var(nc_id, varname.c_str(), NC_DOUBLE, 2, dim_ids.data(), &var_id);
  error(status, "ncmpi_def_var_id", filename, varname);

  status = ncmpi_enddef(nc_id);
  error(status, "ncmpi_enddef", filename);

  status = ncmpi_close(nc_id);
  error(status, "ncmpi_close", filename);
}

//
// Write all
//
template <size_t D>
inline void write(const std::string &filename, const std::string &varname, const Utils::DomainDecomposition<D> &dd,
                  const Array<double, D> &arr) {
  MPI_Info info;
  MPI_Info_create(&info);
  int nc_id = -1;
  auto status = ncmpi_open(dd.comm, filename.c_str(), NC_WRITE, info, &nc_id);
  error(status, "ncmpi_open", filename);
  MPI_Info_free(&info);

  int var_id = -1;
  status = ncmpi_inq_varid(nc_id, varname.c_str(), &var_id);
  error(status, "ncmpi_inq_varid", filename, varname);

  status = ncmpi_put_vara_double_all(nc_id, var_id, dd.start.data(), dd.n_local.data(), (double *)arr.data());
  error(status, "ncmpi_put_vara_double", filename, varname);
}

} // namespace IO
} // namespace ELM

#endif
