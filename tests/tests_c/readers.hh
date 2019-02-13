//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_NETCDF_HH_
#define ELM_KERNEL_TEST_NETCDF_HH_

#include "netcdf.h"

#define NC_HANDLE_ERROR( status, what )         \
do {                                         \
   if ( status )                             \
   {                                         \
      std::cout                              \
         << __FILE__                         \
         << ':'                              \
         << __LINE__                         \
         << ':'                              \
         << what                             \
         << " failed with rc = "             \
         << status                           \
         << ':'                              \
         << nc_strerror( status )            \
         << '\n' ;                           \
      abort() ;                              \
   }                                         \
} while ( 0 )


namespace ELM {
namespace Utils {

//
// Reads n_grid_cells worth of phenology data from NetCDF files, each with
// n_pft PFTs, into LAI and SAI matrices with a potential offset (in grid
// cells) given by offset.
// -----------------------------------------------------------------------------
template<typename Matrix_t>
void read_phenology(const std::string& fname,
                    size_t n_grid_cells, size_t n_pfts,
                    size_t offset,
                    Matrix_t& lai, Matrix_t& sai) {
  int ncid = -1;
  auto status = nc_open(fname.c_str(), NC_NOWRITE, &ncid);
  NC_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname+"\"");

  auto start = std::array<size_t,4>{0,0,0,0};
  auto count = std::array<size_t,4>{n_grid_cells, n_pfts, 1, 1};

  // get id and read LAI
  int varid = -1;
  status = nc_inq_varid(ncid, "MONTHLY_LAI", &varid);
  NC_HANDLE_ERROR(status, "nc_inq_varid for LAI");

  std::vector<double> data(n_grid_cells * n_pfts);
  status = nc_get_vara_double(ncid, varid, start.data(), count.data(), data.data());
  NC_HANDLE_ERROR(status, "nc_get_vara_double monthly_lai");

  // copy into array -- note this may thrash cache depending upon order of
  // LAI.  Good for C order, bad for Fortran
  for (int i=offset; i!=offset+n_grid_cells; ++i) {
    for (int j=0; j!=n_pfts; ++j) {
      lai(i,j) = data[(i-offset)*n_pfts + j];
    }
  }
  
  
  // get id and read SAI directly into the array
  status = nc_inq_varid(ncid, "MONTHLY_SAI", &varid);
  NC_HANDLE_ERROR( status, "nc_inq_varid for SAI" ) ;
  
  status = nc_get_vara_double(ncid, varid, start.data(), count.data(), data.data());
  NC_HANDLE_ERROR(status, "nc_get_vara_double monthly_sai");
  // copy into array -- note this may thrash cache depending upon order of
  // LAI.  Good for C order, bad for Fortran
  for (int i=offset; i!=offset+n_grid_cells; ++i) {
    for (int j=0; j!=n_pfts; ++j) {
      sai(i,j) = data[(i-offset)*n_pfts + j];
    }
  }

  status = nc_close(ncid);
  NC_HANDLE_ERROR( status, "nc_close" ) ;
}


//
// Reads NetCDF met forcing data, returning the min number of times read.
// -----------------------------------------------------------------------------
template<typename Matrix_t>
int read_forcing(const std::string& fname,
                 size_t n_times, size_t start_grid_cell, size_t n_grid_cells,
                 Matrix_t& rain, Matrix_t& snow, Matrix_t& temp) {

  size_t min_ntimes = n_times;
  for (size_t lcv_gc=0; lcv_gc!=n_grid_cells; ++lcv_gc) {
    std::stringstream fname_full;
    fname_full << fname << lcv_gc+1+start_grid_cell << ".nc";
    
    int ncid = -1;
    auto status = nc_open(fname_full.str().c_str(), NC_NOWRITE, &ncid);
    NC_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname_full.str()+"\"");

    int dimid = -1;
    status = nc_inq_dimid(ncid, "time", &dimid);
    NC_HANDLE_ERROR(status, "nc_inq_dimid");

    size_t ntimes_l;
    status = nc_inq_dimlen(ncid, dimid, &ntimes_l);
    NC_HANDLE_ERROR(status, "nc_inq_dimlen");

    // note the number of local times may be less than the max number of times
    assert(ntimes_l <= n_times);
    min_ntimes = std::min(ntimes_l, min_ntimes);
    auto start = std::array<size_t,3>{0,0,0};
    auto count = std::array<size_t,3>{ntimes_l, 1, 1};
    std::vector<double> data_precip(ntimes_l);
    std::vector<double> data_temp(ntimes_l);

    int varid = -1;
    status = nc_inq_varid(ncid,"PRECTmms", &varid);
    NC_HANDLE_ERROR(status, "nc_inq_varid");

    status = nc_get_vara_double(ncid, varid, start.data(), count.data(), data_precip.data());
    NC_HANDLE_ERROR( status, "nc_get_vara_double total_precip" );

    status = nc_inq_varid(ncid,"TBOT", &varid);
    NC_HANDLE_ERROR(status, "nc_inq_varid");

    status = nc_get_vara_double(ncid, varid, start.data(), count.data(), data_temp.data());
    NC_HANDLE_ERROR( status, "nc_get_vara_double temperature" );
    
    status = nc_close(ncid);
    NC_HANDLE_ERROR( status, "nc_close" ) ;

    // allocate the precip to rain or snow
    for (int lcv_t=0; lcv_t!=ntimes_l; ++lcv_t) {
      if (data_temp[lcv_t] < 273.15) {
        snow[lcv_t][lcv_gc] = data_precip[lcv_t];
        rain[lcv_t][lcv_gc] = 0.;
        temp[lcv_t][lcv_gc] = data_temp[lcv_t];
      } else {
        snow[lcv_t][lcv_gc] = 0.;
        rain[lcv_t][lcv_gc] = data_precip[lcv_t];
        temp[lcv_t][lcv_gc] = data_temp[lcv_t];
      }
    }
  }
  return min_ntimes;
}


} // namespace Utils
} // namespace ELM


#endif
