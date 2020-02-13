#ifndef ELM_UTILS_READERS_HH_
#define ELM_UTILS_READERS_HH_

#include "array.hh"


// Example of what the reader interface might look like?
namespace ELM {
namespace IO {

//
// Returns shape in file, as a pair (N_TIMES, NX_GLOBAL, NY_GLOBAL)
//
std::tuple<int, int, int>
get_dimensions(const std::string& dir,const std::string& basename,
               int start_year, int start_month, int n_months);


//
// Read a variable from the phenology file.
//
// Assumes shape(arr) == (N_MONTHS, N_PFT, N_LAT_LOCAL, N_LON_LOCAL)
//  
void
read_phenology(const MPI_Comm& comm,
               const std::string& dir,const std::string& basename, const std::string& phenology_type,
               int start_year, int start_month,
               int i_beg, int j_beg, ELM::Utils::Array<double,4>& arr);

//
// Read a variable from the phenology file and copy it into the provided array.
//
// Assumes shape(arr) == (N_MONTHS, N_GRID_CELLS_LOCAL, N_PFT)
//  
template<class Array_t>
void
read_and_reshape_phenology(const MPI_Comm& comm,
			   const std::string& dir,const std::string& basename, const std::string& phenology_type,
			   int start_year, int start_month,
			   int i_beg, int j_beg, int n_lat_local, int n_lon_local, Array_t& arr) 
{
  assert(arr.extent(1) == n_lat_local * n_lon_local);
  ELM::Utils::Array<double,4> arr_for_read(arr.extent(0), arr.extent(2), n_lat_local, n_lon_local);
  read_phenology(comm, dir, basename, phenology_type, start_year, start_month, 
		 i_beg, j_beg, arr_for_read);
  for (int i=0; i!=arr.extent(0); ++i) {
    for (int p=0; p!=arr.extent(2); ++p) {
      for (int j=0; j!=n_lat_local; ++j) {
        for (int k=0; k!=n_lon_local; ++k) {
	  arr(i,j*n_lon_local + k,p) = arr_for_read(i,p,j,k);
	}
      }
    }
  }
}


//
// Read a variable from the forcing files.
//
// Assumes shape(arr) == (N_TIMES, N_LAT_LOCAL, N_LON_LOCAL )
//
void
read_forcing(const MPI_Comm& comm,
             const std::string& dir,const std::string& basename, const std::string& forcing_type,
             int start_year, int start_month, int n_months, 
             int i_beg, int j_beg, ELM::Utils::Array<double,3>& arr);

//
// Read a variable from the forcing files.
//
// Assumes shape(arr) == (N_TIMES, N_GRID_CELLS_LOCAL )
//
template<class Array_t>
void
read_and_reshape_forcing(const MPI_Comm& comm,
	     const std::string& dir,const std::string& basename, const std::string& forcing_type,
             int start_year, int start_month, int n_months,
	     int i_beg, int j_beg, int n_lat_local, int n_lon_local, Array_t& arr) 
{
  assert(arr.extent(1) == n_lat_local * n_lon_local);
  ELM::Utils::Array<double,3> arr_for_read(arr.extent(0), n_lat_local, n_lon_local);
  read_forcing(comm, dir, basename, forcing_type, start_year, start_month, n_months,
	       i_beg, j_beg, arr_for_read);
  for (int i=0; i!=arr.extent(0); ++i) {
    for (int j=0; j!=n_lat_local; ++j) {
      for (int k=0; k!=n_lon_local; ++k) {
	arr(i,j*n_lon_local + k) = arr_for_read(i,j,k);
      }
    }
  }
}



// Convert precipitation to rain and snow
template<class Array_t>
void
convert_precip_to_rain_snow(Array_t& rain,
                            Array_t& snow, 
                            const Array_t& temp)
{
  const int nt = rain.extent(0);
  const int ng = rain.extent(1);

  for (int k=0; k<nt; k++) {
    for (int j=0; j<ng; j++) {
      if (temp(k,j) < 273.15) {
        // no need to update snow, both are initially total precip
        rain(k,j) = 0.0; 
      } else {
        // no need to update rain, both are initially total precip
        snow(k,j) = 0.0;
      }
    }
  }
}


//
// Write a variable to disk
//
// Assumes shape(arr) == (N_LAT_LOCAL, N_LON_LOCAL )
//
void
write_grid_cell(const MPI_Comm& comm,
                const std::string& filename, const std::string& varname,
                int i_beg, int j_beg, int n_lat_global, int n_lon_global,
                ELM::Utils::Array<double,2>& arr);


//
// Write a variable to disk
//
// Assumes shape(arr) == (N_GRID_CELLS )
//
template<class Array_t>
void
reshape_and_write_grid_cell(const MPI_Comm& comm,
                            const std::string& filename, const std::string& varname,
                            int i_beg, int j_beg, int n_lat_local, int n_lon_local,
                            int n_lat_global, int n_lon_global,                            
                            const Array_t& arr)
{
  assert(arr.extent(0) == n_lat_local * n_lon_local);
  ELM::Utils::Array<double,2> arr_for_write(n_lat_local, n_lon_local);
  for (int j=0; j!=n_lat_local; ++j) {
    for (int k=0; k!=n_lon_local; ++k) {
      arr_for_write(j,k) = arr(j*n_lon_local + k);
    }
  }
  write_grid_cell(comm, filename, varname, i_beg, j_beg, n_lat_global, n_lon_global, arr_for_write);
}


} // namespace
} // namespace


#endif
