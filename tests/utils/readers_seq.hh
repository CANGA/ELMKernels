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
read_phenology(const std::string& dir,const std::string& basename, const std::string& phenology_type,
               int start_year, int start_month,
               int i_beg, int j_beg, ELM::Utils::Array<double,4>& arr);

//
// Read a variable from the phenology file and copy it into the provided array.
//
// Assumes shape(arr) == (N_MONTHS, N_GRID_CELLS_LOCAL, N_PFT)
//  
template<class Array_t>
void
read_and_reshape_phenology(const std::string& dir,const std::string& basename, const std::string& phenology_type,int start_year, int start_month,
			   int i_beg, int j_beg, int n_lat_local, int n_lon_local, Array_t& arr) 
{
  std::cout << "n_lat,lon_local = " << n_lat_local << "," << n_lon_local << std::endl;
  std::cout << "arr shape = " << arr.extent(0) << "," << arr.extent(1) << "," << arr.extent(2) << std::endl;
  assert(arr.extent(1) == n_lat_local * n_lon_local);
  ELM::Utils::Array<double,4> arr_for_read(arr.extent(0), arr.extent(2), n_lat_local, n_lon_local);
  read_phenology(dir, basename, phenology_type, start_year, start_month, 
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
read_forcing(const std::string& dir,const std::string& basename, const std::string& forcing_type,
             int start_year, int start_month, int n_months, 
             int i_beg, int j_beg, ELM::Utils::Array<double,3>& arr);

//
// Read a variable from the forcing files.
//
// Assumes shape(arr) == (N_TIMES, N_GRID_CELLS_LOCAL )
//
template<class Array_t>
void
read_and_reshape_forcing(const std::string& dir,const std::string& basename, const std::string& forcing_type, int start_year, int start_month, int n_months,
	     int i_beg, int j_beg, int n_lat_local, int n_lon_local, Array_t& arr) 
{
  assert(arr.extent(1) == n_lon_local * n_lon_local);
  ELM::Utils::Array<double,3> arr_for_read(arr.extent(0), n_lat_local, n_lon_local);
  read_forcing(dir, basename, forcing_type, start_year, start_month, n_months,
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
void convert_precip_to_rain_snow(ELM::Utils::Array<double,3>& rain, ELM::Utils::Array<double,3>& snow, 
											ELM::Utils::Array<double,3>& temp);



} // namespace
} // namespace


#endif
