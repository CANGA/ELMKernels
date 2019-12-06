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
// Assumes shape(arr) == (N_MONTHS, NX_LOCAL, NY_LOCAL, N_PFT)
//  
void
read_phenology(const MPI_Comm& comm,
               const std::string& dir,const std::string& basename, const std::string& phenology_type,
               int start_year, int start_month,
               int i_beg, int j_beg, ELM::Utils::Array<double,4>& arr);


//
// Read a variable from the forcing files.
//
// Assumes shape(arr) == (N_TIMES, NX_LOCAL, NY_LOCAL )
//
void
read_forcing(const MPI_Comm& comm,
             const std::string& dir,const std::string& basename, const std::string& forcing_type,
             int start_year, int start_month, int n_months, 
             int i_beg, int j_beg, ELM::Utils::Array<double,3>& arr);


// Convert precipitation to rain and snow
void convert_precip_to_rain_snow(ELM::Utils::Array<double,3>& rain, ELM::Utils::Array<double,3>& snow, 
											ELM::Utils::Array<double,3>& temp);



} // namespace
} // namespace


#endif
