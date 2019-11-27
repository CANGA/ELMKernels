// Example of what the reader interface might look like?

#include "array.hh"

namespace ELM {
namespace IO {

//
// Returns shape in file, as a pair (ntimes, NX_GLOBAL, NY_GLOBAL)
//
std::pair<std::size_t, std::size_t, std::size_t> get_dimensions(const std::string& fname, std::size_t start_year, std::size_t start_months, std::size_t n_months);

//
// Returns number of phenology times (likely always 12?)
//
std::size_t
get_phenology_num_times(const std::string& fname);

//
// Assumes shape:
//  ( get_phenology_num_times(), NX, NY, N_PFTS )
//
// offset provides an offset in the first dimension (num_times) to allow
// reading in multiple years in multiple files via sequential calls.


// I have to think about this one
void read_phenology(const std::string& fname,
                    std::size_t offset,
                    Array<4,double>& lai,
                    Array<4,double>& sai);

//
// Assumes shape:
//   ( get_forcing_num_times(), NX, NY )
//
void read_forcing(const std::string& fname, std::size_t start_year, std::size_t start_months, std::size_t n_months, std::size_t i_beg, std::size_t j_beg);

	

	//I would need nx, ny and times from shape
	//temperature and precipitation (rain, snow) inside this subroutine
} // namespace
} // namespace
