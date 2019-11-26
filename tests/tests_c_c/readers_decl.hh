// Example of what the reader interface might look like?

#include "array.hh"

namespace ELM {
namespace IO {

//
// Returns shape in file, as a pair (NX, NY)
//
std::pair<std::size_t, std::size_t>
get_num_grid_cells(const std::string& fname);

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
void read_phenology(const std::string& fname,
                    std::size_t offset,
                    Array<4,double>& lai,
                    Array<4,double>& sai);

//
// Number of forcing times (likely always 365 * 8?)
std::size_t
get_forcing_num_times(const std::string& fname);

//
// Assumes shape:
//   ( get_forcing_num_times(), NX, NY )
//
void read_forcing(const std::string& fname,
                  const std::string& varname,
                  std::size_t offset,
                  Array<3,double>& forc);

} // namespace
} // namespace
