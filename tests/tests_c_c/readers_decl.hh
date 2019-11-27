// Example of what the reader interface might look like?
namespace ELM {
namespace IO {

//
// Returns shape in file, as a pair (N_TIMES, NX_GLOBAL, NY_GLOBAL)
//
std::tuple<std::size_t, std::size_t, std::size_t>
get_dimensions(const std::string& fname,
               std::size_t start_year, std::size_t start_months, std::size_t n_months);


//
// Read a variable from the phenology file.
//
// Assumes shape(arr) == (N_MONTHS, NX_LOCAL, NY_LOCAL, N_PFT)
//  
template<class Array_t>
void
read_phenology(const std::string& fname, const std::string& phenology_type,
               std::size_t start_year, std::size_t start_month,
               std::size_t i_beg, std::size_t j_beg, Array_t arr);


//
// Read a variable from the forcing files.
//
// Assumes shape(arr) == (N_TIMES, NX_LOCAL, NY_LOCAL )
//
template<class Array_t>
void
read_forcing(const std::string& fname, const std::string& forcing_type,
             std::size_t start_year, std::size_t start_month,
             std::size_t i_beg, std::size_t j_beg, Array_t& forcing);

} // namespace
} // namespace
