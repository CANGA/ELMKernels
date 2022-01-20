
#pragma once

namespace ELM::snicar_data {

template <typename ArrayT, typename T, size_t D>
inline void read_and_fill_array(const Comm_type &comm, const std::string &filename, const std::string &varname,
                 Array<T, D> &arr_for_read, ArrayT& arrout) {
  std::array<size_t, D> start{0};
  std::array<size_t, D> count;
  for (int i = 0; i != D; ++i)
    count[i] = static_cast<size_t>(arr_for_read.extent(i));
  IO::read_netcdf(comm, filename, varname, start, count, arr_for_read.data());
  if constexpr (D == 1) {
    for (int i = 0; i < arr_for_read.extent(0); ++i) 
      arrout(i) = arr_for_read(i);
  } else if constexpr (D == 2) {
    for (int i = 0; i < arr_for_read.extent(0); ++i)
      for (int j = 0; j < arr_for_read.extent(1); ++j)
        arrout(i,j) = arr_for_read(i,j);
  } else if constexpr (D == 3) {
    for (int i = 0; i < arr_for_read.extent(0); ++i)
      for (int j = 0; j < arr_for_read.extent(1); ++j)
        for (int k = 0; k < arr_for_read.extent(2); ++k)
          arrout(i,j,k) = arr_for_read(i,j,k);
  }
  return;
}

} // namespace ELM::snicar_data
