
#pragma once

namespace ELM::read_soil {

template <typename ArrayD2>
constexpr void get_albsat(int& mxsoil_color, ArrayD2 albsat)
{
  if (mxsoil_color == 8) {
    albsat(0, 0) = 12.0;
    albsat(0, 1) = 0.24;
    albsat(1, 0) = 0.11;
    albsat(1, 1) = 0.22;
    albsat(2, 0) = 0.10;
    albsat(2, 1) = 0.20;
    albsat(3, 0) = 0.09;
    albsat(3, 1) = 0.18;
    albsat(4, 0) = 0.08;
    albsat(4, 1) = 0.16;
    albsat(5, 0) = 0.07;
    albsat(5, 1) = 0.14;
    albsat(6, 0) = 0.06;
    albsat(6, 1) = 0.12;
    albsat(7, 0) = 0.05;
    albsat(7, 1) = 0.10;
  } else if (mxsoil_color == 20) {
    albsat(0, 0) = 0.25;
    albsat(0, 1) = 0.50;
    albsat(1, 0) = 0.23;
    albsat(1, 1) = 0.46;
    albsat(2, 0) = 0.21;
    albsat(2, 1) = 0.42;
    albsat(3, 0) = 0.20;
    albsat(3, 1) = 0.40;
    albsat(4, 0) = 0.19;
    albsat(4, 1) = 0.38;
    albsat(5, 0) = 0.18;
    albsat(5, 1) = 0.36;
    albsat(6, 0) = 0.17;
    albsat(6, 1) = 0.34;
    albsat(7, 0) = 0.16;
    albsat(7, 1) = 0.32;
    albsat(8, 0) = 0.15;
    albsat(8, 1) = 0.30;
    albsat(9, 0) = 0.14;
    albsat(9, 1) = 0.28;
    albsat(10, 0) = 0.13;
    albsat(10, 1) = 0.26;
    albsat(11, 0) = 0.12;
    albsat(11, 1) = 0.24;
    albsat(12, 0) = 0.11;
    albsat(12, 1) = 0.22;
    albsat(13, 0) = 0.10;
    albsat(13, 1) = 0.20;
    albsat(14, 0) = 0.09;
    albsat(14, 1) = 0.18;
    albsat(15, 0) = 0.08;
    albsat(15, 1) = 0.16;
    albsat(16, 0) = 0.07;
    albsat(16, 1) = 0.14;
    albsat(17, 0) = 0.06;
    albsat(17, 1) = 0.12;
    albsat(18, 0) = 0.05;
    albsat(18, 1) = 0.10;
    albsat(19, 0) = 0.04;
    albsat(19, 1) = 0.08;
  } else {
    throw std::runtime_error("ELM ERROR: mxsoil_color must be 8 or 20");
  }
}

template <typename ArrayD2>
constexpr void get_albdry(int& mxsoil_color, ArrayD2 albdry)
{
  if (mxsoil_color == 8) {
    albdry(0, 0) = 0.24;
    albdry(0, 1) = 0.48;
    albdry(1, 0) = 0.22;
    albdry(1, 1) = 0.44;
    albdry(2, 0) = 0.20;
    albdry(2, 1) = 0.40;
    albdry(3, 0) = 0.18;
    albdry(3, 1) = 0.36;
    albdry(4, 0) = 0.16;
    albdry(4, 1) = 0.32;
    albdry(5, 0) = 0.14;
    albdry(5, 1) = 0.28;
    albdry(6, 0) = 0.12;
    albdry(6, 1) = 0.24;
    albdry(7, 0) = 0.10;
    albdry(7, 1) = 0.20;
  } else if (mxsoil_color == 20) {
    albdry(0, 0) = 0.36;
    albdry(0, 1) = 0.61;
    albdry(1, 0) = 0.34;
    albdry(1, 1) = 0.57;
    albdry(2, 0) = 0.32;
    albdry(2, 1) = 0.53;
    albdry(3, 0) = 0.31;
    albdry(3, 1) = 0.51;
    albdry(4, 0) = 0.30;
    albdry(4, 1) = 0.49;
    albdry(5, 0) = 0.29;
    albdry(5, 1) = 0.48;
    albdry(6, 0) = 0.28;
    albdry(6, 1) = 0.45;
    albdry(7, 0) = 0.27;
    albdry(7, 1) = 0.43;
    albdry(8, 0) = 0.26;
    albdry(8, 1) = 0.41;
    albdry(9, 0) = 0.25;
    albdry(9, 1) = 0.39;
    albdry(10, 0) = 0.24;
    albdry(10, 1) = 0.37;
    albdry(11, 0) = 0.23;
    albdry(11, 1) = 0.35;
    albdry(12, 0) = 0.22;
    albdry(12, 1) = 0.33;
    albdry(13, 0) = 0.20;
    albdry(13, 1) = 0.31;
    albdry(14, 0) = 0.18;
    albdry(14, 1) = 0.29;
    albdry(15, 0) = 0.16;
    albdry(15, 1) = 0.27;
    albdry(16, 0) = 0.14;
    albdry(16, 1) = 0.25;
    albdry(17, 0) = 0.12;
    albdry(17, 1) = 0.23;
    albdry(18, 0) = 0.10;
    albdry(18, 1) = 0.21;
    albdry(19, 0) = 0.08;
    albdry(19, 1) = 0.16;
  } else {
    throw std::runtime_error("ELM ERROR: mxsoil_color must be 8 or 20");
  }
}

// serial I/O function
template <typename ArrayI1, typename ArrayD2>
void read_soil_colors(const Utils::DomainDecomposition<2>& dd,
                      const std::string& filename, ArrayI1 isoicol,
                      ArrayD2 albsat, ArrayD2 albdry)
{
  // get soil color
  {
    // get file start idx and size to read
    std::array<size_t, 2> start = {dd.start[0], dd.start[1]};
    std::array<size_t, 2> count = {dd.n_local[0], dd.n_local[1]};

    // read data
    Array<int, 2> arr_for_read(dd.n_local[0], dd.n_local[1]);
    IO::read_netcdf(dd.comm, filename, "SOIL_COLOR", start, count, arr_for_read.data());

    // place data into [ncells] order
    for (int i = 0; i != static_cast<int>(dd.n_local[0]); ++i) {
      for (int j = 0; j != static_cast<int>(dd.n_local[1]); ++j) {
        isoicol(i * dd.n_local[1] + j) = arr_for_read(i, j);
      }
    }
  }

  // get mxsoil_color
  int mxsoil_color;
  {
    std::array<size_t, 1> start = {0};
    std::array<size_t, 1> count = {1};
    // IO::read(dd.comm, filename, "mxsoil_color", start, count, mxsoil_color.data());
    IO::read_netcdf(dd.comm, filename, "mxsoil_color", start, count, &mxsoil_color);
  }

  // resize albsat and albdry if needed
  if (albsat.extent(0) != mxsoil_color) {
    NS::resize(albsat, mxsoil_color, numrad);
  }
  if (albdry.extent(0) != mxsoil_color) {
    NS::resize(albdry, mxsoil_color, numrad);
  }

  // get correct albsat and albdry arrays based on mxsoil_color
  get_albsat(mxsoil_color, albsat);
  get_albdry(mxsoil_color, albdry);
}

template <typename ArrayD2>
void read_soil_texture(const Utils::DomainDecomposition<2>& dd, const std::string& filename, ArrayD2 pct_sand,
                       ArrayD2 pct_clay, ArrayD2 organic) {
  // get pct_sand and pct_clay
  // get file start idx and size to read
  std::array<size_t, 3> start = {0, dd.start[0], dd.start[1]};
  std::array<size_t, 3> count = {ELM::nlevsoi, dd.n_local[0], dd.n_local[1]};

  // read pct_sand
  Array<double, 3> arr_for_read(ELM::nlevsoi, dd.n_local[0], dd.n_local[1]);
  IO::read_netcdf(dd.comm, filename, "PCT_SAND", start, count, arr_for_read.data());

  // place data into [ncells, nlevsoi] order
  for (int i = 0; i != static_cast<int>(dd.n_local[0]); ++i) {
    for (int j = 0; j != static_cast<int>(dd.n_local[1]); ++j) {
      for (int k = 0; k != ELM::nlevsoi; ++k) {
        pct_sand(i * dd.n_local[1] + j, k) = arr_for_read(k, j, i);
      }
    }
  }

  // read pct_clay
  IO::read_netcdf(dd.comm, filename, "PCT_CLAY", start, count, arr_for_read.data());
  // place data into [ncells, nlevsoi] order
  for (int i = 0; i != static_cast<int>(dd.n_local[0]); ++i) {
    for (int j = 0; j != static_cast<int>(dd.n_local[1]); ++j) {
      for (int k = 0; k != ELM::nlevsoi; ++k) {
        pct_clay(i * dd.n_local[1] + j, k) = arr_for_read(k, j, i);
      }
    }
  }

  // read organic
  IO::read_netcdf(dd.comm, filename, "ORGANIC", start, count, arr_for_read.data());
  // place data into [ncells, nlevsoi] order
  for (int i = 0; i != static_cast<int>(dd.n_local[0]); ++i) {
    for (int j = 0; j != static_cast<int>(dd.n_local[1]); ++j) {
      for (int k = 0; k != ELM::nlevsoi; ++k) {
        organic(i * dd.n_local[1] + j, k) = arr_for_read(k, j, i);
      }
    }
  }
}

} // namespace ELM::read_soil
