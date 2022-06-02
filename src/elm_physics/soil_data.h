
#pragma once

#include "array.hh"
#include "elm_constants.h"
#include "read_input.hh"
#include "utils.hh"

#include <array>
#include <string>

#include "kokkos_includes.hh"

namespace ELM::read_soil {

template <typename ArrayD2>
constexpr void get_albsat(int& mxsoil_color, ArrayD2 albsat);

template <typename ArrayD2>
constexpr void get_albdry(int& mxsoil_color, ArrayD2 albdry);

// serial I/O function
template <typename h_ArrayI1, typename h_ArrayD2>
void read_soil_colors(const Utils::DomainDecomposition<2>& dd, const std::string& filename, h_ArrayI1 isoicol,
                      h_ArrayD2 albsat, h_ArrayD2 albdry);

template <typename h_ArrayD1, typename h_ArrayD2>
void read_soil_texture(const Utils::DomainDecomposition<2>& dd, const std::string& fname_surfdata,
                       const std::string& fname_param, h_ArrayD1 organic_max, h_ArrayD2 pct_sand,
                       h_ArrayD2 pct_clay, h_ArrayD2 organic);

} // namespace ELM::read_soil

#include "soil_data_impl.hh"
