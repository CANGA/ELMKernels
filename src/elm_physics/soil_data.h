
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
template <typename ArrayI1, typename ArrayD2>
void read_soil_colors(const Utils::DomainDecomposition<2>& dd, const std::string& filename, ArrayI1 isoicol,
                      ArrayD2 albsat, ArrayD2 albdry);

template <typename ArrayD2>
void read_soil_texture(const Utils::DomainDecomposition<2>& dd, const std::string& filename, ArrayD2 pct_sand,
                       ArrayD2 pct_clay, ArrayD2 organic);

} // namespace ELM::read_soil

#include "soil_data_impl.hh"
