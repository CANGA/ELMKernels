
#pragma once 

#include "elm_constants.h"
#include "utils.hh"
#include "array.hh"
#include "read_input.hh"

#include <string>
#include <array>

namespace ELM::read_soil {

template <typename ArrayD2>
constexpr void get_albsat(int mxsoil_color, ArrayD2& albsat);

template <typename ArrayD2>
constexpr void get_albdry(int mxsoil_color, ArrayD2& albdry);

// serial I/O function
template <typename ArrayI1, typename ArrayD2>
void read_soil_colors(const Utils::DomainDecomposition<2>& dd, const std::string& filename, ArrayI1& isoicol, ArrayD2& albsat, ArrayD2& albdry);

template <typename ArrayD2>
void read_soil_texture(const Utils::DomainDecomposition<2>& dd, const std::string& filename, ArrayD2& pct_sand, ArrayD2& pct_clay, ArrayD2& organic);

} // namespace ELM::read_soil

#include "read_soil_colors_impl.hh"
