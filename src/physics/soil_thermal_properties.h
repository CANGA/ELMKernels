
#pragma once

#include "elm_constants.h"
#include "land_data.h"

#include <cmath>

#include "compile_options.hh"


namespace ELM::soil_thermal {

  namespace detail {
    static constexpr double TKICE{2.290}; // thermal conductivity of ice   [W/m/K]
    static constexpr double TKWAT{0.57};  // thermal conductivity of water [W/m/K]
    static constexpr double TKBDRK{3.0};  // thermal conductivity of 'typical' saturated granitic rock (Clauser and Huenges, 1995)(W/m/K)
    static constexpr double THIN_SFCLAYER{1.0e-6};   // Threshold for thin surface layer
  }

/*
const int& ltype,
const ArrayD1 h2osoi_liq  [nlevgrnd + nlevsno]
const ArrayD1 h2osoi_ice [nlevgrnd + nlevsno]
const ArrayD1 t_soisno [nlevgrnd + nlevsno]
const ArrayD1 dz  [nlevgrnd + nlevsno]
const ArrayD1 watsat, [nlevgrnd]
tkmg [nlevgrnd]
thk  [nlevgrnd + nlevsno]
tkdry [nlevgrnd]


*/
template <typename ArrayD2>
ACCELERATE
void calc_soil_tk(const int& c,
                  const int& ltype,
                  const ArrayD2 h2osoi_liq,
                  const ArrayD2 h2osoi_ice,
                  const ArrayD2 t_soisno,
                  const ArrayD2 dz,
                  const ArrayD2 watsat,
                  const ArrayD2 tkmg,
                  const ArrayD2 tkdry,
                  ArrayD2 thk);


template <typename ArrayD2>
ACCELERATE
void calc_snow_tk(const int& c,
                  const int& snl,
                  const double& frac_sno,
                  const ArrayD2 h2osoi_liq,
                  const ArrayD2 h2osoi_ice,
                  const ArrayD2 dz,
                  ArrayD2 thk);


template <typename ArrayD2>
ACCELERATE
void calc_face_tk(const int& c,
                  const int& snl,
                  const ArrayD2 thk,
                  const ArrayD2 z,
                  const ArrayD2 zi,
                  ArrayD2 tk);


template <typename ArrayD2>
ACCELERATE
void calc_soil_heat_capacity(const int& c,
                             const int& ltype,
                             const int& snl,
                             const double& h2osno,
                             const ArrayD2 watsat,
                             const ArrayD2 h2osoi_ice,
                             const ArrayD2 h2osoi_liq,
                             const ArrayD2 dz,
                             const ArrayD2 csol,
                             ArrayD2 cv);


template <typename ArrayD2>
ACCELERATE
void calc_snow_heat_capacity(const int& c,
                             const int& snl,
                             const double& frac_sno,
                             const ArrayD2 h2osoi_ice,
                             const ArrayD2 h2osoi_liq,
                             ArrayD2 cv);


template <typename ArrayD2>
ACCELERATE
double calc_h2osfc_tk(const int& c,
                      const double& h2osfc,
                      const ArrayD2 thk,
                      const ArrayD2 z);


ACCELERATE
double calc_h2osfc_heat_capacity(const int& snl, const double& h2osfc, const double& frac_h2osfc);

// this should be somewhere else
ACCELERATE
double calc_h2osfc_height(const int& snl, const double& h2osfc, const double& frac_h2osfc);

} // namespace ELM::soil_thermal

#include "soil_thermal_properties_impl.hh"
