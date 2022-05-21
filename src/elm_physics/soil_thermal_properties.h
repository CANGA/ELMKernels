
#pragma once

#include "elm_constants.h"

#include <cmath>

#include "kokkos_includes.hh"


namespace ELM::soil_thermal {

  namespace detail {
    static constexpr double TKICE = 2.290; // thermal conductivity of ice   [W/m/K]
    static constexpr double TKWAT = 0.57;  // thermal conductivity of water [W/m/K]
    static constexpr double TKBDRK = 3.0;  // thermal conductivity of 'typical' saturated granitic rock (Clauser and Huenges, 1995)(W/m/K)
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
template <typename ArrayD1>
ACCELERATE
void calc_soil_tk(const int& ltype,
                  const ArrayD1 h2osoi_liq,
                  const ArrayD1 h2osoi_ice,
                  const ArrayD1 t_soisno,
                  const ArrayD1 dz,
                  const ArrayD1 watsat,
                  const ArrayD1 tkmg,
                  const ArrayD1 tkdry,
                  ArrayD1 thk);


template <typename ArrayD1>
ACCELERATE
void calc_snow_tk(const int& snl,
                  const double& frac_sno,
                  const ArrayD1 h2osoi_liq,
                  const ArrayD1 h2osoi_ice,
                  const ArrayD1 dz,
                  ArrayD1 thk);


template <typename ArrayD1>
ACCELERATE
void calc_face_tk(const int& snl,
                  const ArrayD1 thk,
                  const ArrayD1 z,
                  const ArrayD1 zi,
                  ArrayD1 tk);

template <typename ArrayD1>
ACCELERATE
void calc_h2osfc_tk(const double& h2osfc,
                    const ArrayD1 thk,
                    const ArrayD1 z,
                    double& tk_h2osfc);


template <typename ArrayD1>
ACCELERATE
void calc_soil_heat_capacity(const int& ltype,
                             const int& snl,
                             const double& h2osno,
                             const ArrayD1 watsat,
                             const ArrayD1 h2osoi_ice,
                             const ArrayD1 h2osoi_liq,
                             const ArrayD1 dz,
                             const ArrayD1 csol,
                             ArrayD1 cv);


template <typename ArrayD1>
ACCELERATE
void calc_snow_heat_capacity(const int& snl,
                             const double& frac_sno,
                             const ArrayD1 h2osoi_ice,
                             const ArrayD1 h2osoi_liq,
                             ArrayD1 cv);

} // namespace ELM::soil_thermal

#include "soil_thermal_properties_impl.hh"
