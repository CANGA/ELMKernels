
#pragma once

#include "aerosol_data.h"
#include "elm_constants.h"
#include "array.hh"

#include <functional>

namespace ELM::aerosols {

// functor to calculate aerosol deposition
// this is from AerosolMod.F90
// should be run after inter-layer aerosol fluxes are accounted for in 
// SnowHydrologyMod.F90::SnowWater()
template<typename T, typename ArrayI1>
struct ComputeAerosolDeposition {
  ComputeAerosolDeposition(const T& aerosol_forc, const ArrayI1& snl, AerosolMasses& aerosol_masses);

  void operator()(const int i) const;

private:
  ArrayI1 snl_;
  AerosolMasses aerosol_masses_;
  double forc_bcphi_, forc_bcpho_, forc_dst1_, forc_dst2_, forc_dst3_, forc_dst4_;
};


// functor to calculate aerosol mass and concentration in snow layers
// this is from AerosolMod.F90
// gets run in ELM directly after hydrology is called
template<typename ArrayI1, typename ArrayD1, typename ArrayD2>
struct ComputeAerosolConcenAndMass {
  ComputeAerosolConcenAndMass(const bool& do_capsnow, const double& dtime, const ArrayI1& snl, const ArrayD2& h2osoi_liq,
    const ArrayD2& h2osoi_ice, const ArrayD2& snw_rds, const ArrayD1& qflx_snwcp_ice, AerosolMasses& aerosol_masses,
    AerosolConcentrations& aerosol_concentrations);

  void operator()(const int i) const;

private:
  bool do_capsnow_;
  double dtime_;
  ArrayI1 snl_;
  ArrayD2 h2osoi_liq_;
  ArrayD2 h2osoi_ice_;
  ArrayD2 snw_rds_;
  ArrayD1 qflx_snwcp_ice_;
  AerosolMasses aerosol_masses_;
  AerosolConcentrations aerosol_concentrations_;
};

// convenience function to invoke aerosol mass and concen functor
template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
void invoke_aerosol_concen_and_mass(const bool& do_capsnow, const double& dtime, const ArrayI1& snl, const ArrayD2& h2osoi_liq,
  const ArrayD2& h2osoi_ice, const ArrayD2& snw_rds, const ArrayD1& qflx_snwcp_ice, AerosolMasses& aerosol_masses,
  AerosolConcentrations& aerosol_concentrations);

} // namespace ELM::aerosols

#include "aerosol_physics_impl.hh"


