
#pragma once

#include "aerosol_data_old.h"
#include "array.hh"
#include "elm_constants.h"

#include <functional>
#include <tuple>

#include "compile_options.hh"

namespace ELM {

template <typename ArrayD2>
struct AerosolMasses {
  const int nlevsno_{ELMdims::nlevsno()};
  AerosolMasses(const size_t& ncells);
  ~AerosolMasses() = default;
  ArrayD2 mss_bcphi, mss_bcpho, mss_dst1, mss_dst2, mss_dst3, mss_dst4;
};

template <typename ArrayD2>
struct AerosolConcentrations {
  const int nlevsno_{ELMdims::nlevsno()};
  AerosolConcentrations(const size_t& ncells);
  ~AerosolConcentrations() = default;
  ArrayD2 mss_cnc_bcphi, mss_cnc_bcpho, mss_cnc_dst1, mss_cnc_dst2, mss_cnc_dst3, mss_cnc_dst4;
};

}

namespace ELM::aerosols {

// functor to calculate aerosol deposition
// this is from AerosolMod.F90
// should be run after inter-layer aerosol fluxes are accounted for in
// SnowHydrologyMod.F90::SnowWater()
template <typename T, typename ArrayI1, typename ArrayD2>
struct ComputeAerosolDeposition {
  ComputeAerosolDeposition(const T& aerosol_forc, const ArrayI1 snl, AerosolMasses<ArrayD2>& aerosol_masses);

  ACCELERATE
  void operator()(const int i) const;

private:
  ArrayI1 snl_;
  AerosolMasses<ArrayD2> aerosol_masses_;
  double forc_bcphi_, forc_bcpho_, forc_dst1_, forc_dst2_, forc_dst3_, forc_dst4_;
};

// functor to calculate aerosol mass and concentration in snow layers
// this is from AerosolMod.F90
// gets run in ELM directly after hydrology is called
template <typename ArrayB1, typename ArrayI1, typename ArrayD1, typename ArrayD2>
struct ComputeAerosolConcenAndMass {
  ComputeAerosolConcenAndMass(const double& dtime, const ArrayB1 do_capsnow, const ArrayI1 snl,
                              const ArrayD2 h2osoi_liq, const ArrayD2 h2osoi_ice, const ArrayD2 snw_rds,
                              const ArrayD1 qflx_snwcp_ice, AerosolMasses<ArrayD2>& aerosol_masses,
                              AerosolConcentrations<ArrayD2>& aerosol_concentrations);

  ACCELERATE
  void operator()(const int i) const;

private:
  double dtime_;
  ArrayB1 do_capsnow_;
  ArrayI1 snl_;
  ArrayD2 h2osoi_liq_;
  ArrayD2 h2osoi_ice_;
  ArrayD2 snw_rds_;
  ArrayD1 qflx_snwcp_ice_;
  AerosolMasses<ArrayD2> aerosol_masses_;
  AerosolConcentrations<ArrayD2> aerosol_concentrations_;
};

} // namespace ELM::aerosols

#include "aerosol_physics_old_impl.hh"
