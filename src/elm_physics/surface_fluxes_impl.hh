
#pragma once

namespace ELM::surface_fluxes {

// Soil Energy balance check
template<typename ArrayD1>
double soil_energy_balance(const bool& ctype, const int& snl, const double& eflx_soil_grnd, const double& xmf, const double& xmf_h2osfc, const double& frac_h2osfc, const double& t_h2osfc,
  const double& t_h2osfc_bef, const double& dtime, const double& eflx_h2osfc_to_snow, const double& eflx_building_heat, const double&frac_sno_eff, const ArrayD1 t_soisno, const ArrayD1 tssbef, const ArrayD1 fact) {

  double errsoi = eflx_soil_grnd - xmf - xmf_h2osfc - frac_h2osfc * (t_h2osfc - t_h2osfc_bef) * (t_h2osfc / dtime);
  errsoi += eflx_h2osfc_to_snow;
  // For urban sunwall, shadewall, and roof columns, the "soil" energy balance check
  // must include the heat flux from the interior of the building.
  if (ctype == icol_sunwall || ctype == icol_shadewall || ctype == icol_roof) {
    errsoi += eflx_building_heat;
  }
  for (int j = 0; j < ELM::nlevgrnd+ELM::nlevsno; ++j) {
    if ((ctype != icol_sunwall && ctype != icol_shadewall && ctype != icol_roof) || ( j < ELM::nlevurb)) {
      // area weight heat absorbed by snow layers
      if (j >= ELM::nlevsno-snl && j < ELM::nlevsno) { errsoi -= frac_sno_eff * (t_soisno(j) - tssbef(j)) / fact(j); }
      if (j >= ELM::nlevsno) { errsoi -= (t_soisno(j) - tssbef(j)) / fact(j); }
    }
  }
return errsoi;
}

} // namespace ELM::surface_fluxes
