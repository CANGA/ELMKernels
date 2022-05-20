

#pragma once

namespace ELM::soil_temperature::detail {

  ACCELERATE
  double calc_lwrad_emit(const double& emg, const double& t_grnd)
  {
    return emg * ELMconst::STEBOL * pow(t_grnd, 4.0);
  }

  ACCELERATE
  double calc_dlwrad_emit(const double& emg, const double& t_grnd)
  {
    return 4.0 * emg * ELMconst::STEBOL * pow(t_grnd, 3.0);
  }

  // t_soisno is t_soisno(nlevsno-snl)
  ACCELERATE
  double calc_lwrad_emit_snow(const double& emg, const double& t_soisno)
  {
    return emg * ELMconst::STEBOL * pow(t_soisno, 4.0);
  }

  // t_soisno is t_soisno(nlevsno)
  ACCELERATE
  double calc_lwrad_emit_soil(const double& emg, const double& t_soisno)
  {
    return emg * ELMconst::STEBOL * pow(t_soisno, 4.0);
  }

  ACCELERATE
  double calc_lwrad_emit_h2osfc(const double& emg, const double& t_h2osfc)
  {
    return emg * ELMconst::STEBOL * pow(t_h2osfc, 4.0);
  }

} // namespace ELM::soil_temperature::detail

namespace ELM::soil_temperature {

  // t_soisno is t_soisno(nlevsno)
  // return value is also ELM's elfx_gnet_soil without any pft weighting
  ACCELERATE
  double calc_hs_soil(const int& frac_veg_nosno,
                      const double& t_soisno,
                      const double& sabg_soil,
                      const double& dlrad,
                      const double& emg,
                      const double& forc_lwrad,
                      const double& eflx_sh_soil,
                      const double& qflx_ev_soil,
                      const double& htvp)
  {
    return sabg_soil + dlrad + (1.0 - frac_veg_nosno) * emg * forc_lwrad -
      detail::calc_lwrad_emit_soil(emg, t_soisno) - (eflx_sh_soil + qflx_ev_soil * htvp);
  }

  // return value is also ELM's elfx_gnet_h2osfc without any pft weighting
  ACCELERATE
  double calc_hs_h2osfc(const int& frac_veg_nosno,
                        const double& t_h2osfc,
                        const double& sabg_soil,
                        const double& dlrad,
                        const double& emg,
                        const double& forc_lwrad,
                        const double& eflx_sh_soil,
                        const double& qflx_ev_soil,
                        const double& htvp)
  {
    return sabg_soil + dlrad + (1.0 - frac_veg_nosno) * emg * forc_lwrad -
      detail::calc_lwrad_emit_h2osfc(emg, t_h2osfc) - (eflx_sh_soil + qflx_ev_soil * htvp);
  }

  // sabg_lyr is sabg_lyr(nlevsno-snl)
  ACCELERATE
  double calc_hs_top(const int& frac_veg_nosno,
                     const double& t_grnd,
                     const double& sabg_lyr,
                     const double& dlrad,
                     const double& emg,
                     const double& forc_lwrad,
                     const double& eflx_sh_grnd,
                     const double& qflx_evap_soi,
                     const double& htvp)
  {
    return sabg_lyr + dlrad + (1.0 - frac_veg_nosno) * emg * forc_lwrad -
      detail::calc_lwrad_emit(emg, t_grnd) - (eflx_sh_grnd + qflx_evap_soi * htvp);
  }

  // t_soisno is t_soisno(nlevsno-snl)
  // sabg_lyr is sabg_lyr(nlevsno-snl)
  // return value is also ELM's elfx_gnet_snow without any pft weighting
  ACCELERATE
  double calc_hs_top_snow(const int& frac_veg_nosno,
                          const double& t_soisno,
                          const double& sabg_lyr,
                          const double& dlrad,
                          const double& emg,
                          const double& forc_lwrad,
                          const double& eflx_sh_snow,
                          const double& qflx_ev_snow,
                          const double& htvp)
  {
    return sabg_lyr + dlrad + (1.0 - frac_veg_nosno) * emg * forc_lwrad -
      detail::calc_lwrad_emit_snow(emg, t_soisno) - (eflx_sh_snow + qflx_ev_snow * htvp);
  }

  ACCELERATE
  double calc_dhsdT (const double& cgrnd, const double& emg, const double& t_grnd)
  {
    return -cgrnd - detail::calc_dlwrad_emit(emg, t_grnd);
  }

  ACCELERATE
  double check_absorbed_solar(const double& frac_sno_eff, const double& sabg_snow, const double& sabg_soil)
  {
    return frac_sno_eff * sabg_snow + (1.0 - frac_sno_eff ) * sabg_soil;
  }


  template <typename ArrayI1, typename ArrayD1>
  ACCELERATE
  void calc_diffusive_heat_flux(const int& snl,
                                const ArrayD1 tk,
                                const ArrayD1 t_soisno,
                                const ArrayD1 z,
                                ArrayD1 fn)
  {
    using ELMdims::nlevgrnd;
    using ELMdims::nlevsno;

    for (int i = 0; i < nlevgrnd + nlevsno; ++i) {
      if (i >= nlevsno - snl) { // at or below top of snow/surface
        if (i < nlevgrnd + nlevsno - 1) { // not the bottom
          fn(i) = tk(i) * (t_soisno(i+1) - t_soisno(i)) / (z(i+1)-z(i));
        } else { // the bottom
          fn(i) = 0.0; // hardwired as 0 for now, maybe add eflx_bot in the future (it's default 0 in ELM)
        }
      }
    }
  }

  template <typename ArrayI1, typename ArrayD1>
  ACCELERATE
  void calc_heat_flux_matrix_factor(const int& snl,
                                    const double& dtime,
                                    const ArrayD1 cv,
                                    const ArrayD1 dz,
                                    const ArrayD1 z,
                                    const ArrayD1 zi,
                                    ArrayD1 fact)
  {
    static constexpr double capr{0.34}; // Tuning factor to turn first layer T into surface T
    using ELMdims::nlevgrnd;
    using ELMdims::nlevsno;

    for (int i = 0; i < nlevgrnd + nlevsno; ++i) {
      if (i == nlevsno - snl) {
        fact(i) = dtime / cv(i) * dz(i) / (0.5 * (z(i) - zi(i-1) + capr * (z(i+1) - zi(i-1))));
      } else {
        fact(i) = dtime / cv(i);
      }
    }
  }

} // namespace ELM::soil_temperature
