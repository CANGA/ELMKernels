

#pragma once

#include "pentadiagonal_solver.h"
#include "soil_temp_rhs.h"
#include "soil_temp_lhs.h"
#include "soil_thermal_properties.h"
#include "helper_functions.hh"
#include "phase_change.h"
#include "invoke_kernel.hh"

namespace ELM::soil_temp {

  ACCELERATE
  double calc_surface_heat_flux(const int& frac_veg_nosno,
                                const double& dlrad,
                                const double& emg,
                                const double& forc_lwrad,
                                const double& htvp,
                                const double& solar_abg,
                                const double& temp,
                                const double& eflx_sh,
                                const double& qflx_ev)
  {
    return solar_abg + dlrad + (1.0 - frac_veg_nosno) * emg * forc_lwrad -
      calc_lwrad_emit(emg, temp) - (eflx_sh + qflx_ev * htvp);
  }

  ACCELERATE
  double calc_dhsdT(const double& cgrnd, const double& emg, const double& t_grnd)
  {
    return -cgrnd - calc_dlwrad_emit(emg, t_grnd);
  }

  ACCELERATE
  double check_absorbed_solar(const double& frac_sno_eff, const double& sabg_snow, const double& sabg_soil)
  {
    return frac_sno_eff * sabg_snow + (1.0 - frac_sno_eff ) * sabg_soil;
  }


  // calculates fn, the diffusive heat flux through layer interfaces
  // fn[nlevgrnd()+nlevsno()]
  // fn(i) is the interface between cells i and i+1
  // this is different than zi, where zi(i) is between cells i-1 and i 
  template <typename ArrayD1>
  ACCELERATE
  void calc_diffusive_heat_flux(const int& snl,
                                const ArrayD1 tk,
                                const ArrayD1 t_soisno,
                                const ArrayD1 z,
                                ArrayD1 fn)
  {
    using ELMdims::nlevgrnd;
    using ELMdims::nlevsno;

    // zero out inactive layer interfaces
    const int top = nlevsno() - snl;
    for (int i = 0; i < top; ++i) {
      fn(i) = 0.0;
    }

    // all active layer interfaces above bottom interface
    for (int i = top; i < nlevgrnd() + nlevsno() - 1; ++i) {
      fn(i) = tk(i) * (t_soisno(i+1) - t_soisno(i)) / (z(i+1)-z(i));
    }

    // bottom layer flux
    // hardwired as 0 for now, both here and in ELM
    // could add external eflx_bot in the future
    // eg to represent geothermal heat or
    // coupling condition from driving model  
    fn(nlevgrnd() + nlevsno() - 1) = 0.0;
  }


  ACCELERATE
  double calc_lwrad_emit(const double& emg, const double& temp)
  {
    return emg * ELMconst::STEBOL() * pow(temp, 4.0);
  }


  ACCELERATE
  double calc_dlwrad_emit(const double& emg, const double& t_grnd)
  {
    return 4.0 * emg * ELMconst::STEBOL() * pow(t_grnd, 3.0);
  }


  template <typename ArrayD1>
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

    // zero out inactive layers
    const int top = nlevsno() - snl;
    for (int i = 0; i < top; ++i) {
      fact(i) = 0.0;
    }

    // top active layer (snow or top subsurface layer if no snow)
    fact(top) = dtime / cv(top) * dz(top) / (0.5 * (z(top) -
        zi(top) + capr * (z(top+1) - zi(top))));

    // all layers below top active layer
    for (int i = top + 1; i < nlevgrnd() + nlevsno(); ++i) {
      fact(i) = dtime / cv(i);
    }
  }


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void set_tvector(const int& c,
                   const ArrayI1 snl,
                   const ArrayD1 t_h2osfc,
                   const ArrayD2 t_soisno,
                   ArrayD2 tvector)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;
    const int top = nlevsno() - snl(c);
    
    // zero inactive layers
    for (int i = 0; i < top; ++i)
      tvector(c, i) = 0.0;

    // snow layers
    for (int i = top; i < nlevsno(); ++i)
      tvector(c, i) = t_soisno(c, i);

    // surface water
    tvector(c, nlevsno()) = t_h2osfc(c);

    // soil
    for (int i = nlevsno(); i < nlevsno() + nlevgrnd(); ++i)
      tvector(c, i + 1) = t_soisno(c, i);
  }


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void update_temperature(const int& c,
                          const ArrayI1 snl,
                          const ArrayD1 frac_h2osfc,
                          const ArrayD2 tvector,
                          ArrayD1 t_h2osfc,
                          ArrayD2 t_soisno)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;
    const int top = nlevsno() - snl(c);

    // snow layers
    for (int i = top; i < nlevsno(); ++i)
      t_soisno(c, i) = tvector(c, i);

    // soil
    for (int i = nlevsno(); i < nlevsno() + nlevgrnd(); ++i)
      t_soisno(c, i) = tvector(c, i + 1);

    // surface water
    // ? first subsurface cell : surface water
    t_h2osfc(c) = (frac_h2osfc(c) != 0.0) ? tvector(c, nlevsno()) : t_soisno(c, nlevsno());
  }

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void update_t_grnd(const int& c,
                     const ArrayI1 snl,
                     const ArrayD1 frac_h2osfc,
                     const ArrayD1 frac_sno_eff,
                     const ArrayD1 t_h2osfc,
                     const ArrayD2 t_soisno,
                     ArrayD1 t_grnd)
  {
    using ELMdims::nlevsno;

    if (snl(c) > 0) {
      const int top = nlevsno() - snl(c);
      if(frac_h2osfc(c) != 0.0) {
        t_grnd(c) = frac_sno_eff(c) * t_soisno(c, top) + (1.0 - frac_sno_eff(c) - frac_h2osfc(c)) *
          t_soisno(c, nlevsno()) + frac_h2osfc(c) * t_h2osfc(c);
      } else {
        t_grnd(c) = frac_sno_eff(c) * t_soisno(c, top) + (1.0 - frac_sno_eff(c)) * t_soisno(c, nlevsno());
      }
    } else {
      if(frac_h2osfc(c) != 0.0) {
        t_grnd(c) = (1.0 - frac_h2osfc(c)) * t_soisno(c, nlevsno()) + frac_h2osfc(c) * t_h2osfc(c);
      } else {
        t_grnd(c) = t_soisno(c, nlevsno());
      }
    }
  }

} // namespace ELM::soil_temp
