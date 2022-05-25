

#pragma once

#include "pentadiagonal_solver.h"
#include "pentadiagonal_solver.h"
#include "soil_temp_rhs.h"
#include "soil_temp_lhs.h"
#include "soil_thermal_properties.h"
#include "helper_functions.hh"

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
      detail::calc_lwrad_emit(emg, temp) - (eflx_sh + qflx_ev * htvp);
  }

  ACCELERATE
  double calc_dhsdT(const double& cgrnd, const double& emg, const double& t_grnd)
  {
    return -cgrnd - detail::calc_dlwrad_emit(emg, t_grnd);
  }

  ACCELERATE
  double check_absorbed_solar(const double& frac_sno_eff, const double& sabg_snow, const double& sabg_soil)
  {
    return frac_sno_eff * sabg_snow + (1.0 - frac_sno_eff ) * sabg_soil;
  }


  // calculates fn, the diffusive heat flux through layer interfaces
  // fn[nlevgrnd+nlevsno]
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
    const int top = nlevsno - snl;
    for (int i = 0; i < top; ++i) {
      fn(i) = 0.0;
    }

    // all active layer interfaces above bottom interface
    for (int i = top; i < nlevgrnd + nlevsno - 1; ++i) {
      fn(i) = tk(i) * (t_soisno(i+1) - t_soisno(i)) / (z(i+1)-z(i));
    }

    // bottom layer flux
    // hardwired as 0 for now, both here and in ELM
    // could add external eflx_bot in the future
    // eg to represent geothermal heat or
    // coupling condition from driving model  
    fn(nlevgrnd + nlevsno - 1) = 0.0;
  }

  template <typename ArrayD3, typename ArrayI1, typename ArrayD1, typename ArrayD2>
  void solve_temperature(const double& dtime,
                         const ArrayI1 snl,
                         const ArrayI1 frac_veg_nosno,
                         const ArrayD1 dlrad,
                         const ArrayD1 emg,
                         const ArrayD1 forc_lwrad,
                         const ArrayD1 htvp,
                         const ArrayD1 cgrnd,
                         const ArrayD1 eflx_sh_soil,
                         const ArrayD1 qflx_ev_soil,
                         const ArrayD1 eflx_sh_h2osfc,
                         const ArrayD1 qflx_ev_h2osfc,
                         const ArrayD1 eflx_sh_grnd,
                         const ArrayD1 qflx_evap_soi,
                         const ArrayD1 eflx_sh_snow,
                         const ArrayD1 qflx_ev_snow,
                         const ArrayD1 frac_sno_eff,
                         const ArrayD1 frac_sno,
                         const ArrayD1 frac_h2osfc,
                         const ArrayD1 h2osno,
                         const ArrayD1 h2osfc,
                         const ArrayD1 sabg_snow,
                         const ArrayD1 sabg_soil,
                         const ArrayD2 sabg_lyr,
                         const ArrayD2 h2osoi_liq,
                         const ArrayD2 h2osoi_ice,
                         const ArrayD2 watsat,
                         const ArrayD2 tkmg,
                         const ArrayD2 tkdry,
                         const ArrayD2 csol,
                         const ArrayD2 dz,
                         const ArrayD2 zsoi,
                         const ArrayD2 zisoi,
                         ArrayD1 t_h2osfc,
                         ArrayD1 t_grnd,
                         ArrayD2 t_soisno,
                         ArrayD2 fact)
  {
    using ELMdims::nlevgrnd;
    using ELMdims::nlevsno;
    using Utils::create;
    using Utils::assign;

    const int ncells = snl.extent(0);

    // dummy ltype for now
    auto ltype = create<ArrayI1>("ltype", ncells);
    assign(ltype, 1);


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // call soil thermal properties kernels
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    auto thk = create<ArrayD2>("tk", ncells, nlevgrnd + nlevsno); // thermal conductivity of layer
    auto tk = create<ArrayD2>("tk", ncells, nlevgrnd + nlevsno); // thermal conductivity at layer interface
    auto cv = create<ArrayD2>("cv", ncells, nlevgrnd + nlevsno);
    auto tk_h2osfc = create<ArrayD1>("tk_h2osfc", ncells);
    auto c_h2osfc = create<ArrayD1>("c_h2osfc", ncells);
    auto dz_h2osfc = create<ArrayD1>("dz_h2osfc", ncells);
    auto soil_thermal_props = [=] (const int& c) {
      soil_thermal::calc_soil_tk(c, ltype(c), h2osoi_liq, h2osoi_ice, t_soisno, dz, watsat, tkmg, tkdry, thk);
      soil_thermal::calc_snow_tk(c, snl(c), frac_sno(c), h2osoi_liq, h2osoi_ice, dz, thk);
      soil_thermal::calc_face_tk(c, snl(c), thk, zsoi, zisoi, tk);
      soil_thermal::calc_soil_heat_capacity(c, ltype(c), snl(c), h2osno(c), watsat, h2osoi_ice, h2osoi_liq, dz, csol, cv);
      soil_thermal::calc_snow_heat_capacity(c, snl(c), frac_sno(c), h2osoi_ice, h2osoi_liq, cv);
      tk_h2osfc(c) = soil_thermal::calc_h2osfc_tk(c, h2osfc(c), thk, zsoi);
      c_h2osfc(c) = ELM::soil_thermal::calc_h2osfc_heat_capacity(snl(c), h2osfc(c), frac_h2osfc(c));
      dz_h2osfc(c) = ELM::soil_thermal::calc_h2osfc_height(snl(c), h2osfc(c), frac_h2osfc(c));
    };
    invoke_kernel(soil_thermal_props, std::make_tuple(ncells), "soil_thermal_props");


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // call pre-solve soil temperature kernels
    // surface heat fluxes
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    auto hs_soil = create<ArrayD1>("hs_soil", ncells); // [W/m2] soil heat flux
    auto hs_snow = create<ArrayD1>("hs_snow", ncells); // [W/m2] snow heat flux - maybe use for coupling?
    auto hs_h2osfc = create<ArrayD1>("hs_h2osfc", ncells); // [W/m2] standing water heat flux
    auto hs_top = create<ArrayD1>("hs_top", ncells); // [W/m2] net heat flux into surface layer
    auto hs_top_snow = create<ArrayD1>("hs_top_snow", ncells); // [W/m2] net heat flux into snow surface layer
    auto dhsdT = create<ArrayD1>("dhsdT", ncells); // derivative of heat flux wrt temperature
    const auto& soitop = nlevsno;
    auto surface_heat_fluxes = [=] (const int& c) {

      const int snotop = nlevsno-snl(c);

      hs_soil(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_soil(c), t_soisno(c, soitop), eflx_sh_soil(c), qflx_ev_soil(c));

      hs_h2osfc(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_soil(c), t_h2osfc(c), eflx_sh_h2osfc(c), qflx_ev_h2osfc(c));

      hs_top(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_lyr(c, snotop), t_grnd(c), eflx_sh_grnd(c), qflx_evap_soi(c));

      hs_top_snow(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_lyr(c, snotop), t_soisno(c, snotop), eflx_sh_snow(c), qflx_ev_snow(c));

      dhsdT(c) = ELM::soil_temp::calc_dhsdT(cgrnd(c), emg(c), t_grnd(c));
    };
    invoke_kernel(surface_heat_fluxes, std::make_tuple(ncells), "surface_heat_fluxes");


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // call pre-solve soil temperature kernels
    // diffusive heat fluxes and matrix factor used in solve
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    auto fn = create<ArrayD2>("fn", ncells, nlevgrnd + nlevsno); // heat diffusion through the layer interface [W/m2]
    //auto fact = create<ArrayD2>("fact", ncells, nlevgrnd + nlevsno); // factors used in computing tridiagonal matrix - needed outside
    auto diffusive_heat_flux = [=] (const int& c) {

      calc_diffusive_heat_flux(snl(c),
          Kokkos::subview(tk, c, Kokkos::ALL),
          Kokkos::subview(t_soisno, c, Kokkos::ALL),
          Kokkos::subview(zsoi, c, Kokkos::ALL),
          Kokkos::subview(fn, c, Kokkos::ALL));

      detail::calc_heat_flux_matrix_factor(snl(c),
          dtime,
          Kokkos::subview(cv, c, Kokkos::ALL),
          Kokkos::subview(dz, c, Kokkos::ALL),
          Kokkos::subview(zsoi, c, Kokkos::ALL),
          Kokkos::subview(zisoi, c, Kokkos::ALL),
          Kokkos::subview(fact, c, Kokkos::ALL));
    };
    invoke_kernel(diffusive_heat_flux, std::make_tuple(ncells), "diffusive_heat_flux");


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // set RHS vector and LHS matrix for temperature solve
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    auto rhs_vector = create<ArrayD2>("rhs_vector", ncells, nlevgrnd + nlevsno + 1); // RHS for soil temp solve
    auto lhs_matrix = create<ArrayD3>("lhs_matrix", ncells, nlevgrnd + nlevsno + 1, ELM::ELMdims::nband); // LHS for soil temp solve
    // these launch their own parallel loops
    set_RHS(dtime, snl, hs_top_snow, dhsdT, hs_soil, frac_sno_eff, t_soisno, fact, fn, sabg_lyr, zsoi,
      tk_h2osfc, t_h2osfc, dz_h2osfc, c_h2osfc, hs_h2osfc, rhs_vector);
    set_LHS(dtime, snl, dz_h2osfc, c_h2osfc, tk_h2osfc, frac_h2osfc, frac_sno_eff, dhsdT, zsoi, fact, tk, lhs_matrix);


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // solve for temperature
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // maximum size of system
    const int N = nlevgrnd + nlevsno + 1;

    auto A = create<ArrayD2>("A", ncells, N - 1); // A(Ai, ..., An-1)
    auto B = create<ArrayD2>("B", ncells, N - 2); // B(Bi, ..., Bn-2)
    auto Z = create<ArrayD2>("Z", ncells, N); // Z(Zi, ..., Zn)

    auto solver = [=] (const int& c) {
      // rhs_vector will contain solution
      solver::PDMA(c, snl, lhs_matrix, A, B, Z, rhs_vector);
    };
    invoke_kernel(solver, std::make_tuple(snl.extent(0)), "soil_temp::solve_temp");


    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // place updated temperature into t_soisno, t_h2osfc
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    auto update_temp = [=] (const int& c) {
      detail::update_temperature(c, snl, frac_h2osfc, rhs_vector, t_h2osfc, t_soisno);
    };
    invoke_kernel(update_temp, std::make_tuple(ncells), "soil_temp::update_temperature");
  }

} // namespace ELM::soil_temp


namespace ELM::soil_temp::detail {

  ACCELERATE
  double calc_lwrad_emit(const double& emg, const double& temp)
  {
    return emg * ELMconst::STEBOL * pow(temp, 4.0);
  }


  ACCELERATE
  double calc_dlwrad_emit(const double& emg, const double& t_grnd)
  {
    return 4.0 * emg * ELMconst::STEBOL * pow(t_grnd, 3.0);
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
    const int top = nlevsno - snl;
    for (int i = 0; i < top; ++i) {
      fact(i) = 0.0;
    }

    // top active layer (snow or top subsurface layer if no snow)
    fact(top) = dtime / cv(top) * dz(top) / (0.5 * (z(top) -
        zi(top) + capr * (z(top+1) - zi(top))));

    // all layers below top active layer
    for (int i = top + 1; i < nlevgrnd + nlevsno; ++i) {
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
    const int top = nlevsno - snl(c);
    
    // zero inactive layers
    for (int i = 0; i < top; ++i)
      tvector(c, i) = 0.0;

    // snow layers
    for (int i = top; i < nlevsno; ++i)
      tvector(c, i) = t_soisno(c, i);

    // surface water
    tvector(c, nlevsno) = t_h2osfc(c);

    // soil
    for (int i = nlevsno; i < nlevsno + nlevgrnd; ++i)
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
    const int top = nlevsno - snl(c);

    // snow layers
    for (int i = top; i < nlevsno; ++i)
      t_soisno(c, i) = tvector(c, i);

    // soil
    for (int i = nlevsno; i < nlevsno + nlevgrnd; ++i)
      t_soisno(c, i) = tvector(c, i + 1);

    // surface water
    // ? first subsurface cell : surface water
    t_h2osfc(c) = (frac_h2osfc(c) == 0.) ? t_soisno(c, nlevsno) : tvector(c, nlevsno);
  }

} // namespace ELM::soil_temp::detail
