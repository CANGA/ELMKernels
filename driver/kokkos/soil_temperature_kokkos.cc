
#include "invoke_kernel.hh"
#include "soil_temperature.h"
#include "soil_temperature_kokkos.hh"

void ELM::kokkos_soil_temperature(ELMStateType& S,
                                  const double& dtime)
{
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call soil_temperature kernels
  // parallel loops are invoked further into the call tree
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

  // ELM dimensions
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;
  using Utils::create;
  using Utils::assign;

  // copy view pointers from ELMState
  auto snl = S.snl;
  auto frac_veg_nosno = S.frac_veg_nosno;
  auto dlrad = S.dlrad;
  auto emg = S.emg;
  auto forc_lwrad = S.forc_lwrad;
  auto htvp = S.htvp;
  auto cgrnd = S.cgrnd;
  auto eflx_sh_soil = S.eflx_sh_soil;
  auto qflx_ev_soil = S.qflx_ev_soil;
  auto eflx_sh_h2osfc = S.eflx_sh_h2osfc;
  auto qflx_ev_h2osfc = S.qflx_ev_h2osfc;
  auto eflx_sh_grnd = S.eflx_sh_grnd;
  auto qflx_evap_soi = S.qflx_evap_soi;
  auto eflx_sh_snow = S.eflx_sh_snow;
  auto qflx_ev_snow = S.qflx_ev_snow;
  auto frac_sno_eff = S.frac_sno_eff;
  auto frac_sno = S.frac_sno;
  auto frac_h2osfc = S.frac_h2osfc;
  auto sabg_snow = S.sabg_snow;
  auto sabg_soil = S.sabg_soil;
  auto sabg_lyr = S.sabg_lyr;
  auto watsat = S.watsat;
  auto sucsat = S.sucsat;
  auto bsw = S.bsw;
  auto tkmg = S.tkmg;
  auto tkdry = S.tkdry;
  auto csol = S.csol;
  auto dz = S.dz;
  auto zsoi = S.zsoi;
  auto zisoi = S.zisoi;
  auto h2osfc = S.h2osfc;
  auto h2osno = S.h2osno;
  auto snow_depth = S.snow_depth;
  auto int_snow = S.int_snow;
  auto t_h2osfc = S.t_h2osfc;
  auto t_grnd = S.t_grnd;
  auto xmf_h2osfc = S.xmf_h2osfc;
  auto xmf = S.xmf;
  auto qflx_h2osfc_to_ice = S.qflx_h2osfc_ice;
  auto eflx_h2osfc_to_snow = S.eflx_h2osfc_snow;
  auto qflx_snofrz = S.qflx_snofrz;
  auto qflx_snow_melt = S.qflx_snow_melt;
  auto qflx_snomelt = S.qflx_snomelt;
  auto eflx_snomelt = S.eflx_snomelt;
  auto imelt = S.imelt;
  auto h2osoi_liq = S.h2osoi_liq;
  auto h2osoi_ice = S.h2osoi_ice;
  auto qflx_snofrz_lyr = S.qflx_snofrz_lyr;
  auto t_soisno = S.t_soisno;
  auto fact = S.fact;

  size_t ncols = snl.extent(0);

  // dummy ltype for now
  ViewI1 ltype("ltype", ncols);
  assign(ltype, 1);


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call soil thermal properties kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

  ViewD2 tk("tk", ncols, nlevgrnd + nlevsno); // thermal conductivity at layer interface
  ViewD2 cv("cv", ncols, nlevgrnd + nlevsno);
  ViewD1 tk_h2osfc("tk_h2osfc", ncols);
  ViewD1 c_h2osfc("c_h2osfc", ncols);
  ViewD1 dz_h2osfc("dz_h2osfc", ncols);
  {
    ViewD2 thk("thk", ncols, nlevgrnd + nlevsno); // thermal conductivity of layer
    auto soil_thermal_props = ELM_LAMBDA (const int& c) {
      soil_thermal::calc_soil_tk(c, ltype(c), h2osoi_liq, h2osoi_ice, t_soisno, dz, watsat, tkmg, tkdry, thk);
      soil_thermal::calc_snow_tk(c, snl(c), frac_sno(c), h2osoi_liq, h2osoi_ice, dz, thk);
      soil_thermal::calc_face_tk(c, snl(c), thk, zsoi, zisoi, tk);
      soil_thermal::calc_soil_heat_capacity(c, ltype(c), snl(c), h2osno(c), watsat, h2osoi_ice, h2osoi_liq, dz, csol, cv);
      soil_thermal::calc_snow_heat_capacity(c, snl(c), frac_sno(c), h2osoi_ice, h2osoi_liq, cv);
      tk_h2osfc(c) = soil_thermal::calc_h2osfc_tk(c, h2osfc(c), thk, zsoi);
      c_h2osfc(c) = ELM::soil_thermal::calc_h2osfc_heat_capacity(snl(c), h2osfc(c), frac_h2osfc(c));
      dz_h2osfc(c) = ELM::soil_thermal::calc_h2osfc_height(snl(c), h2osfc(c), frac_h2osfc(c));
    };
    invoke_kernel(soil_thermal_props, std::make_tuple(ncols), "soil_thermal_props");
  }


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call pre-solve soil temperature kernels
  // surface heat fluxes
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  ViewD1 hs_soil("hs_soil", ncols); // [W/m2] soil heat flux
  //ViewD1 hs_snow("hs_snow", ncols); // [W/m2] snow heat flux - maybe use for coupling?
  ViewD1 hs_h2osfc("hs_h2osfc", ncols); // [W/m2] standing water heat flux
  //ViewD1 hs_top("hs_top", ncols); // [W/m2] net heat flux into surface layer - maybe use for coupling?
  ViewD1 hs_top_snow("hs_top_snow", ncols); // [W/m2] net heat flux into snow surface layer
  ViewD1 dhsdT("dhsdT", ncols); // derivative of heat flux wrt temperature
  const auto& soitop = nlevsno;
  {
    auto surface_heat_fluxes = ELM_LAMBDA (const int& c) {

      const int snotop = nlevsno-snl(c);

      hs_soil(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_soil(c), t_soisno(c, soitop), eflx_sh_soil(c), qflx_ev_soil(c));

      hs_h2osfc(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_soil(c), t_h2osfc(c), eflx_sh_h2osfc(c), qflx_ev_h2osfc(c));

      //hs_top(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
      //    sabg_lyr(c, snotop), t_grnd(c), eflx_sh_grnd(c), qflx_evap_soi(c));

      hs_top_snow(c) = ELM::soil_temp::calc_surface_heat_flux(frac_veg_nosno(c), dlrad(c), emg(c), forc_lwrad(c), htvp(c),
          sabg_lyr(c, snotop), t_soisno(c, snotop), eflx_sh_snow(c), qflx_ev_snow(c));

      dhsdT(c) = ELM::soil_temp::calc_dhsdT(cgrnd(c), emg(c), t_grnd(c));
    };
  invoke_kernel(surface_heat_fluxes, std::make_tuple(ncols), "surface_heat_fluxes");
  }


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call pre-solve soil temperature kernels
  // diffusive heat fluxes and matrix factor used in solve
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  ViewD2 fn("fn", ncols, nlevgrnd + nlevsno); // heat diffusion through the layer interface [W/m2]
  {
    auto diffusive_heat_flux = ELM_LAMBDA (const int& c) {

      ELM::soil_temp::calc_diffusive_heat_flux(snl(c),
          Kokkos::subview(tk, c, Kokkos::ALL),
          Kokkos::subview(t_soisno, c, Kokkos::ALL),
          Kokkos::subview(zsoi, c, Kokkos::ALL),
          Kokkos::subview(fn, c, Kokkos::ALL));

      ELM::soil_temp::calc_heat_flux_matrix_factor(snl(c),
          dtime,
          Kokkos::subview(cv, c, Kokkos::ALL),
          Kokkos::subview(dz, c, Kokkos::ALL),
          Kokkos::subview(zsoi, c, Kokkos::ALL),
          Kokkos::subview(zisoi, c, Kokkos::ALL),
          Kokkos::subview(fact, c, Kokkos::ALL));
    };
    invoke_kernel(diffusive_heat_flux, std::make_tuple(ncols), "diffusive_heat_flux");
  }


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // set RHS vector and LHS matrix for temperature solve
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  ViewD2 rhs_vector("rhs_vector", ncols, nlevgrnd + nlevsno + 1); // RHS for soil temp solve
  ViewD3 lhs_matrix("lhs_matrix", ncols, nlevgrnd + nlevsno + 1, ELM::ELMdims::nband); // LHS for soil temp solve
  // these launch their own parallel loops
  ELM::soil_temp::set_RHS(dtime, snl, hs_top_snow, dhsdT, hs_soil, frac_sno_eff, t_soisno, fact, fn, sabg_lyr, zsoi,
    tk_h2osfc, t_h2osfc, dz_h2osfc, c_h2osfc, hs_h2osfc, rhs_vector);
  ELM::soil_temp::set_LHS(dtime, snl, dz_h2osfc, c_h2osfc, tk_h2osfc, frac_h2osfc, frac_sno_eff, dhsdT, zsoi, fact, tk, lhs_matrix);

    // print matrix
    //#include <iostream>
    //#include <iomanip>
    //auto print_matrix = [] (const double& xx) { return (xx == 0) ? " - " : " X ";};
    //std::cout << "LHS MATRIX:" << std::endl;
    //std::cout << "             0  1  2  3  4" << std::endl;
    //for (size_t i = 0; i < lhs_matrix.extent(1); ++i) {
    //  std::cout << std::setw(2) << i;
    //  std::cout << "          " <<
    //    print_matrix(lhs_matrix(0, i, 0)) <<
    //    print_matrix(lhs_matrix(0, i, 1)) <<
    //    print_matrix(lhs_matrix(0, i, 2)) <<
    //    print_matrix(lhs_matrix(0, i, 3)) <<
    //    print_matrix(lhs_matrix(0, i, 4)) << std::endl;
    //}

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // solve for temperature
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // maximum size of system
  const int N = nlevgrnd + nlevsno + 1;

  ViewD2 A("A", ncols, N - 1); // A(Ai, ..., An-1)
  ViewD2 B("B", ncols, N - 2); // B(Bi, ..., Bn-2)
  ViewD2 Z("Z", ncols, N); // Z(Zi, ..., Zn)

  {
    auto solver = ELM_LAMBDA (const int& c) {
      // rhs_vector will contain solution
      solver::PDMA(c, snl, lhs_matrix, A, B, Z, rhs_vector);
    };
    invoke_kernel(solver, std::make_tuple(snl.extent(0)), "soil_temp::solve_temp");
  }

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // place updated temperature into t_soisno, t_h2osfc
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  {
    auto update_temp = ELM_LAMBDA (const int& c) {
      ELM::soil_temp::update_temperature(c, snl, frac_h2osfc, rhs_vector, t_h2osfc, t_soisno);
    };
    invoke_kernel(update_temp, std::make_tuple(ncols), "soil_temp::update_temperature");
  }

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // phase change kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  //auto fn1 = create<ArrayD2>("fn1", ncols, nlevgrnd+nlevsno);
  {
    auto phase_change = ELM_LAMBDA (const int& c) {

      //  phase_change_correction(snl(c), Kokkos::subview(tk, c, Kokkos::ALL),
      //    Kokkos::subview(t_soisno, c, Kokkos::ALL), Kokkos::subview(zsoi, c, Kokkos::ALL),
      //    Kokkos::subview(fn1, c, Kokkos::ALL));

      ELM::soil_temp::phase_change_h2osfc(snl(c), dtime, frac_sno(c), frac_h2osfc(c), dhsdT(c), c_h2osfc(c),
        fact(c, nlevsno - 1), t_h2osfc(c), h2osfc(c), xmf_h2osfc(c), qflx_h2osfc_to_ice(c),
        eflx_h2osfc_to_snow(c), h2osno(c), int_snow(c), snow_depth(c), h2osoi_ice(c, nlevsno - 1),
        t_soisno(c, nlevsno - 1));

      ELM::soil_temp::phase_change_soisno(snl(c), ltype(c), dtime, dhsdT(c), frac_h2osfc(c), frac_sno_eff(c),
        Kokkos::subview(fact, c, Kokkos::ALL), Kokkos::subview(watsat, c, Kokkos::ALL),
        Kokkos::subview(sucsat, c, Kokkos::ALL), Kokkos::subview(bsw, c, Kokkos::ALL),
        Kokkos::subview(dz, c, Kokkos::ALL), h2osno(c), snow_depth(c), xmf(c), qflx_snofrz(c),
        qflx_snow_melt(c), qflx_snomelt(c), eflx_snomelt(c), Kokkos::subview(imelt, c, Kokkos::ALL),
        Kokkos::subview(qflx_snofrz_lyr, c, Kokkos::ALL), Kokkos::subview(h2osoi_ice, c, Kokkos::ALL),
        Kokkos::subview(h2osoi_liq, c, Kokkos::ALL), Kokkos::subview(t_soisno, c, Kokkos::ALL));
    };
    invoke_kernel(phase_change, std::make_tuple(ncols), "soil_temp::phase_change");
  }


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // update t_grnd
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  {
    auto update_tgrnd = ELM_LAMBDA (const int& c) {
      ELM::soil_temp::update_t_grnd(c, snl, frac_h2osfc, frac_sno_eff, t_h2osfc, t_soisno, t_grnd);
    };
    invoke_kernel(update_tgrnd, std::make_tuple(ncols), "soil_temp::update_t_grnd");
  }

}
