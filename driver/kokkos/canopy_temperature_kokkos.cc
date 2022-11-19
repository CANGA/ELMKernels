
#include "invoke_kernel.hh"
#include "canopy_temperature.h"
#include "canopy_temperature_kokkos.hh"

void ELM::kokkos_canopy_temperature(ELMStateType& S)
{
  // local views
  size_t ncols = S.snl.extent(0);
  ViewD1 qred("qred", ncols); // soil surface relative humidity
  ViewD1 hr("hr", ncols); // relative humidity
  ViewD1 soilalpha("soilalpha", ncols); // evaporation resistance -- not currently used
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call canopy_temperature kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto cantemp_kernels = ELM_LAMBDA (const int& idx) {

    ELM::canopy_temperature::old_ground_temp(
        S.Land,
        S.t_h2osfc(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        S.t_h2osfc_bef(idx),
        Kokkos::subview(S.tssbef, idx, Kokkos::ALL));

    ELM::canopy_temperature::ground_temp(
        S.Land,
        S.snl(idx),
        S.frac_sno_eff(idx),
        S.frac_h2osfc(idx),
        S.t_h2osfc(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        S.t_grnd(idx));

    ELM::canopy_temperature::calc_soilalpha(
        S.Land,
        S.frac_sno(idx),
        S.frac_h2osfc(idx),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.watsat, idx, Kokkos::ALL),
        Kokkos::subview(S.sucsat, idx, Kokkos::ALL),
        Kokkos::subview(S.bsw, idx, Kokkos::ALL),
        Kokkos::subview(S.watdry, idx, Kokkos::ALL),
        Kokkos::subview(S.watopt, idx, Kokkos::ALL),
        qred(idx), hr(idx),
        soilalpha(idx));

    ELM::canopy_temperature::calc_soilbeta(
        S.Land,
        S.frac_sno(idx),
        S.frac_h2osfc(idx),
        Kokkos::subview(S.watsat, idx, Kokkos::ALL),
        Kokkos::subview(S.watfc, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        S.soilbeta(idx));

    ELM::canopy_temperature::humidities(
        S.Land,
        S.snl(idx),
        S.forc_qbot(idx),
        S.forc_pbot(idx),
        S.t_h2osfc(idx),
        S.t_grnd(idx),
        S.frac_sno(idx),
        S.frac_sno_eff(idx),
        S.frac_h2osfc(idx),
        qred(idx),
        hr(idx),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        S.qg_snow(idx),
        S.qg_soil(idx),
        S.qg(idx),
        S.qg_h2osfc(idx),
        S.dqgdT(idx));

    ELM::canopy_temperature::ground_properties(
        S.Land,
        S.snl(idx),
        S.frac_sno(idx),
        S.forc_thbot(idx),
        S.forc_qbot(idx),
        S.elai(idx),
        S.esai(idx),
        S.htop(idx),
        S.pft_data->displar,
        S.pft_data->z0mr,
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        S.emg(idx),
        S.emv(idx),
        S.htvp(idx),
        S.z0mg(idx),
        S.z0hg(idx),
        S.z0qg(idx),
        S.z0mv(idx),
        S.z0hv(idx),
        S.z0qv(idx),
        S.thv(idx),
        S.z0m(idx),
        S.displa(idx));

    ELM::canopy_temperature::forcing_height(
        S.Land,
        S.veg_active(idx),
        S.frac_veg_nosno(idx),
        S.z0m(idx),
        S.z0mg(idx),
        S.forc_tbot(idx),
        S.displa(idx),
        S.forc_hgt_u_patch(idx),
        S.forc_hgt_t_patch(idx),
        S.forc_hgt_q_patch(idx),
        S.thm(idx));

    ELM::canopy_temperature::init_energy_fluxes(
        S.Land,
        S.eflx_sh_tot(idx),
        S.eflx_lh_tot(idx),
        S.eflx_sh_veg(idx),
        S.qflx_evap_tot(idx),
        S.qflx_evap_veg(idx),
        S.qflx_tran_veg(idx));
  }; // end cantemp lambda
  invoke_kernel(cantemp_kernels, std::make_tuple(S.snl.extent(0)), "kokkos_canopy_temperature");
}
