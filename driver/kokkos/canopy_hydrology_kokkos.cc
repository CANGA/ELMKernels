
#include "invoke_kernel.hh"
#include "canopy_hydrology.h"

#include "canopy_hydrology_kokkos.hh"

void ELM::kokkos_canopy_hydrology(ELMStateType& S, AtmDataManager<ViewD1, ViewD2, AtmForcType::PREC>& forc_PREC,
                                 const double& model_dt_secs, const Utils::Date& time_plus_half_dt_secs)
{
  constexpr double dewmx{0.1};
  constexpr int oldfflag{1};

  size_t ncells = S.snl.extent(0);

  // get forc_rain and forc_snow
  ViewD1 forc_rain("forc_rain", ncells);
  ViewD1 forc_snow("forc_snow", ncells);
  forc_PREC.get_atm_forcing(model_dt_secs/86400.0, time_plus_half_dt_secs, S.forc_tbot, forc_rain, forc_snow);

  ViewD1 qflx_candrip("qflx_candrip", ncells);
  ViewD1 qflx_through_snow("qflx_through_snow", ncells);
  ViewD1 qflx_through_rain("qflx_through_rain", ncells);
  ViewD1 fracsnow("fracsnow", ncells);
  ViewD1 fracrain("fracrain", ncells);

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call canopy_hydrology kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto canhydro_kernels = ELM_LAMBDA (const int& idx) {
    // local vars - these need to be thread local in parallel runs
    double qflx_irrig = 0.0; // hardwired here

    ELM::canopy_hydrology::interception(
        S.Land,
        S.frac_veg_nosno(idx),
        forc_rain(idx),
        forc_snow(idx),
        dewmx,
        S.elai(idx),
        S.esai(idx),
        model_dt_secs,
        S.h2ocan(idx),
        qflx_candrip(idx),
        qflx_through_snow(idx),
        qflx_through_rain(idx),
        fracsnow(idx),
        fracrain(idx));

    ELM::canopy_hydrology::ground_flux(
        S.Land,
        S.do_capsnow(idx),
        S.frac_veg_nosno(idx),
        forc_rain(idx),
        forc_snow(idx),
        qflx_irrig,
        qflx_candrip(idx),
        qflx_through_snow(idx),
        qflx_through_rain(idx),
        fracsnow(idx),
        fracrain(idx),
        S.qflx_snwcp_liq(idx),
        S.qflx_snwcp_ice(idx),
        S.qflx_snow_grnd(idx),
        S.qflx_rain_grnd(idx));

    ELM::canopy_hydrology::fraction_wet(
        S.Land,
        S.frac_veg_nosno(idx),
        dewmx,
        S.elai(idx),
        S.esai(idx),
        S.h2ocan(idx),
        S.fwet(idx),
        S.fdry(idx));

    ELM::canopy_hydrology::snow_init(
        S.Land,
        model_dt_secs,
        S.do_capsnow(idx),
        oldfflag,
        S.forc_tbot(idx),
        S.t_grnd(idx),
        S.qflx_snow_grnd(idx),
        S.qflx_snow_melt(idx),
        S.n_melt(idx),
        S.snow_depth(idx),
        S.h2osno(idx),
        S.int_snow(idx),
        Kokkos::subview(S.swe_old, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S.t_soisno, idx, Kokkos::ALL),
        Kokkos::subview(S.frac_iceold, idx, Kokkos::ALL),
        S.snl(idx),
        Kokkos::subview(S.dz, idx, Kokkos::ALL),
        Kokkos::subview(S.zsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.zisoi, idx, Kokkos::ALL),
        Kokkos::subview(S.snw_rds, idx, Kokkos::ALL),
        S.frac_sno_eff(idx),
        S.frac_sno(idx));

    ELM::canopy_hydrology::fraction_h2osfc(
        S.Land,
        S.micro_sigma(idx),
        S.h2osno(idx),
        S.h2osfc(idx),
        Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
        S.frac_sno(idx),
        S.frac_sno_eff(idx),
        S.frac_h2osfc(idx));
  }; // end canhydro lambda
  invoke_kernel(canhydro_kernels, std::make_tuple(S.snl.extent(0)), "kokkos_canopy_hydrology");
}
