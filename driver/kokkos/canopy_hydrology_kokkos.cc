
#include "invoke_kernel.hh"
#include "atm_data.h"
#include "canopy_hydrology.h"
#include "canopy_hydrology_kokkos.hh"

void ELM::kokkos_canopy_hydrology(ELMStateType& S,
                                 const double& model_dt_secs)
{
  const int ncols = S.num_columns;
  ViewD1 qflx_candrip("qflx_candrip", ncols);
  ViewD1 qflx_through_snow("qflx_through_snow", ncols);
  ViewD1 qflx_through_rain("qflx_through_rain", ncols);
  ViewD1 fracsnow("fracsnow", ncols);
  ViewD1 fracrain("fracrain", ncols);

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
        S.forc_rain(idx),
        S.forc_snow(idx),
        S.dewmx,
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
        S.forc_rain(idx),
        S.forc_snow(idx),
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

    ELM::canopy_hydrology::snow_init(
        S.Land,
        model_dt_secs,
        S.do_capsnow(idx),
        S.oldfflag,
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
  apply_parallel_for(canhydro_kernels, "kokkos_canopy_hydrology", ncols);
}

void ELM::kokkos_frac_wet(ELMStateType& S)
{
  auto fracwet_kernel = ELM_LAMBDA (const int& idx) {
    ELM::canopy_hydrology::fraction_wet(
        S.Land,
        S.frac_veg_nosno(idx),
        S.dewmx,
        S.elai(idx),
        S.esai(idx),
        S.h2ocan(idx),
        S.fwet(idx),
        S.fdry(idx));
  }; // end fracwet lambda
  apply_parallel_for(fracwet_kernel, "kokkos_canhydro_fracwet_kernel", S.num_columns);
}
