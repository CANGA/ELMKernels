
#include "soil_temperature.h"
#include "soil_temperature_kokkos.hh"

void ELM::kokkos_soil_temperature(const std::shared_ptr<ELMStateType>& S,
                                  const double& dtime)
{
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call soil_temperature kernels
  // parallel loops are invoked further into the call tree
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  ELM::soil_temp::solve_temperature<ViewD3>(
    dtime,
    S->snl,
    S->frac_veg_nosno,
    S->dlrad,
    S->emg,
    S->forc_lwrad,
    S->htvp,
    S->cgrnd,
    S->eflx_sh_soil,
    S->qflx_ev_soil,
    S->eflx_sh_h2osfc,
    S->qflx_ev_h2osfc,
    S->eflx_sh_grnd,
    S->qflx_evap_soi,
    S->eflx_sh_snow,
    S->qflx_ev_snow,
    S->frac_sno_eff,
    S->frac_sno,
    S->frac_h2osfc,
    S->sabg_snow,
    S->sabg_soil,
    S->sabg_lyr,
    S->watsat,
    S->sucsat,
    S->bsw,
    S->tkmg,
    S->tkdry,
    S->csol,
    S->dz,
    S->zsoi,
    S->zisoi,
    S->h2osfc,
    S->h2osno,
    S->snow_depth,
    S->int_snow,
    S->t_h2osfc,
    S->t_grnd,
    S->xmf_h2osfc_dummy,
    S->xmf_dummy,
    S->qflx_h2osfc_ice_dummy,
    S->eflx_h2osfc_snow_dummy,
    S->qflx_snofrz,
    S->qflx_snow_melt,
    S->qflx_snomelt,
    S->eflx_snomelt,
    S->imelt,
    S->h2osoi_liq,
    S->h2osoi_ice,
    S->qflx_snofrz_lyr,
    S->t_soisno,
    S->fact);
}
