
#include "surface_albedo.h"
#include "snow_snicar.h"
#include "invoke_kernel.hh"

#include "albedo_kokkos.hh"


void ELM::kokkos_albedo_snicar(const std::shared_ptr<ELMStateType>& S,
                               const std::shared_ptr<ELM::AerosolConcentrations<ViewD2>>& aerosol_concentrations,
                               const std::shared_ptr<ELM::SnicarData<ViewD1, ViewD2, ViewD3>>& snicar_data,
                               const std::shared_ptr<ELM::PFTData<ViewD1, ViewD2>>& pft_data)
{

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call surface albedo and SNICAR kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto albedo_kernels = ELM_LAMBDA (const int& idx) {
    // parse pft data for Land.vtype
    ELM::PFTDataAlb alb_pft = pft_data->get_pft_alb(S->vtype(idx));

    ELM::surface_albedo::init_timestep(
        S->Land.urbpoi,
        S->elai(idx),
        Kokkos::subview(aerosol_concentrations->mss_cnc_bcphi, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations->mss_cnc_bcpho, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations->mss_cnc_dst1, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations->mss_cnc_dst2, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations->mss_cnc_dst3, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations->mss_cnc_dst4, idx, Kokkos::ALL),
        S->vcmaxcintsun(idx),
        S->vcmaxcintsha(idx),
        Kokkos::subview(S->albsod, idx, Kokkos::ALL),
        Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S->albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S->albgri, idx, Kokkos::ALL),
        Kokkos::subview(S->albd, idx, Kokkos::ALL),
        Kokkos::subview(S->albi, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sun, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sha, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sun, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sha, idx, Kokkos::ALL),
        Kokkos::subview(S->ftdd, idx, Kokkos::ALL),
        Kokkos::subview(S->ftid, idx, Kokkos::ALL),
        Kokkos::subview(S->ftii, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absdv, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absdn, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absiv, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absin, idx, Kokkos::ALL),
        Kokkos::subview(S->mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL));

    ELM::surface_albedo::soil_albedo(
        S->Land,
        S->snl(idx),
        S->t_grnd(idx),
        S->coszen(idx),
        Kokkos::subview(S->h2osoi_vol, idx, Kokkos::ALL),
        Kokkos::subview(S->albsat, S->isoicol(idx), Kokkos::ALL),
        Kokkos::subview(S->albdry, S->isoicol(idx), Kokkos::ALL),
        Kokkos::subview(S->albsod, idx, Kokkos::ALL),
        Kokkos::subview(S->albsoi, idx, Kokkos::ALL));

    {
      int flg_slr_in = 1; // direct-beam

      ELM::snow_snicar::init_timestep (
          S->Land.urbpoi,
          flg_slr_in,
          S->coszen(idx),
          S->h2osno(idx),
          S->snl(idx),
          Kokkos::subview(S->h2osoi_liq, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_ice, idx, Kokkos::ALL),
          Kokkos::subview(S->snw_rds, idx, Kokkos::ALL),
          S->snl_top(idx),
          S->snl_btm(idx),
          Kokkos::subview(S->flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL),
          S->flg_nosnl(idx),
          Kokkos::subview(S->h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_liq_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->snw_rds_lcl, idx, Kokkos::ALL),
          S->mu_not(idx),
          Kokkos::subview(S->flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_slri_lcl, idx, Kokkos::ALL));

      ELM::snow_snicar::snow_aerosol_mie_params(
          S->Land.urbpoi,
          flg_slr_in,
          S->snl_top(idx),
          S->snl_btm(idx),
          S->coszen(idx),
          S->h2osno(idx),
          Kokkos::subview(S->snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_liq_lcl, idx, Kokkos::ALL),
          snicar_data->ss_alb_oc1,
          snicar_data->asm_prm_oc1,
          snicar_data->ext_cff_mss_oc1,
          snicar_data->ss_alb_oc2,
          snicar_data->asm_prm_oc2,
          snicar_data->ext_cff_mss_oc2,
          snicar_data->ss_alb_dst1,
          snicar_data->asm_prm_dst1,
          snicar_data->ext_cff_mss_dst1,
          snicar_data->ss_alb_dst2,
          snicar_data->asm_prm_dst2,
          snicar_data->ext_cff_mss_dst2,
          snicar_data->ss_alb_dst3,
          snicar_data->asm_prm_dst3,
          snicar_data->ext_cff_mss_dst3,
          snicar_data->ss_alb_dst4,
          snicar_data->asm_prm_dst4,
          snicar_data->ext_cff_mss_dst4,
          snicar_data->ss_alb_snw_drc,
          snicar_data->asm_prm_snw_drc,
          snicar_data->ext_cff_mss_snw_drc,
          snicar_data->ss_alb_snw_dfs,
          snicar_data->asm_prm_snw_dfs,
          snicar_data->ext_cff_mss_snw_dfs,
          snicar_data->ss_alb_bc1,
          snicar_data->asm_prm_bc1,
          snicar_data->ext_cff_mss_bc1,
          snicar_data->ss_alb_bc2,
          snicar_data->asm_prm_bc2,
          snicar_data->ext_cff_mss_bc2,
          snicar_data->bcenh,
          Kokkos::subview(S->mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->tau_star, idx, Kokkos::ALL, Kokkos::ALL));

      ELM::snow_snicar::snow_radiative_transfer_solver(
          S->Land.urbpoi,
          flg_slr_in,
          S->flg_nosnl(idx),
          S->snl_top(idx),
          S->snl_btm(idx),
          S->coszen(idx),
          S->h2osno(idx),
          S->mu_not(idx),
          Kokkos::subview(S->flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_slri_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
          Kokkos::subview(S->g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->tau_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL));    

      ELM::snow_snicar::snow_albedo_radiation_factor(
          S->Land.urbpoi,
          flg_slr_in,
          S->snl_top(idx),
          S->coszen(idx),
          S->mu_not(idx),
          S->h2osno(idx),
          Kokkos::subview(S->snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
          Kokkos::subview(S->albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->albsnd, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL));
    }


    {
      int flg_slr_in = 2; // diffuse

      ELM::snow_snicar::init_timestep (
          S->Land.urbpoi,
          flg_slr_in,
          S->coszen(idx),
          S->h2osno(idx),
          S->snl(idx),
          Kokkos::subview(S->h2osoi_liq, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_ice, idx, Kokkos::ALL),
          Kokkos::subview(S->snw_rds, idx, Kokkos::ALL),
          S->snl_top(idx),
          S->snl_btm(idx),
          Kokkos::subview(S->flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL),
          S->flg_nosnl(idx),
          Kokkos::subview(S->h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_liq_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->snw_rds_lcl, idx, Kokkos::ALL),
          S->mu_not(idx),
          Kokkos::subview(S->flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_slri_lcl, idx, Kokkos::ALL));

      ELM::snow_snicar::snow_aerosol_mie_params(
          S->Land.urbpoi,
          flg_slr_in,
          S->snl_top(idx),
          S->snl_btm(idx),
          S->coszen(idx),
          S->h2osno(idx),
          Kokkos::subview(S->snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->h2osoi_liq_lcl, idx, Kokkos::ALL),
          snicar_data->ss_alb_oc1,
          snicar_data->asm_prm_oc1,
          snicar_data->ext_cff_mss_oc1,
          snicar_data->ss_alb_oc2,
          snicar_data->asm_prm_oc2,
          snicar_data->ext_cff_mss_oc2,
          snicar_data->ss_alb_dst1,
          snicar_data->asm_prm_dst1,
          snicar_data->ext_cff_mss_dst1,
          snicar_data->ss_alb_dst2,
          snicar_data->asm_prm_dst2,
          snicar_data->ext_cff_mss_dst2,
          snicar_data->ss_alb_dst3,
          snicar_data->asm_prm_dst3,
          snicar_data->ext_cff_mss_dst3,
          snicar_data->ss_alb_dst4,
          snicar_data->asm_prm_dst4,
          snicar_data->ext_cff_mss_dst4,
          snicar_data->ss_alb_snw_drc,
          snicar_data->asm_prm_snw_drc,
          snicar_data->ext_cff_mss_snw_drc,
          snicar_data->ss_alb_snw_dfs,
          snicar_data->asm_prm_snw_dfs,
          snicar_data->ext_cff_mss_snw_dfs,
          snicar_data->ss_alb_bc1,
          snicar_data->asm_prm_bc1,
          snicar_data->ext_cff_mss_bc1,
          snicar_data->ss_alb_bc2,
          snicar_data->asm_prm_bc2,
          snicar_data->ext_cff_mss_bc2,
          snicar_data->bcenh,
          Kokkos::subview(S->mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->tau_star, idx, Kokkos::ALL, Kokkos::ALL));

      ELM::snow_snicar::snow_radiative_transfer_solver(
          S->Land.urbpoi,
          flg_slr_in,
          S->flg_nosnl(idx),
          S->snl_top(idx),
          S->snl_btm(idx),
          S->coszen(idx),
          S->h2osno(idx),
          S->mu_not(idx),
          Kokkos::subview(S->flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_slri_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
          Kokkos::subview(S->g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->tau_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL));

      ELM::snow_snicar::snow_albedo_radiation_factor(
          S->Land.urbpoi,
          flg_slr_in,
          S->snl_top(idx),
          S->coszen(idx),
          S->mu_not(idx),
          S->h2osno(idx),
          Kokkos::subview(S->snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
          Kokkos::subview(S->albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S->albsni, idx, Kokkos::ALL),
          Kokkos::subview(S->flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL));
    }

    ELM::surface_albedo::ground_albedo(
        S->Land.urbpoi,
        S->coszen(idx),
        S->frac_sno(idx),
        Kokkos::subview(S->albsod, idx, Kokkos::ALL),
        Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S->albsnd, idx, Kokkos::ALL),
        Kokkos::subview(S->albsni, idx, Kokkos::ALL),
        Kokkos::subview(S->albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S->albgri, idx, Kokkos::ALL));

    ELM::surface_albedo::flux_absorption_factor(
        S->Land,
        S->coszen(idx),
        S->frac_sno(idx),
        Kokkos::subview(S->albsod, idx, Kokkos::ALL),
        Kokkos::subview(S->albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S->albsnd, idx, Kokkos::ALL),
        Kokkos::subview(S->albsni, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL),
        Kokkos::subview(S->flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL),
        Kokkos::subview(S->flx_absdv, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absdn, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absiv, idx, Kokkos::ALL),
        Kokkos::subview(S->flx_absin, idx, Kokkos::ALL));

    ELM::surface_albedo::canopy_layer_lai(
        S->Land.urbpoi,
        S->elai(idx),
        S->esai(idx),
        S->tlai(idx),
        S->tsai(idx),
        S->nrad(idx),
        S->ncan(idx),
        Kokkos::subview(S->tlai_z, idx, Kokkos::ALL),
        Kokkos::subview(S->tsai_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sha_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sha_z, idx, Kokkos::ALL));

    ELM::surface_albedo::two_stream_solver(
        S->Land,
        S->nrad(idx),
        S->coszen(idx),
        S->t_veg(idx),
        S->fwet(idx),
        S->elai(idx),
        S->esai(idx),
        Kokkos::subview(S->tlai_z, idx, Kokkos::ALL),
        Kokkos::subview(S->tsai_z, idx, Kokkos::ALL),
        Kokkos::subview(S->albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S->albgri, idx, Kokkos::ALL),
        alb_pft,
        S->vcmaxcintsun(idx),
        S->vcmaxcintsha(idx),
        Kokkos::subview(S->albd, idx, Kokkos::ALL),
        Kokkos::subview(S->ftid, idx, Kokkos::ALL),
        Kokkos::subview(S->ftdd, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sun, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sha, idx, Kokkos::ALL),
        Kokkos::subview(S->albi, idx, Kokkos::ALL),
        Kokkos::subview(S->ftii, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sun, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sha, idx, Kokkos::ALL),
        Kokkos::subview(S->fsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabd_sha_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S->fabi_sha_z, idx, Kokkos::ALL));
  }; // end albedo lambda
  invoke_kernel(albedo_kernels, std::make_tuple(S->snl.extent(0)), "kokkos_albedo_and_snicar");
}



