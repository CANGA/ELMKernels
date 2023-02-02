
#include "invoke_kernel.hh"
#include "surface_albedo.h"
#include "snow_snicar.h"
#include "aerosol_physics.h"
#include "snicar_data.h"
#include "albedo_kokkos.hh"


void ELM::kokkos_albedo_snicar(ELMStateType& S)
{
  // get data managers from state
  auto& aerosol_concentrations = *S.aerosol_concentrations.get();
  auto& snicar_data = *S.snicar_data.get();
  auto& pft_data = *S.pft_data.get();

  // local variables
  const int ncols = S.num_columns;
  ViewI2 snw_rds_lcl("snw_rds_lcl", ncols, ELM::ELMdims::nlevsno());
  ViewD2 h2osoi_ice_lcl("h2osoi_ice_lcl", ncols, ELM::ELMdims::nlevsno());
  ViewD2 h2osoi_liq_lcl("h2osoi_liq_lcl", ncols, ELM::ELMdims::nlevsno());
  ViewD2 albout_lcl("albout_lcl", ncols, ELM::ELMdims::numrad_snw());
  ViewD2 flx_slrd_lcl("flx_slrd_lcl", ncols, ELM::ELMdims::numrad_snw());
  ViewD2 flx_slri_lcl("flx_slri_lcl", ncols, ELM::ELMdims::numrad_snw());
  ViewD2 tsai_z("tsai_z", ncols, ELM::ELMdims::nlevcan());

  ViewD2 fabd_sun("fabd_sun", ncols, ELM::ELMdims::numrad());
  ViewD2 fabd_sha("fabd_sha", ncols, ELM::ELMdims::numrad());
  ViewD2 fabi_sun("fabi_sun", ncols, ELM::ELMdims::numrad());
  ViewD2 fabi_sha("fabi_sha", ncols, ELM::ELMdims::numrad());

  ViewD3 flx_abs_lcl("flx_abs_lcl", ncols, ELM::ELMdims::nlevsno()+1, ELM::ELMdims::numrad_snw());
  ViewD3 mss_cnc_aer_in_fdb("mss_cnc_aer_in_fdb", ncols, ELM::ELMdims::nlevsno(), ELM::ELMdims::sno_nbr_aer());
  ViewD3 g_star("g_star", ncols, ELM::ELMdims::numrad_snw(), ELM::ELMdims::nlevsno());
  ViewD3 omega_star("omega_star", ncols, ELM::ELMdims::numrad_snw(), ELM::ELMdims::nlevsno());
  ViewD3 tau_star("tau_star", ncols, ELM::ELMdims::numrad_snw(), ELM::ELMdims::nlevsno());
  ViewD3 flx_absd_snw("flx_absd_snw", ncols, ELM::ELMdims::nlevsno()+1, ELM::ELMdims::numrad());
  ViewD3 flx_absi_snw("flx_absi_snw", ncols, ELM::ELMdims::nlevsno()+1, ELM::ELMdims::numrad());

  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call surface albedo and SNICAR kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto albedo_kernels = ELM_LAMBDA (const int& idx) {

    // thread-local
    int snl_top, snl_btm, flg_nosnl;
    double mu_not;
    // parse pft data for Land.vtype
    ELM::PFTDataAlb alb_pft = pft_data.get_pft_alb(S.vtype(idx));

    ELM::surface_albedo::init_timestep(
        S.Land.urbpoi,
        S.elai(idx),
        Kokkos::subview(aerosol_concentrations.mss_cnc_bcphi, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations.mss_cnc_bcpho, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations.mss_cnc_dst1, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations.mss_cnc_dst2, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations.mss_cnc_dst3, idx, Kokkos::ALL),
        Kokkos::subview(aerosol_concentrations.mss_cnc_dst4, idx, Kokkos::ALL),
        S.vcmaxcintsun(idx),
        S.vcmaxcintsha(idx),
        Kokkos::subview(S.albsod, idx, Kokkos::ALL),
        Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S.albgri, idx, Kokkos::ALL),
        Kokkos::subview(S.albd, idx, Kokkos::ALL),
        Kokkos::subview(S.albi, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd, idx, Kokkos::ALL),
        Kokkos::subview(fabd_sun, idx, Kokkos::ALL),
        Kokkos::subview(fabd_sha, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sun, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sha, idx, Kokkos::ALL),
        Kokkos::subview(S.ftdd, idx, Kokkos::ALL),
        Kokkos::subview(S.ftid, idx, Kokkos::ALL),
        Kokkos::subview(S.ftii, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absdv, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absdn, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absiv, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absin, idx, Kokkos::ALL),
        Kokkos::subview(mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL));

    ELM::surface_albedo::soil_albedo(
        S.Land,
        S.snl(idx),
        S.t_grnd(idx),
        S.coszen(idx),
        Kokkos::subview(S.h2osoi_vol, idx, Kokkos::ALL),
        Kokkos::subview(S.albsat, S.isoicol(idx), Kokkos::ALL),
        Kokkos::subview(S.albdry, S.isoicol(idx), Kokkos::ALL),
        Kokkos::subview(S.albsod, idx, Kokkos::ALL),
        Kokkos::subview(S.albsoi, idx, Kokkos::ALL));

    {
      int flg_slr_in = 1; // direct-beam

      ELM::snow_snicar::init_timestep (
          S.Land.urbpoi,
          flg_slr_in,
          S.coszen(idx),
          S.h2osno(idx),
          S.snl(idx),
          Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
          Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
          Kokkos::subview(S.snw_rds, idx, Kokkos::ALL),
          snl_top,
          snl_btm,
          Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL),
          flg_nosnl,
          Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
          Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
          mu_not,
          Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL));

      ELM::snow_snicar::snow_aerosol_mie_params(
          S.Land.urbpoi,
          flg_slr_in,
          snl_top,
          snl_btm,
          S.coszen(idx),
          S.h2osno(idx),
          Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
          snicar_data.ss_alb_oc1,
          snicar_data.asm_prm_oc1,
          snicar_data.ext_cff_mss_oc1,
          snicar_data.ss_alb_oc2,
          snicar_data.asm_prm_oc2,
          snicar_data.ext_cff_mss_oc2,
          snicar_data.ss_alb_dst1,
          snicar_data.asm_prm_dst1,
          snicar_data.ext_cff_mss_dst1,
          snicar_data.ss_alb_dst2,
          snicar_data.asm_prm_dst2,
          snicar_data.ext_cff_mss_dst2,
          snicar_data.ss_alb_dst3,
          snicar_data.asm_prm_dst3,
          snicar_data.ext_cff_mss_dst3,
          snicar_data.ss_alb_dst4,
          snicar_data.asm_prm_dst4,
          snicar_data.ext_cff_mss_dst4,
          snicar_data.ss_alb_snw_drc,
          snicar_data.asm_prm_snw_drc,
          snicar_data.ext_cff_mss_snw_drc,
          snicar_data.ss_alb_snw_dfs,
          snicar_data.asm_prm_snw_dfs,
          snicar_data.ext_cff_mss_snw_dfs,
          snicar_data.ss_alb_bc1,
          snicar_data.asm_prm_bc1,
          snicar_data.ext_cff_mss_bc1,
          snicar_data.ss_alb_bc2,
          snicar_data.asm_prm_bc2,
          snicar_data.ext_cff_mss_bc2,
          snicar_data.bcenh,
          Kokkos::subview(mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL));

      ELM::snow_snicar::snow_radiative_transfer_solver(
          S.Land.urbpoi,
          flg_slr_in,
          flg_nosnl,
          snl_top,
          snl_btm,
          S.coszen(idx),
          S.h2osno(idx),
          mu_not,
          Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
          Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL));    

      ELM::snow_snicar::snow_albedo_radiation_factor(
          S.Land.urbpoi,
          flg_slr_in,
          snl_top,
          S.coszen(idx),
          mu_not,
          S.h2osno(idx),
          Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
          Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S.albsnd, idx, Kokkos::ALL),
          Kokkos::subview(flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL));
    }


    {
      int flg_slr_in = 2; // diffuse

      ELM::snow_snicar::init_timestep (
          S.Land.urbpoi,
          flg_slr_in,
          S.coszen(idx),
          S.h2osno(idx),
          S.snl(idx),
          Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
          Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
          Kokkos::subview(S.snw_rds, idx, Kokkos::ALL),
          snl_top,
          snl_btm,
          Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL),
          flg_nosnl,
          Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
          Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
          mu_not,
          Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL));

      ELM::snow_snicar::snow_aerosol_mie_params(
          S.Land.urbpoi,
          flg_slr_in,
          snl_top,
          snl_btm,
          S.coszen(idx),
          S.h2osno(idx),
          Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(h2osoi_ice_lcl, idx, Kokkos::ALL),
          Kokkos::subview(h2osoi_liq_lcl, idx, Kokkos::ALL),
          snicar_data.ss_alb_oc1,
          snicar_data.asm_prm_oc1,
          snicar_data.ext_cff_mss_oc1,
          snicar_data.ss_alb_oc2,
          snicar_data.asm_prm_oc2,
          snicar_data.ext_cff_mss_oc2,
          snicar_data.ss_alb_dst1,
          snicar_data.asm_prm_dst1,
          snicar_data.ext_cff_mss_dst1,
          snicar_data.ss_alb_dst2,
          snicar_data.asm_prm_dst2,
          snicar_data.ext_cff_mss_dst2,
          snicar_data.ss_alb_dst3,
          snicar_data.asm_prm_dst3,
          snicar_data.ext_cff_mss_dst3,
          snicar_data.ss_alb_dst4,
          snicar_data.asm_prm_dst4,
          snicar_data.ext_cff_mss_dst4,
          snicar_data.ss_alb_snw_drc,
          snicar_data.asm_prm_snw_drc,
          snicar_data.ext_cff_mss_snw_drc,
          snicar_data.ss_alb_snw_dfs,
          snicar_data.asm_prm_snw_dfs,
          snicar_data.ext_cff_mss_snw_dfs,
          snicar_data.ss_alb_bc1,
          snicar_data.asm_prm_bc1,
          snicar_data.ext_cff_mss_bc1,
          snicar_data.ss_alb_bc2,
          snicar_data.asm_prm_bc2,
          snicar_data.ext_cff_mss_bc2,
          snicar_data.bcenh,
          Kokkos::subview(mss_cnc_aer_in_fdb, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL));

      ELM::snow_snicar::snow_radiative_transfer_solver(
          S.Land.urbpoi,
          flg_slr_in,
          flg_nosnl,
          snl_top,
          snl_btm,
          S.coszen(idx),
          S.h2osno(idx),
          mu_not,
          Kokkos::subview(flx_slrd_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_slri_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
          Kokkos::subview(g_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(omega_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(tau_star, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL));

      ELM::snow_snicar::snow_albedo_radiation_factor(
          S.Land.urbpoi,
          flg_slr_in,
          snl_top,
          S.coszen(idx),
          mu_not,
          S.h2osno(idx),
          Kokkos::subview(snw_rds_lcl, idx, Kokkos::ALL),
          Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
          Kokkos::subview(albout_lcl, idx, Kokkos::ALL),
          Kokkos::subview(flx_abs_lcl, idx, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(S.albsni, idx, Kokkos::ALL),
          Kokkos::subview(flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL));
    }

    ELM::surface_albedo::ground_albedo(
        S.Land.urbpoi,
        S.coszen(idx),
        S.frac_sno(idx),
        Kokkos::subview(S.albsod, idx, Kokkos::ALL),
        Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.albsnd, idx, Kokkos::ALL),
        Kokkos::subview(S.albsni, idx, Kokkos::ALL),
        Kokkos::subview(S.albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S.albgri, idx, Kokkos::ALL));

    ELM::surface_albedo::flux_absorption_factor(
        S.Land,
        S.coszen(idx),
        S.frac_sno(idx),
        Kokkos::subview(S.albsod, idx, Kokkos::ALL),
        Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.albsnd, idx, Kokkos::ALL),
        Kokkos::subview(S.albsni, idx, Kokkos::ALL),
        Kokkos::subview(flx_absd_snw, idx, Kokkos::ALL, Kokkos::ALL),
        Kokkos::subview(flx_absi_snw, idx, Kokkos::ALL, Kokkos::ALL),
        Kokkos::subview(S.flx_absdv, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absdn, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absiv, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absin, idx, Kokkos::ALL));

    ELM::surface_albedo::canopy_layer_lai(
        S.Land.urbpoi,
        S.elai(idx),
        S.esai(idx),
        S.tlai(idx),
        S.tsai(idx),
        S.nrad(idx),
        Kokkos::subview(S.tlai_z, idx, Kokkos::ALL),
        Kokkos::subview(tsai_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd_sha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sha_z, idx, Kokkos::ALL));

    ELM::surface_albedo::two_stream_solver(
        S.Land,
        S.nrad(idx),
        S.coszen(idx),
        S.t_veg(idx),
        S.fwet(idx),
        S.elai(idx),
        S.esai(idx),
        Kokkos::subview(S.tlai_z, idx, Kokkos::ALL),
        Kokkos::subview(tsai_z, idx, Kokkos::ALL),
        Kokkos::subview(S.albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S.albgri, idx, Kokkos::ALL),
        alb_pft,
        S.vcmaxcintsun(idx),
        S.vcmaxcintsha(idx),
        Kokkos::subview(S.albd, idx, Kokkos::ALL),
        Kokkos::subview(S.ftid, idx, Kokkos::ALL),
        Kokkos::subview(S.ftdd, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd, idx, Kokkos::ALL),
        Kokkos::subview(fabd_sun, idx, Kokkos::ALL),
        Kokkos::subview(fabd_sha, idx, Kokkos::ALL),
        Kokkos::subview(S.albi, idx, Kokkos::ALL),
        Kokkos::subview(S.ftii, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sun, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sha, idx, Kokkos::ALL),
        Kokkos::subview(S.fsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd_sha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sha_z, idx, Kokkos::ALL));
  }; // end albedo lambda
  apply_parallel_for(albedo_kernels, "kokkos_albedo_and_snicar", ncols);
}



