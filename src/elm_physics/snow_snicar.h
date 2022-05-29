/*! \file snow_snicar.h
\brief Functions derived from subroutine SNICAR_AD_RT() in SnowSNICARMod.F90

Calculates snow albedo and per snow layer absorbed flux

Call sequence:
init_timestep() -> snow_aerosol_mie_params() -> snow_radiative_transfer_solver() -> snow_albedo_radiation_factor()
Must be called twice, for both direct (flg_slr_in == 1) and diffuse (flg_slr_in == 2) radiation fluxes
*/
#pragma once

#include "elm_constants.h"
#include <cmath>
#include <stdexcept>

#include "kokkos_includes.hh"

namespace ELM::snow_snicar {


// idx_*_min and idx_*_max vars have been reduced by one to account for 0/1 indexing - idx_bc*, idx_T*, idx_rhos*
// some variables are statically allocated based on the idx_*_max values
// need to remember to add one back when using idx_*_max for array sizing
namespace detail {
static constexpr double min_snw{1.0e-30};      // minimum snow mass required for SNICAR RT calculation [kg m-2]
static constexpr int idx_bc_nclrds_min{0};     // minimum index for BC particle size in optics lookup table
static constexpr int idx_bc_nclrds_max{9};     // maximum index for BC particle size in optics lookup table
static constexpr int idx_bcint_icerds_min{0};  // minimum index for snow grain size in optics lookup table for within-ice BC
static constexpr int idx_bcint_icerds_max{7};  // maximum index for snow grain size in optics lookup table for within-ice BC
static constexpr int idx_Mie_snw_mx{1471};     // number of effective radius indices used in Mie lookup table [idx]
static constexpr int snw_rds_max_tbl{1500};    // maximum effective radius defined in Mie lookup table [microns]
static constexpr int snw_rds_min_tbl{30};      // minimium effective radius defined in Mie lookup table [microns]

static constexpr int idx_T_max{10};        // maxiumum temperature index used in aging lookup table [idx]
static constexpr int idx_Tgrd_max{30};     // maxiumum temperature gradient index used in aging lookup table [idx]
static constexpr int idx_rhos_max{7};      // maxiumum snow density index used in aging lookup table [idx]
static constexpr int idx_T_min{0};         // minimum temperature index used in aging lookup table [idx]
static constexpr int idx_Tgrd_min{0};      // minimum temperature gradient index used in aging lookup table [idx]
static constexpr int idx_rhos_min{0};      // minimum snow density index used in aging lookup table [idx]
} // namespace detail

// constants used in algorithm
namespace alg {
static constexpr double c0{0.0};
static constexpr double c1{1.0};
static constexpr double c3{3.0};
static constexpr double c4{4.0};
static constexpr double c6{6.0};
static constexpr double cp01{0.01};
static constexpr double cp5{0.5};
static constexpr double cp75{0.75};
static constexpr double c1p5{1.5};
static constexpr double trmin{0.001};
} // namespace alg

using ELMdims::nlevsno;
using ELMdims::numrad;
using ELMdims::numrad_snw;
using ELMdims::sno_nbr_aer;

/*! Initialize variables for SNICAR kernels
\param[in]urbpoi                              [bool] true if urban point, false otherwise
\param[in]flg_slr_in                          [int] flag: ==1 for direct-beam incident flux, ==2 for diffuse incident flux
\param[in]coszen                              [double] solar zenith angle factor
\param[in]h2osno                              [double] snow water (mm H2O)
\param[in]snl                                 [int] number of snow layers
\param[in]h2osoi_liq[nlevgrnd+nlevsno]        [double] liquid water content (col,lyr) [kg/m2]
\param[in]h2osoi_ice[nlevgrnd+nlevsno]        [double] ice lens content (col,lyr) [kg/m2]
\param[in]snw_rds[nlevsno]                    [double] snow grain radius (col,lyr) [microns]

\param[out]snl_top                            [int] top snow layer index [idx]
\param[out]snl_btm                            [int] bottom snow layer index [idx]
\param[out]flx_abs_lcl[nlevsno+1][numrad_snw] [double] absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]
\param[out]flx_abs[nlevsno+1][numrad]         [double] absorbed flux in each layer per unit flux incident [frc]
\param[out]flg_nosnl                          [int] flag: =1 if there is snow, but zero snow layers
\param[out]h2osoi_ice_lcl[nlevsno]            [double] liquid water mass [kg/m2]
\param[out]h2osoi_liq_lcl[nlevsno]            [double] ice mass [kg/m2]
\param[out]snw_rds_lcl[nlevsno]               [int] snow effective radius [m^-6]
\param[out]mu_not                             [double] cosine solar zenith angle above the fresnel level
\param[out]flx_slrd_lcl[numrad_snw]           [double] direct beam incident irradiance [W/m2]
\param[out]flx_slri_lcl[numrad_snw]           [double] diffuse incident irradiance [W/m2]
*/
template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
ACCELERATE
void init_timestep(const bool& urbpoi, const int& flg_slr_in, const double& coszen, const double& h2osno, const int& snl,
                   const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 snw_rds, int& snl_top,
                   int& snl_btm, ArrayD2 flx_abs_lcl, ArrayD2 flx_abs, int& flg_nosnl, ArrayD1 h2osoi_ice_lcl,
                   ArrayD1 h2osoi_liq_lcl, ArrayI1 snw_rds_lcl, double& mu_not, ArrayD1 flx_slrd_lcl,
                   ArrayD1 flx_slri_lcl);

/*! Get appropriate variables from input lookup tables and use them to calculate Delta-Eddington parameters
\param[in]urbpoi                                                         [bool] true if urban point, false otherwise
\param[in]flg_slr_in                                                     [int] flag: ==1 for direct-beam incident flux, ==2 for diffuse incident flux
\param[in]snl_top                                                        [int] top snow layer index [idx]
\param[in]snl_btm                                                        [int] bottom snow layer index [idx]
\param[in]coszen                                                         [double] solar zenith angle factor
\param[in]h2osno                                                         [double] snow water (mm H2O)
\param[in]snw_rds_lcl[nlevsno]                                           [int] snow effective radius [m^-6]
\param[in]h2osoi_ice_lcl[nlevsno]                                        [double] liquid water mass [kg/m2]
\param[in]h2osoi_liq_lcl[nlevsno]                                        [double] ice mass [kg/m2]
\param[in]ss_alb_oc1[numrad_snw]                                         [double] Mie single scatter albedos for hydrophillic OC [frc]
\param[in]asm_prm_oc1[numrad_snw]                                        [double] Mie asymmetry parameters for hydrophillic OC [frc]
\param[in]ext_cff_mss_oc1[numrad_snw]                                    [double] Mie mass extinction coefficients for hydrophillic OC [frc]
\param[in]ss_alb_oc2[numrad_snw]                                         [double] Mie single scatter albedos for hydrophobic OC [frc]
\param[in]asm_prm_oc2[numrad_snw]                                        [double] Mie asymmetry parameters for hydrophobic OC [frc]
\param[in]ext_cff_mss_oc2[numrad_snw]                                    [double] Mie mass extinction coefficients for hydrophobic OC [frc]
\param[in]ss_alb_dst1[numrad_snw]                                        [double] Mie single scatter albedos for dust species 1 [frc]
\param[in]asm_prm_dst1[numrad_snw]                                       [double] Mie asymmetry parameters for dust species 1 [frc]
\param[in]ext_cff_mss_dst1[numrad_snw]                                   [double] Mie mass extinction coefficients for dust species 1 [frc]
\param[in]ss_alb_dst2[numrad_snw]                                        [double] Mie single scatter albedos for dust species 2 [frc]
\param[in]asm_prm_dst2[numrad_snw]                                       [double] Mie asymmetry parameters for dust species 2 [frc]
\param[in]ext_cff_mss_dst2[numrad_snw]                                   [double] Mie mass extinction coefficients for dust species 2 [frc]
\param[in]ss_alb_dst3[numrad_snw]                                        [double] Mie single scatter albedos for dust species 3 [frc]
\param[in]asm_prm_dst3[numrad_snw]                                       [double] Mie asymmetry parameters for dust species 3 [frc]
\param[in]ext_cff_mss_dst3[numrad_snw]                                   [double] Mie mass extinction coefficients for dust species 3 [frc]
\param[in]ss_alb_dst4[numrad_snw]                                        [double] Mie single scatter albedos for dust species 4 [frc]
\param[in]asm_prm_dst4[numrad_snw]                                       [double] Mie asymmetry parameters for dust species 4 [frc]
\param[in]ext_cff_mss_dst4[numrad_snw]                                   [double] Mie mass extinction coefficients for dust species 4 [frc
\param[in]ss_alb_snw_drc[numrad_snw][idx_Mie_snw_mx]                     [double] Mie single scatter albedos for direct-beam ice [frc]
\param[in]asm_prm_snw_drc[numrad_snw][idx_Mie_snw_mx]                    [double] Mie asymmetry parameters for direct-beam ice [frc]
\param[in]ext_cff_mss_snw_drc[numrad_snw][idx_Mie_snw_mx]                [double] Mie mass extinction coefficients for direct-beam ice [frc]
\param[in]ss_alb_snw_dfs[numrad_snw][idx_Mie_snw_mx]                     [double] Mie single scatter albedos for diffuse ice [frc]
\param[in]asm_prm_snw_dfs[numrad_snw][idx_Mie_snw_mx]                    [double] Mie asymmetry parameters for diffuse ice [frc]
\param[in]ext_cff_mss_snw_dfs[numrad_snw][idx_Mie_snw_mx]                [double] Mie mass extinction coefficients for diffuse ice [frc
\param[in]ss_alb_bc1[idx_bc_nclrds_max+1][numrad_snw]                    [double] Mie single scatter albedos for within-ice BC [frc]
\param[in]asm_prm_bc1[idx_bc_nclrds_max+1][numrad_snw]                   [double] Mie asymmetry parameters for within-ice BC [frc]
\param[in]ext_cff_mss_bc1[idx_bc_nclrds_max+1][numrad_snw]               [double] Mie mass extinction coefficients for within-ice BC [frc]
\param[in]ss_alb_bc2[idx_bc_nclrds_max+1][numrad_snw]                    [double] Mie single scatter albedos for external BC [frc]
\param[in]asm_prm_bc2[idx_bc_nclrds_max+1][numrad_snw]                   [double] Mie asymmetry parameters for external BC [frc]
\param[in]ext_cff_mss_bc2[idx_bc_nclrds_max+1][numrad_snw]               [double] Mie mass extinction coefficients for external BC [frc]
\param[in]mss_cnc_aer_in[nlevsno][sno_nbr_aer]                           [double] mass concentration of all aerosol species [kg/kg]
\param[in]bcenh[idx_bcint_icerds_max+1][idx_bc_nclrds_max+1][numrad_snw] [double] Absorption enhancement factors for within-ice BC

\param[out]g_star[numrad_snw][nlevsno]                                   [double] transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
\param[out]omega_star[numrad_snw][nlevsno]                               [double] transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer [frc]
\param[out]tau_star[numrad_snw][nlevsno]                                 [double] transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer [-]
*/
template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3, typename SubviewD1, typename SubviewD2>
ACCELERATE
void snow_aerosol_mie_params(const bool& urbpoi, const int& flg_slr_in, const int& snl_top, const int& snl_btm,
                             const double& coszen, const double& h2osno, const ArrayI1 snw_rds_lcl,
                             const SubviewD1 h2osoi_ice_lcl, const SubviewD1 h2osoi_liq_lcl, const ArrayD1 ss_alb_oc1,
                             const ArrayD1 asm_prm_oc1, const ArrayD1 ext_cff_mss_oc1, const ArrayD1 ss_alb_oc2,
                             const ArrayD1 asm_prm_oc2, const ArrayD1 ext_cff_mss_oc2, const ArrayD1 ss_alb_dst1,
                             const ArrayD1 asm_prm_dst1, const ArrayD1 ext_cff_mss_dst1, const ArrayD1 ss_alb_dst2,
                             const ArrayD1 asm_prm_dst2, const ArrayD1 ext_cff_mss_dst2, const ArrayD1 ss_alb_dst3,
                             const ArrayD1 asm_prm_dst3, const ArrayD1 ext_cff_mss_dst3, const ArrayD1 ss_alb_dst4,
                             const ArrayD1 asm_prm_dst4, const ArrayD1 ext_cff_mss_dst4,
                             const ArrayD2 ss_alb_snw_drc, const ArrayD2 asm_prm_snw_drc,
                             const ArrayD2 ext_cff_mss_snw_drc, const ArrayD2 ss_alb_snw_dfs,
                             const ArrayD2 asm_prm_snw_dfs, const ArrayD2 ext_cff_mss_snw_dfs,
                             const ArrayD2 ss_alb_bc1, const ArrayD2 asm_prm_bc1, const ArrayD2 ext_cff_mss_bc1,
                             const ArrayD2 ss_alb_bc2, const ArrayD2 asm_prm_bc2, const ArrayD2 ext_cff_mss_bc2,
                             const ArrayD3 bcenh, const SubviewD2 mss_cnc_aer_in, SubviewD2 g_star, SubviewD2 omega_star,
                             SubviewD2 tau_star);

/*! Snow radiative transfer solver

\param[in]urbpoi                              [bool] true if urban point, false otherwise
\param[in]flg_slr_in                          [int] flag: ==1 for direct-beam incident flux, ==2 for diffuse incident flux 
\param[in]flg_nosnl                           [int] flag: =1 if there is snow, but zero snow layers
\param[in]snl_top                             [int] top snow layer index [idx]
\param[in]snl_btm                             [int] bottom snow layer index [idx]
\param[in]coszen                              [double] solar zenith angle factor
\param[in]h2osno                              [double] snow water (mm H2O)
\param[in]mu_not                              [double] cosine solar zenith angle above the fresnel level
\param[in]flx_slrd_lcl[numrad_snw]            [double] direct beam incident irradiance [W/m2]
\param[in]flx_slri_lcl[numrad_snw]            [double] diffuse incident irradiance [W/m2]
\param[in]albsoi[numrad]                      [double] albedo of surface underlying snow [frc]
\param[in]g_star[numrad_snw][nlevsno]         [double] transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
\param[in]omega_star[numrad_snw][nlevsno]     [double] transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer [frc]
\param[in]tau_star[numrad_snw][nlevsno]       [double] transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer [-]
\param[out]albout_lcl[numrad_snw]             [double] snow albedo by band [frc]
\param[out]flx_abs_lcl[nlevsno+1][numrad_snw] [double] absorbed flux per unit incident flux at top of snowpack [frc]
*/
template <typename ArrayD1, typename ArrayD2>
ACCELERATE
void snow_radiative_transfer_solver(const bool& urbpoi, const int& flg_slr_in, const int& flg_nosnl, const int& snl_top,
                                    const int& snl_btm, const double& coszen, const double& h2osno,
                                    const double& mu_not, const ArrayD1 flx_slrd_lcl, const ArrayD1 flx_slri_lcl,
                                    const ArrayD1 albsoi, const ArrayD2 g_star, const ArrayD2 omega_star,
                                    const ArrayD2 tau_star, ArrayD1 albout_lcl, ArrayD2 flx_abs_lcl);

/*! Calculates final snow albedo and snow absorbed radiative flux
\param[in]urbpoi                             [bool] true if urban point, false otherwise
\param[in]flg_slr_in                         [int] flag: ==1 for direct-beam incident flux, ==2 for diffuse incident flux
\param[in]snl_top                            [int] top snow layer index [idx]
\param[in]coszen                             [double] solar zenith angle factor
\param[in]mu_not                             [double] cosine solar zenith angle above the fresnel level
\param[in]h2osno                             [double] snow water (mm H2O)
\param[in]snw_rds_lcl[nlevsno]               [int] snow effective radius [m^-6]
\param[in]albsoi[numrad]                     [double] albedo of surface underlying snow [frc]
\param[in]albout_lcl[numrad_snw]             [double] snow albedo by band [frc]
\param[in]flx_abs_lcl[nlevsno+1][numrad_snw] [double] absorbed flux per unit incident flux at top of snowpack [frc]
\param[out]albout[numrad]                    [double] snow albedo, averaged into 2 bands (=0 if no sun or no snow) [frc]
\param[out]flx_abs[nlevsno+1][numrad]        [double] absorbed flux in each layer per unit flux incident, averaged into 2 bands  [frc]
*/
template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
ACCELERATE
void snow_albedo_radiation_factor(const bool& urbpoi, const int& flg_slr_in, const int& snl_top, const double& coszen,
                                  const double& mu_not, const double& h2osno, const ArrayI1 snw_rds_lcl,
                                  const ArrayD1 albsoi, const ArrayD1 albout_lcl, const ArrayD2 flx_abs_lcl,
                                  ArrayD1 albout, ArrayD2 flx_abs);

} // namespace ELM::snow_snicar

#include "snow_snicar_impl.hh"
