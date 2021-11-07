/*! \file SnowSNICAR.h
\brief Functions derived from subroutine SNICAR_AD_RT() in SnowSNICARMod.F90

Calculates snow albedo and per snow layer absorbed flux 

Call sequence:
InitTimestep() -> SnowAerosolMieParams() -> SnowRadiativeTransfer() -> SnowAlbedoRadiationFlux()
Must be called twice, for both direct (flg_slr_in == 1) and diffuse (flg_slr_in == 2) radiation fluxes
*/
#pragma once

#include "ELMConstants.h"

namespace ELM {
namespace SNICAR {

constexpr int numrad_snw = 5; // number of spectral bands used in snow model [nbr]
constexpr int nir_bnd_bgn = 1; // first band index in near-IR spectrum [idx]
constexpr int nir_bnd_end = 4; // ending near-IR band index [idx]
constexpr int sno_nbr_aer = 8; // number of aerosol species in snowpack
constexpr double min_snw = 1.0e-30; // minimum snow mass required for SNICAR RT calculation [kg m-2]
constexpr int ngmax = 8; // gaussian integration index

// idx_*_min and idx_*_max vars have been reduced by one to account for 0/1 indexing - idx_bc*, idx_T*, idx_rhos*
// some variables are statically allocated based on the idx_*_max values
// need to remember to add one back when using idx_*_max for array sizing
constexpr int idx_bc_nclrds_min = 0; // minimum index for BC particle size in optics lookup table
constexpr int idx_bc_nclrds_max = 9; // maximum index for BC particle size in optics lookup table
constexpr int idx_bcint_icerds_min = 0; // minimum index for snow grain size in optics lookup table for within-ice BC
constexpr int idx_bcint_icerds_max = 7; // maximum index for snow grain size in optics lookup table for within-ice BC
constexpr int snw_rds_max_tbl = 1500;// maximum effective radius defined in Mie lookup table [microns]
constexpr int snw_rds_min_tbl = 30; // minimium effective radius defined in Mie lookup table [microns]
constexpr double snw_rds_max = 1500.0; // maximum allowed snow effective radius [microns]
constexpr double snw_rds_refrz  = 1000.0; // effective radius of re-frozen snow [microns]

constexpr double rds_bcint_lcl = 100.0; // effective radius of within-ice BC [nm]
constexpr double rds_bcext_lcl = 100.0; // effective radius of external BC [nm]

constexpr int idx_Mie_snw_mx = 1471; // number of effective radius indices used in Mie lookup table [idx]
constexpr int idx_T_max      = 10; // maxiumum temperature index used in aging lookup table [idx]
constexpr int idx_T_min      = 0; // minimum temperature index used in aging lookup table [idx]
constexpr int idx_Tgrd_max   = 30; // maxiumum temperature gradient index used in aging lookup table [idx]
constexpr int idx_Tgrd_min   = 0; // minimum temperature gradient index used in aging lookup table [idx]
constexpr int idx_rhos_max   = 7; // maxiumum snow density index used in aging lookup table [idx]
constexpr int idx_rhos_min   = 0; // minimum snow density index used in aging lookup table [idx]


// Gaussian integration angle and coefficients for diffuse radiation
constexpr double difgauspt[8] // gaussian angles (radians)
= { 0.9894009,  0.9445750,
    0.8656312,  0.7554044,
    0.6178762,  0.4580168,
    0.2816036,  0.0950125 };
constexpr double difgauswt[8] //gaussian weights
= { 0.0271525,  0.0622535,
    0.0951585,  0.1246290,
    0.1495960,  0.1691565,
    0.1826034,  0.1894506 };

// constants used in algorithm
constexpr double c0      = 0.0;
constexpr double c1      = 1.0;
constexpr double c3      = 3.0;
constexpr double c4      = 4.0;
constexpr double c6      = 6.0;
constexpr double cp01    = 0.01;
constexpr double cp5     = 0.5;
constexpr double cp75    = 0.75;
constexpr double c1p5    = 1.5;
constexpr double trmin   = 0.001;
constexpr double argmax  = 10.0; // maximum argument of exponential

// cconstant coefficients used for SZA parameterization
constexpr double sza_a0 =  0.085730;
constexpr double sza_a1 = -0.630883;
constexpr double sza_a2 =  1.303723;
constexpr double sza_b0 =  1.467291;
constexpr double sza_b1 = -3.338043;
constexpr double sza_b2 =  6.807489;
constexpr double puny   =  1.0e-11;
constexpr double mu_75  =  0.2588; // cosine of 75 degree

constexpr double exp_min = exp(-argmax);

// shorter name for pi
const double pi = ELM_PI;
// always use Delta approximation for snow
constexpr int DELTA = 1;


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
template <class ArrayI1, class ArrayD1, class ArrayD2>
void InitTimestep (const int &urbpoi, const int &flg_slr_in, const double &coszen,
  const double &h2osno, const int &snl, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice,
  const ArrayD1 snw_rds, int &snl_top, int &snl_btm, ArrayD2 flx_abs_lcl, ArrayD2 flx_abs, int &flg_nosnl,
  ArrayD1 h2osoi_ice_lcl, ArrayD1 h2osoi_liq_lcl, ArrayI1 snw_rds_lcl,
  double &mu_not, ArrayD1 flx_slrd_lcl, ArrayD1 flx_slri_lcl);



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
\param[in]ext_cff_mss_snw_dfs[numrad_snw][idx_Mie_snw_mx]                [double] Mie mass extinction coefficients for diffuse ice [frc]
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
template <class ArrayI1, class ArrayD1, class ArrayD2, class ArrayD3>
void SnowAerosolMieParams(const int &urbpoi, const int &flg_slr_in, const int &snl_top, const int &snl_btm,
  const double &coszen, const double &h2osno, const ArrayI1 snw_rds_lcl, 
  const ArrayD1 h2osoi_ice_lcl, const ArrayD1 h2osoi_liq_lcl,
  const ArrayD1 ss_alb_oc1, 
  const ArrayD1 asm_prm_oc1, const ArrayD1 ext_cff_mss_oc1, const ArrayD1 ss_alb_oc2, const ArrayD1 asm_prm_oc2, 
  const ArrayD1 ext_cff_mss_oc2, const ArrayD1 ss_alb_dst1, const ArrayD1 asm_prm_dst1,
  const ArrayD1 ext_cff_mss_dst1, const ArrayD1 ss_alb_dst2, const ArrayD1 asm_prm_dst2, 
  const ArrayD1 ext_cff_mss_dst2, const ArrayD1 ss_alb_dst3, const ArrayD1 asm_prm_dst3, 
  const ArrayD1 ext_cff_mss_dst3, const ArrayD1 ss_alb_dst4, const ArrayD1 asm_prm_dst4, 
  const ArrayD1 ext_cff_mss_dst4, const ArrayD2 ss_alb_snw_drc,
  const ArrayD2 asm_prm_snw_drc, const ArrayD2 ext_cff_mss_snw_drc, const ArrayD2 ss_alb_snw_dfs,
  const ArrayD2 asm_prm_snw_dfs, const ArrayD2 ext_cff_mss_snw_dfs, const ArrayD2 ss_alb_bc1, 
  const ArrayD2 asm_prm_bc1, const ArrayD2 ext_cff_mss_bc1, const ArrayD2 ss_alb_bc2, const ArrayD2 asm_prm_bc2, 
  const ArrayD2 ext_cff_mss_bc2, const ArrayD3 bcenh, const ArrayD2 mss_cnc_aer_in,
  ArrayD2 g_star, ArrayD2 omega_star, ArrayD2 tau_star);




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
template <class ArrayD1, class ArrayD2>
void SnowRadiativeTransfer(const int &urbpoi, const int &flg_slr_in, const int &flg_nosnl, const int &snl_top, const int &snl_btm, 
  const double &coszen, const double &h2osno, const double &mu_not,
  const ArrayD1 flx_slrd_lcl, const ArrayD1 flx_slri_lcl, const ArrayD1 albsoi, const ArrayD2 g_star, 
  const ArrayD2 omega_star, const ArrayD2 tau_star, ArrayD1 albout_lcl, 
  ArrayD2 flx_abs_lcl);





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
template <class ArrayI1, class ArrayD1, class ArrayD2>
void SnowAlbedoRadiationFlux(const bool &urbpoi, const int &flg_slr_in, const int &snl_top, const double &coszen,
  const double &mu_not,
const double &h2osno, 
const ArrayI1 snw_rds_lcl, 
const ArrayD1 albsoi, const ArrayD1 albout_lcl, const ArrayD2 flx_abs_lcl, 
ArrayD1 albout, ArrayD2 flx_abs);

} // namespace ELM
} // namespace SnowSnicar

#include "SnowSNICAR_impl.hh"
