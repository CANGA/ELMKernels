
#include <cmath>

//#ifdef MODAL_AER evals to true
//logical,  public, parameter :: DO_SNO_OC =    .false.   ! parameter to include organic carbon (OC)
//                                                          ! in snowpack radiative calculations
//logical,  public, parameter :: DO_SNO_AER =   .true.    ! parameter to include aerosols in snowpack radiative calculations
//flg_snw_ice == 1
//num_snowc filter is just if (snl > 0) {}
namespace ELM {
namespace SNICAR {

//  subroutine SnowAge_grain(bounds, &
//       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
//       waterflux_vars, waterstate_vars, temperature_vars)
    //
    // !DESCRIPTION:
    // Updates the snow effective grain size (radius). 
    // Contributions to grain size evolution are from:
    //   1. vapor redistribution (dry snow) 
    //   2. liquid water redistribution (wet snow)
    //   3. re-freezing of liquid water
    // 
    // Vapor redistribution: Method is to retrieve 3 best-bit parameters that
    // depend on snow temperature, temperature gradient, and density,
    // that are derived from the microphysical model described in: 
    // Flanner and Zender (2006), Linking snowpack microphysics and albedo
    // evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834. 
    // The parametric equation has the form: 
    // dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where:
    //   r is the effective radius,
    //   tau and kappa are best-fit parameters,
    //   drdt_0 is the initial rate of change of effective radius, and
    //   dr_fresh is the difference between the current and fresh snow states 
    //  (r_current - r_fresh).
    //
    // Liquid water redistribution: Apply the grain growth function from:
    //   Brun, E. (1989), Investigation of wet-snow metamorphism in respect of 
    //   liquid-water content, Annals of Glaciology, 13, 22-26.
    //   There are two parameters that describe the grain growth rate as 
    //   a function of snow liquid water content (LWC). The "LWC=0" parameter
    //   is zeroed here because we are accounting for dry snowing with a 
    //   different representation
    //
    // Re-freezing of liquid water: Assume that re-frozen liquid water clumps
    //   into an arbitrarily large effective grain size (snw_rds_refrz). 
    //   The phenomenon is observed (Grenfell), but so far unquantified, as far as 
    //   I am aware.
void SnowAge_grain() {
  if (snl > 0) {

    snl_btm = nlevsno - 1;
    snl_top = nlevsno - snl;

    for (int i = snl_top; i <= snl_btm; ++i) { cdz[i] = frac_sno * dz[i]; }

    // loop over snow layers
    for (int i = snl_top; i <= snl_btm; ++i) {
      //
      //**********  1. DRY SNOW AGING  ***********
      //
      h2osno_lyr = h2osoi_liq[i] + h2osoi_ice[i];

      // temperature gradient
      if (i == snl_top) { 
        // top layer
        t_snotop = t_soisno[snl_top];
        t_snobtm = (t_soisno[i+1]*dz[i] + t_soisno[i]*dz[i+1]) / (dz[i]+dz[i+1]);
      } else {
        t_snotop = (t_soisno[i-1]*dz[i] + t_soisno[i]*dz[i-1]) / (dz[i]+dz[i-1]);
        t_snobtm = (t_soisno[i+1]*dz[i] + t_soisno[i]*dz[i+1]) / (dz[i]+dz[i+1]);
      }

      dTdz[i] = abs((t_snotop - t_snobtm) / cdz[i]);

      // snow density
      rhos = (h2osoi_liq[i] + h2osoi_ice[i]) / cdz[i];

      // make sure rhos doesn't drop below 50 (see rhos_idx below)
      rhos = std::max(50.0, rhos);

      // best-fit table indecies
      T_idx    = round((t_soisno[i] - 223) / 5) + 1;
      Tgrd_idx = round(dTdz[i] / 10) + 1;
      rhos_idx = round((rhos-50) / 50) + 1;

      // boundary check:
      if (T_idx < idx_T_min) { T_idx = idx_T_min; }
      if (T_idx > idx_T_max) { T_idx = idx_T_max; }
      if (Tgrd_idx < idx_Tgrd_min) { Tgrd_idx = idx_Tgrd_min; }
      if (Tgrd_idx > idx_Tgrd_max) { Tgrd_idx = idx_Tgrd_max; }
      if (rhos_idx < idx_rhos_min) { rhos_idx = idx_rhos_min; }
      if (rhos_idx > idx_rhos_max) { rhos_idx = idx_rhos_max; }

      // best-fit parameters
      bst_tau   = snowage_tau[rhos_idx][Tgrd_idx][T_idx];
      bst_kappa = snowage_kappa[rhos_idx][Tgrd_idx][T_idx];     
      bst_drdt0 = snowage_drdt0[rhos_idx][Tgrd_idx][T_idx];


      // change in snow effective radius, using best-fit parameters
      // added checks suggested by mgf. --HW 10/15/2015 
      dr_fresh = snw_rds[i] - snw_rds_min;
      if (abs(dr_fresh) < 1.0e-8) {
         dr_fresh = 0.0;
      } else if (dr_fresh < 0.0) {
        throw std::runtime_error("ELM ERROR: SnowAge dr_fresh < 0.0.");
      }

      dr = (bst_drdt0 * pow(bst_tau / (dr_fresh+bst_tau), 1.0 / bst_kappa)) * (dtime/3600.0);
            
      //
      //**********  2. WET SNOW AGING  ***********
      //
      // We are assuming wet and dry evolution occur simultaneously, and 
      // the contributions from both can be summed. 
      // This is justified by setting the linear offset constant C1_liq_Brun89 to zero [Brun, 1989]

      // liquid water faction
      frc_liq = std::min(0.1, (h2osoi_liq[i] / (h2osoi_liq[i] + h2osoi_ice[i])));

      // dr_wet = 1E6_r8*(dtime*(C1_liq_Brun89 + C2_liq_Brun89*(frc_liq**(3))) / (4*SHR_CONST_PI*(snw_rds(c_idx,i)/1E6)**(2)))
      // simplified, units of microns:
      dr_wet = 1.0e18 * (dtime * (C2_liq_Brun89 * pow(frc_liq, 3.0)) / (4.0 * ELM_PI * pow(snw_rds, 2.0)));
      dr += dr_wet;

      //
      //**********  3. SNOWAGE SCALING (TURNED OFF BY DEFAULT)  *************
      //
      // Multiply rate of change of effective radius by some constant, xdrdt
      if (flg_snoage_scl) { dr = dr*xdrdt; }

      //
      //**********  4. INCREMENT EFFECTIVE RADIUS, ACCOUNTING FOR:  ***********
      //               DRY AGING
      //               WET AGING
      //               FRESH SNOW
      //               RE-FREEZING
      //
      // new snowfall [kg/m2]
      if (do_capsnow) {
        newsnow = std::max(0.0, (qflx_snwcp_ice * dtime));
      } else {
        newsnow = std::max(0.0, (qflx_snow_grnd_col * dtime));
      }

      // snow that has re-frozen [kg/m2]
      refrzsnow = std::max(0.0, (qflx_snofrz_lyr[i] * dtime));

      // fraction of layer mass that is re-frozen
      frc_refrz = refrzsnow / h2osno_lyr;

      // fraction of layer mass that is new snow
      if (i == snl_top) {
        frc_newsnow = newsnow / h2osno_lyr;
      } else {
        frc_newsnow = 0.0;
      }

      if ((frc_refrz + frc_newsnow) > 1.0) {
        frc_refrz = frc_refrz / (frc_refrz + frc_newsnow);
        frc_newsnow = 1.0 - frc_refrz;
        frc_oldsnow = 0.0;
      } else {
        frc_oldsnow = 1.0 - frc_refrz - frc_newsnow;
      }

      // mass-weighted mean of fresh snow, old snow, and re-frozen snow effective radius
      snw_rds(c_idx,i) = (snw_rds(c_idx,i)+dr)*frc_oldsnow + snw_rds_min*frc_newsnow + snw_rds_refrz*frc_refrz
      //
      //**********  5. CHECK BOUNDARIES   ***********
      //
      // boundary check
      if (snw_rds[i] < snw_rds_min) {
         snw_rds[i] = snw_rds_min;
      }

      if (snw_rds[i] > snw_rds_max) {
         snw_rds[i] = snw_rds_max;
      }

      // set top layer variables for history files
      if (i == snl_top) {
        snot_top = t_soisno[i];
        dTdz_top = dTdz[i];
        snw_rds_top = snw_rds[i];
        sno_liq_top = h2osoi_liq[i] / (h2osoi_liq[i]+h2osoi_ice[i]);
      }
    } // for i = snl_top .. snl_btm
  } // if snl > 0

  // Special case: snow on ground, but not enough to have defined a snow layer:
  //   set snw_rds to fresh snow grain size:
  if (snl == 0) {
    if (h2osno > 0.0) { snw_rds[nlevsno-1] = snw_rds_min; }
  }
} // SnowAge_grain



idx_bcint_icerds                  ! index of ice effective radius for optical properties lookup table
idx_bcint_nclrds                  ! index of within-ice BC effective radius for optical properties lookup 
idx_bcext_nclrds                  ! index of external BC effective radius for optical properties lookup ta
enh_fct                           ! extinction/absorption enhancement factor for within-ice BC
tmp1                              ! temporary variable




/*

urbpoi                             [bool]    true if urban point, false otherwise
flg_slr_in                         [int]     flag: ==1 for direct-beam incident flux, ==2 for diffuse incident flux
coszen                             [double]  solar zenith angle factor 
h2osno                         [double]   snow water (mm H2O)
snl                            [int]      number of snow layers
h2osoi_liq[nlevgrnd+nlevsno]   [double]  liquid water content (col,lyr) [kg/m2]
h2osoi_ice[nlevgrnd+nlevsno]   [double]  ice lens content (col,lyr) [kg/m2]    
snw_rds[nlevsno]               [double] snow grain radius (col,lyr) [microns] 



snl_top                            [int]     top snow layer index [idx]
snl_btm                            [int]      bottom snow layer index [idx]
flx_abs_lcl[nlevsno+1][numrad_snw] [double]  absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]
flg_nosnl                          [int]     flag: =1 if there is snow, but zero snow layers,
h2osoi_ice_lcl[nlevsno]            [double] liquid water mass [kg/m2]
h2osoi_liq_lcl[nlevsno]            [double] ice mass [kg/m2]
snw_rds_lcl[nlevsno]               [double]  effective radius [m^-6]
rds_bcint_lcl[nlevsno]             [double]    effective radius of within-ice BC [nm]
rds_bcext_lcl[nlevsno]             [double]    effective radius of external BC [nm]
mu_not                             [double] cosine solar zenith angle above the fresnel level
flx_slrd_lcl[numrad_snw]    [double] direct beam incident irradiance [W/m2] (set to 1)
flx_slri_lcl[numrad_snw]    [double] diffuse incident irradiance [W/m2] (set to 1)



*/
const void SNICARInitTimestep (const int & urbpoi, const int & flg_slr_in, const double &coszen,
  const double &h2osno, const int & snl, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice,
  const ArrayD1 snw_rds, int &snl_top, int &snl_btm, ArrayD2 flx_abs_lcl, int &flg_nosnl,
  ArrayD1 h2osoi_ice_lcl, ArrayD1 h2osoi_liq_lcl, ArrayD1 snw_rds_lcl, ArrayD1 rds_bcint_lcl,
  ArrayD1 rds_bcext_lcl, double &mu_not, flx_slrd_lcl, ArrayD1 flx_slri_lcl) {
  
  if (!Land.urbpoi) {

    // Zero absorbed radiative fluxes:
    for (int i = 0; i <= nlevsno; ++i) {
      for (int ib = 0; i < numrad; ++i) {
        flx_abs[i][ib] = 0.0;
      }
      for (int ib = 0; i < numrad_snw; ++i) {
        flx_abs_lcl[i][ib] = 0.0;
      }
    }

    // Qualifier for computing snow RT:
    //  1) sunlight from atmosphere model
    //  2) minimum amount of snow on ground.
    //     Otherwise, set snow albedo to zero
    if ((coszen > 0.0) && (h2osno > min_snw)) {

      // If there is snow, but zero snow layers, we must create a layer locally.
      // This layer is presumed to have the fresh snow effective radius.
      int snl_lcl;
      if (snl == 0) {
        flg_nosnl = 1;
        snl_lcl = 1;
        h2osoi_ice_lcl[nlevsno-1] = h2osno;
        h2osoi_liq_lcl[nlevsno-1] = 0.0;
        snw_rds_lcl = round(snw_rds_min);
      } else {
        flg_nosnl = 0;
        snl_lcl = snl;
        for (int i = 0; i < nlevsno; ++i) {
          h2osoi_liq_lcl[i] = h2osoi_liq[i];
          h2osoi_ice_lcl[i] = h2osoi_ice[i];
          snw_rds_lcl[i] = snw_rds[i];
        }
      }

      snl_btm = nlevsno - 1; // index of bottom snow layer
      snl_top = nlevsno - snl_lcl; // index of top snow layer

      // Assume fixed BC effective radii of 100nm. This is close to
      // the effective radius of 95nm (number median radius of
      // 40nm) assumed for freshly-emitted BC in MAM.  Future
      // implementations may prognose the BC effective radius in
      // snow.
      for (int i = 0; i < nlevsno; ++i) {
        rds_bcint_lcl[i] = 100.0;
        rds_bcext_lcl[i] = 100.0;
      }

      // Error check for snow grain size:
      for (int i = snl_top; i <= snl_btm; ++i) {
        if ((snw_rds_lcl[i] < snw_rds_min_tbl) || (snw_rds_lcl[i] > snw_rds_max_tbl)) {
          throw std::runtime_error("ELM ERROR: SNICAR snow grain radius out of bounds.");
        }
      }

      // mu_not is cosine solar zenith angle above the fresnel level; make
      // sure mu_not is large enough for stable and meaningful radiation
      // solution: .01 is like sun just touching horizon with its lower edge
      // equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
      mu_not = std::max(coszen, cp01);

      // Set direct or diffuse incident irradiance to 1
      // (This has to be within the bnd loop because mu_not is adjusted in rare cases)
      if (flg_slr_in == 1) {
        for (int bnd_idx = 0; bnd_idx < numrad_snw; ++bnd_idx) {
          flx_slrd_lcl[bnd_idx] = 1.0 / (mu_not * pi); // this corresponds to incident irradiance of 1.0
          flx_slri_lcl[bnd_idx] = 0.0;
        }
      } else if (flg_slr_in == 2) {
        for (int bnd_idx = 0; bnd_idx < numrad_snw; ++bnd_idx) {
          flx_slrd_lcl[bnd_idx] = 0.0;
          flx_slri_lcl[bnd_idx] = 1.0;
        }
      } else {
        throw std::runtime_error("ELM ERROR: SNICAR solar waveband flag out of bounds - only 2 bands supported.");
      }
    }
  }
}





/*



outputs:
!!--!! make these aware of bnd_idx - will need persistence between this call and the next
g_star[numrad_snw][nlevsno]          [double]
omega_star[numrad_snw][nlevsno]      [double]
tau_star[numrad_snw][nlevsno]        [double]


*/
void SnowAerosolMieParams() {

  // Define local Mie parameters based on snow grain size and aerosol species,
  //  retrieved from a lookup table.
  if (!Land.urbpoi) {
    if ((coszen > 0.0) && (h2osno > min_snw)) {

      // Set local aerosol array
      for (int i = 0; i < nlevsno; ++i) {
        for (int j = 0; j < sno_nbr_aer; ++j) { 
          mss_cnc_aer_lcl[i][j] = mss_cnc_aer_in[i][j];
        }
      }

      for (int bnd_idx = 0; bnd_idx < numrad_snw; ++bnd_idx) {

        if ( (numrad_snw == 5) && ((bnd_idx == 4) || (bnd_idx == 3)) ) {
          for (int i = 0; i < nlevsno; ++i) {
            for (int j = 0; j < sno_nbr_aer; ++j) { 
              mss_cnc_aer_lcl[i][j] = 0.0;
            }
          }
        }

        if (flg_slr_in == 1) {
          for (int i = snl_top; i <= snl_btm; ++i) {
            rds_idx = snw_rds_lcl[i] - snw_rds_min_tbl + 1;
            // snow optical properties (direct radiation)
            ss_alb_snw_lcl[i] = ss_alb_snw_drc[rds_idx][bnd_idx];
            asm_prm_snw_lcl[i] = asm_prm_snw_drc[rds_idx][bnd_idx];
            ext_cff_mss_snw_lcl[i] = ext_cff_mss_snw_drc[rds_idx][bnd_idx];
          }
        } else if (flg_slr_in == 2) {
          for (int i = snl_top; i <= snl_btm; ++i) {
            rds_idx = snw_rds_lcl[i] - snw_rds_min_tbl + 1
            // snow optical properties (diffuse radiation)
            ss_alb_snw_lcl[i] = ss_alb_snw_dfs[rds_idx][bnd_idx];
            asm_prm_snw_lcl[i] = asm_prm_snw_dfs[rds_idx][bnd_idx];
            ext_cff_mss_snw_lcl[i] = ext_cff_mss_snw_dfs[rds_idx][bnd_idx];
          }
        }

        //H. Wang
        // aerosol species 3 optical properties
        ss_alb_aer_lcl[2] = ss_alb_oc1[bnd_idx];
        asm_prm_aer_lcl[2] = asm_prm_oc1[bnd_idx];
        ext_cff_mss_aer_lcl[2] = ext_cff_mss_oc1[bnd_idx];

        //  aerosol species 4 optical properties
        ss_alb_aer_lcl[3] = ss_alb_oc2[bnd_idx];
        asm_prm_aer_lcl[3] = asm_prm_oc2[bnd_idx];
        ext_cff_mss_aer_lcl[3] = ext_cff_mss_oc2[bnd_idx];

        // aerosol species 5 optical properties
        ss_alb_aer_lcl[4] = ss_alb_dst1[bnd_idx];
        asm_prm_aer_lcl[4] = asm_prm_dst1[bnd_idx];
        ext_cff_mss_aer_lcl[4] = ext_cff_mss_dst1[bnd_idx];

        // aerosol species 6 optical properties
        ss_alb_aer_lcl[5] = ss_alb_dst2[bnd_idx];
        asm_prm_aer_lcl[5] = asm_prm_dst2[bnd_idx];
        ext_cff_mss_aer_lcl[5] = ext_cff_mss_dst2[bnd_idx];

        // aerosol species 7 optical properties
        ss_alb_aer_lcl[6] = ss_alb_dst3[bnd_idx];
        asm_prm_aer_lcl[6] = asm_prm_dst3[bnd_idx];
        ext_cff_mss_aer_lcl[6] = ext_cff_mss_dst3[bnd_idx];

        // aerosol species 8 optical properties
        ss_alb_aer_lcl[7] = ss_alb_dst4[bnd_idx];
        asm_prm_aer_lcl[7] = asm_prm_dst4[bnd_idx];
        ext_cff_mss_aer_lcl[7] = ext_cff_mss_dst4[bnd_idx];


        // 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
        // 2. optical Depths (tau_snw, tau_aer)
        // 3. weighted Mie properties (tau, omega, g)

        // Weighted Mie parameters of each layer
        for (int i = snl_top; i <= snl_btm; ++i) {
          // mgf++ within-ice and external BC optical properties
          // Lookup table indices for BC optical properties,
          // dependent on snow grain size and BC particle
          // size.

          // valid for 25 < snw_rds < 1625 um:
          if (snw_rds_lcl[i] < 125) {
             tmp1 = snw_rds_lcl[i] / 50;
             idx_bcint_icerds = round(tmp1);
          } else if (snw_rds_lcl[i] < 175) {
             idx_bcint_icerds = 2;
          } else {
             tmp1 = (snw_rds_lcl[i] / 250) + 2;
             idx_bcint_icerds = round(tmp1);
          }

          // valid for 25 < bc_rds < 525 nm
          idx_bcint_nclrds = round(rds_bcint_lcl[i] / 50);
          idx_bcext_nclrds = round(rds_bcext_lcl[i] / 50);

          // check bounds:
          if (idx_bcint_icerds < idx_bcint_icerds_min) idx_bcint_icerds = idx_bcint_icerds_min;
          if (idx_bcint_icerds > idx_bcint_icerds_max) idx_bcint_icerds = idx_bcint_icerds_max;
          if (idx_bcint_nclrds < idx_bc_nclrds_min) idx_bcint_nclrds = idx_bc_nclrds_min;
          if (idx_bcint_nclrds > idx_bc_nclrds_max) idx_bcint_nclrds = idx_bc_nclrds_max;
          if (idx_bcext_nclrds < idx_bc_nclrds_min) idx_bcext_nclrds = idx_bc_nclrds_min;
          if (idx_bcext_nclrds > idx_bc_nclrds_max) idx_bcext_nclrds = idx_bc_nclrds_max;

          // retrieve absorption enhancement factor for within-ice BC
          enh_fct = bcenh[bnd_idx][idx_bcint_nclrds][idx_bcint_icerds];

          // get BC optical properties (moved from above)
          // aerosol species 1 optical properties (within-ice BC)
          ss_alb_aer_lcl[0] = ss_alb_bc1[bnd_idx][idx_bcint_nclrds];
          asm_prm_aer_lcl[0] = asm_prm_bc1[bnd_idx][idx_bcint_nclrds];
          ext_cff_mss_aer_lcl[0] = ext_cff_mss_bc1[bnd_idx][idx_bcint_nclrds] * enh_fct;

          // aerosol species 2 optical properties (external BC)
          ss_alb_aer_lcl[1] = ss_alb_bc2[bnd_idx][idx_bcext_nclrds];
          asm_prm_aer_lcl[1] = asm_prm_bc2[bnd_idx][idx_bcext_nclrds];
          ext_cff_mss_aer_lcl[1] = ext_cff_mss_bc2[bnd_idx][idx_bcext_nclrds];

          L_snw[i] = h2osoi_ice_lcl[i] + h2osoi_liq_lcl[i];
          tau_snw[i] = L_snw[i] * ext_cff_mss_snw_lcl[i];

          for (int j = 0; j < sno_nbr_aer; ++j) {
            L_aer[i][j] = L_snw[i] * mss_cnc_aer_lcl[i][j];
            tau_aer[i][j] = L_aer[i][j] * ext_cff_mss_aer_lcl[j];
          }

          tau_sum   = 0.0;
          omega_sum = 0.0;
          g_sum     = 0.0;

          for (int j = 0; j < sno_nbr_aer; ++j) {
            tau_sum += tau_aer[i][j];
            omega_sum += (tau_aer[i][j] * ss_alb_aer_lcl[j]);
            g_sum += (tau_aer[i][j] * ss_alb_aer_lcl[j] * asm_prm_aer_lcl[j]);
          }

          tau[i] = tau_sum + tau_snw[i];
          omega[i] = (1.0 / tau[i]) * (omega_sum + (ss_alb_snw_lcl[i] * tau_snw[i]));
          g[i] = (1.0 / (tau[i] * omega[i])) * (g_sum + (asm_prm_snw_lcl[i] * ss_alb_snw_lcl[i] * tau_snw[i]));

        } // end Weighted Mie snl loop

        // DELTA transformations, if requested
        if (DELTA == 1) {
          for (int i = snl_top; i <= snl_btm; ++i) {
            g_star[bnd_idx][i] = g[i] / (1.0 + g[i]);
            omega_star[bnd_idx][i] = ((1.0 - pow(g[i], 2.0)) * omega[i]) / (1.0 - (omega[i] * pow(g[i], 2.0)));
            tau_star[bnd_idx][i] = (1.0 - (omega[i] * pow(g[i], 2.0))) * tau[i];
          }
        } else {
          for (int i = snl_top; i <= snl_btm; ++i) {
            g_star[bnd_idx][i] = g[i];
            omega_star[bnd_idx][i] = omega[i];
            tau_star[bnd_idx][i] = tau[i];
          }
        }
      } // end bnd_idx waveband loop
    } // end if coszen > 0 && h2osno > min_snow 
  } // end !urbpoi
}




/*

inputs:
snl_top         [int]      top snow layer index [idx]
snl_btm         [int]      bottom snow layer index [idx]
snl_btm_itf     [int]      index of bottom snow layer interface (nlevsno) [idx]
coszen          [double]   solar zenith angle factor 
h2osno            [double]   snow water (mm H2O)




*/
void SnowRadiativeTransfer() {
  // Begin radiative transfer solver
  // Given input vertical profiles of optical properties, evaluate the
  // monochromatic Delta-Eddington adding-doubling solution

  // note that trndir, trntdr, trndif, rupdir, rupdif, rdndif
  // are variables at the layer interface,
  // for snow with layers rangeing from snl_top to snl_btm
  // there are snl_top to snl_btm+1 layer interface


  /*        interface mapping
  original ELM         this kernel
  
  grid  interface   interface  grid
   --   -4                  0  ---
  |-4|                         |0|
   --   -3                  1  ---
  |-3|                         |1|
   --   -2                  2  ---
  |-2|                         |2|
   --   -1                  3  ---
  |-1|                         |3|
   --    0                  4  ---  [nlevsno - 1]
  | 0|                         |4|  
  -----  1 ground interface 5 ----- [nlevsno]
  | 1|                         |5|  
   --    2                  6  ---  [nlevsno + 1]  */
  

  if (!Land.urbpoi) {
    if ((coszen > 0.0) && (h2osno > min_snw)) {
      // local interface reflect/transmit vars
      double trndir[nlevsno+1];
      double trntdr[nlevsno+1];
      double trndif[nlevsno+1];
      double rupdir[nlevsno+1];
      double rupdif[nlevsno+1];
      double rdndif[nlevsno+1];
      // local layer reflect/transmit vars
      double rdir[nlevsno];
      double rdif_a[nlevsno];
      double rdif_b[nlevsno];
      double tdir[nlevsno];
      double tdif_a[nlevsno];
      double tdif_b[nlevsno];
      double trnlay[nlevsno];
      // net absorbed radiative energy (lyr) [W/m^2]
      double F_abs[nlevsno];

      const int snl_btm_itf = nlevsno; // index of ground/snow interface (same as snl_btm + 1)

      for (int bnd_idx = 0; bnd_idx < numrad_snw; ++bnd_idx) {

        for (int i = snl_top; i <= snl_btm_itf; ++i) {
          trndir[i] = c0;
          trntdr[i] = c0;
          trndif[i] = c0;
          rupdir[i] = c0;
          rupdif[i] = c0;
          rdndif[i] = c0;
        }

        // initialize top interface of top layer
        trndir[snl_top] = c1;
        trntdr[snl_top] = c1;
        trndif[snl_top] = c1;
        rdndif[snl_top] = c0;

        // begin main level loop
        // for layer interfaces except for the very bottom
        for (int i = snl_top; i <= snl_btm; ++i) {
          // initialize all layer apparent optical properties to 0
          rdir[i] = c0;
          rdif_a[i] = c0;
          rdif_b[i] = c0;
          tdir[i] = c0;
          tdif_a[i] = c0;
          tdif_b[i] = c0;
          trnlay[i] = c0;

          // compute next layer Delta-eddington solution only if total transmission
          // of radiation to the interface just above the layer exceeds trmin.

          if (trntdr[i] > trmin ) {

            // calculation over layers with penetrating radiation

            // delta-transformed single-scattering properties
            // of this layer
            const double ts = tau_star[bnd_idx][i];
            const double ws = omega_star[bnd_idx][i];
            const double gs = g_star[bnd_idx][i];

            // Delta-Eddington solution expressions
            // n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
            // u(w,gg,e)        = c1p5*(c1 - w*gg)/e
            // el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
            const double lm = std::sqrt(c3 * (c1 - ws) * (c1 - ws * gs)); // lm = el(ws,gs)
            const double ue = c1p5 * (c1 - ws * gs) / lm; //ue = u(ws,gs,lm)
            const double extins = std::max(exp_min, exp(-lm * ts));
            const double ne = ((ue + c1) * (ue + c1) / extins) - ((ue - c1) * (ue - c1) * extins); // ne = n(ue,extins)

            // first calculation of rdif, tdif using Delta-Eddington formulas
            // rdif_a(k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
            rdif_a[i] = (pow(ue, 2.0)-c1) * (c1/extins - extins) / ne;
            tdif_a[i] = c4 * ue / ne;

            // evaluate rdir,tdir for direct beam
            trnlay[i] = std::max(exp_min, exp(-ts/mu_not));

            // Delta-Eddington solution expressions
            // alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
            // agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
            // alp = alpha(ws,mu_not,gs,lm)
            // gam = agamm(ws,mu_not,gs,lm)
            double alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not));
            double gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not));
            double apg = alp + gam;
            double amg = alp - gam;

            rdir[i] = apg * rdif_a[i] +  amg * (tdif_a[i] * trnlay[i] - c1);
            tdir[i] = apg * tdif_a[i] + (amg *  rdif_a[i] - apg + c1) * trnlay[i];

            // recalculate rdif,tdif using direct angular integration over rdir,tdir,
            // since Delta-Eddington rdif formula is not well-behaved (it is usually
            // biased low and can even be negative); use ngmax angles and gaussian
            // integration for most accuracy:
            const double R1 = rdif_a[i]; // use R1 as temporary
            const double T1 = tdif_a[i]; // use T1 as temporary
            double swt = c0;
            double smr = c0;
            double smt = c0;

            for (int ng = 0; ng < ngmax; ++ng) {
              const double mu  = difgauspt[ng];
              const double gwt = difgauswt[ng];
              swt = swt + mu*gwt;
              const double trn = std::max(exp_min, exp(-ts/mu));
              // alp = alpha(ws,mu,gs,lm)
              // gam = agamm(ws,mu,gs,lm)
              alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu));
              gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu));
              apg = alp + gam;
              amg = alp - gam;
              const double rdr = apg*R1 + amg*T1*trn - amg;
              const double tdr = apg*T1 + amg*R1*trn - apg*trn + trn;
              smr = smr + mu*rdr*gwt;
              smt = smt + mu*tdr*gwt;
            }

            rdif_a[i] = smr/swt;
            tdif_a[i] = smt/swt;

            // homogeneous layer
            rdif_b[i] = rdif_a[i];
            tdif_b[i] = tdif_a[i];     
          } // end if trntdr(i) > trmin

          // Calculate the solar beam transmission, total transmission, and
          // reflectivity for diffuse radiation from below at interface i,
          // the top of the current layer i:
          //
          //              layers       interface
          //
          //       ---------------------  i-1
          //                i-1
          //       ---------------------  i
          //                 i
          //       ---------------------

          trndir[i+1] = trndir[i] * trnlay[i];
          const double refkm1 = c1/(c1 - rdndif[i] * rdif_a[i]);
          const double tdrrdir = trndir[i] * rdir[i];
          const double tdndif = trntdr[i] - trndir[i];
          trntdr[i+1] = trndir[i] * tdir[i] + (tdndif + tdrrdir * rdndif[i]) * refkm1 * tdif_a[i];
          rdndif[i+1] = rdif_b[i] + (tdif_b[i] * rdndif[i] * refkm1 * tdif_a[i]);
          trndif[i+1] = trndif[i] * refkm1 * tdif_a[i];
        } // end main level loop

        // compute reflectivity to direct and diffuse radiation for layers
        // below by adding succesive layers starting from the underlying
        // ground and working upwards:
        //
        //              layers       interface
        //
        //       ---------------------  i
        //                 i
        //       ---------------------  i+1
        //                i+1
        //       ---------------------

        // set the underlying ground albedo == albedo of near-IR
        // unless bnd_idx == 1, for visible
        rupdir[snl_btm_itf] = albsoi[1];
        rupdif[snl_btm_itf] = albsoi[1];
        if (bnd_idx == 0) {
          rupdir[snl_btm_itf] = albsoi[0];
          rupdif[snl_btm_itf] = albsoi[0];
        }

        for (int i = snl_btm; i >= snl_top; --i) {
          // interface scattering
          const double refkp1 = c1 / ( c1 - rdif_b[i] * rupdif[i+1]);
          // dir from top layer plus exp tran ref from lower layer, interface
          // scattered and tran thru top layer from below, plus diff tran ref
          // from lower layer with interface scattering tran thru top from below
          rupdir[i] = rdir[i] + (trnlay[i]  * rupdir[i+1] + (tdir[i] - trnlay[i]) * rupdif[i+1]) * refkp1 * tdif_b[i];
          // dif from top layer from above, plus dif tran upwards reflected and
          // interface scattered which tran top from below
          rupdif[i] = rdif_a[i] + tdif_a[i]*rupdif[i+1]*refkp1*tdif_b[i];
        }

        // net flux (down-up) at each layer interface from the
        // snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
        // the interface reflectivities and transmissivities required
        // to evaluate interface fluxes are returned from solution_dEdd;
        // now compute up and down fluxes for each interface, using the
        // combined layer properties at each interface:
        //
        //              layers       interface
        //
        //       ---------------------  i
        //                 i
        //       ---------------------

        double refk, dfdir, dfdif;
        for (int i = snl_top; i <= snl_btm_itf; ++i) {
          // interface scattering
          refk = c1 / (c1 - rdndif[i] * rupdif[i]);
          // dir tran ref from below times interface scattering, plus diff
          // tran and ref from below times interface scattering
          // fdirup(i) = (trndir(i)*rupdir(i) + &
          //                 (trntdr(i)-trndir(i))  &
          //                 *rupdif(i))*refk
          // dir tran plus total diff trans times interface scattering plus
          // dir tran with up dir ref and down dif ref times interface scattering
          // fdirdn(i) = trndir(i) + (trntdr(i) &
          //               - trndir(i) + trndir(i)  &
          //               *rupdir(i)*rdndif(i))*refk
          // diffuse tran ref from below times interface scattering
          // fdifup(i) = trndif(i)*rupdif(i)*refk
          // diffuse tran times interface scattering
          // fdifdn(i) = trndif(i)*refk

          // netflux, down - up
          // dfdir = fdirdn - fdirup
          dfdir[i] = trndir[i] + (trntdr[i]-trndir[i]) * (c1 - rupdif[i]) * refk -  trndir[i]*rupdir[i]  * (c1 - rdndif[i]) * refk;
          if (dfdir[i] < puny) dfdir[i] = c0;
          // dfdif = fdifdn - fdifup
          dfdif[i] = trndif[i] * (c1 - rupdif[i]) * refk;
          if (dfdif[i] < puny) dfdif[i] = c0;
        }

        // SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
        // direct incident
        double albedo, dftmp, F_sfc_pls;
        if (flg_slr_in == 1) {
          albedo = rupdir[snl_top];
          dftmp  = dfdir;
          refk   = c1 / (c1 - rdndif[snl_top] * rupdif[snl_top]);
          F_sfc_pls = (trndir[snl_top] * rupdir[snl_top] + (trntdr[snl_top] - trndir[snl_top]) * rupdif(snl_top))*refk;
        // diffuse incident
        } else {
          albedo = rupdif[snl_top];
          dftmp  = dfdif;
          refk   = c1 / (c1 - rdndif[snl_top] * rupdif[snl_top]);
          F_sfc_pls = trndif[snl_top] * rupdif[snl_top] * refk;
        }

        // Absorbed flux in each layer
        for (int i = snl_top; i <= snl_btm; ++i) {
          F_abs[i] = dftmp[i] - dftmp[i+1];
          flx_abs_lcl[i][bnd_idx] = F_abs[i];

          // ERROR check: negative absorption
          if (flx_abs_lcl[i][bnd_idx] < -0.00001) {
            throw std::runtime_error("ELM ERROR: SNICAR negative absoption.");
          }
        }

        // absobed flux by the underlying ground
        const double F_btm_net = dftmp[snl_btm_itf];

        // note here, snl_btm_itf = 1 by snow column set up in ELM
        flx_abs_lcl[nlevsno][bnd_idx] = F_btm_net;

        if (flg_nosnl == 1) {
          // If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
          // flx_abs_lcl(:,bnd_idx) = 0._r8
          // flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

          // changed on 20070408:
          // OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
          // handles the case of no snow layers. Then, if a snow layer is addded between now and
          // SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
          flx_abs_lcl[nlevsno-1][bnd_idx] = F_abs[nlevsno-1];
          flx_abs_lcl[nlevsno][bnd_idx] = F_btm_net;
        }

        //Underflow check (we've already tripped the error condition above)
        do i=snl_top,1,1
        for (int i = snl_top; i <= nlevsno; ++i) {
          if (flx_abs_lcl[i][bnd_idx] < 0.0) {
            flx_abs_lcl[i][bnd_idx] = 0.0;
          }
        }

        double F_abs_sum = 0.0;
        for (int i = snl_top; i <= snl_btm; ++i) {
           F_abs_sum = F_abs_sum + F_abs[i];
        }


        // Energy conservation check:
        // Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
        const double energy_sum = (mu_not*pi*flx_slrd_lcl[bnd_idx]) + flx_slri_lcl[bnd_idx] - (F_abs_sum + F_btm_net + F_sfc_pls);
        if (std::abs(energy_sum) > 0.00001) {
          throw std::runtime_error("ELM ERROR: SNICAR Energy conservation error.");
        }

        albout_lcl[bnd_idx] = albedo;
        // Check that albedo is less than 1
        if (albout_lcl[bnd_idx] > 1.0) {
          throw std::runtime_error("ELM ERROR: SNICAR Albedo > 1.0.");
        }
      } // end bnd_idx waveband loop
    } // end if coszen > 0 && h2osno > min_snow 
  } // end !urbpoi
}





/*

inputs:
urbpoi                             [bool]    true if urban point, false otherwise
flg_slr_in                         [int]     flag: ==1 for direct-beam incident flux, ==2 for diffuse incident flux
snl_top                            [int]     top snow layer index [idx]
coszen                             [double]  solar zenith angle factor 
h2osno            [double]   snow water (mm H2O)
albsoi[numrad]                     [double]  albedo of surface underlying snow [frc]
albout_lcl[numrad_snw]             [double]  snow albedo by band [frc]
flx_abs_lcl[nlevsno+1][numrad_snw] [double]  absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]

outputs:
albout[numrad]                [double] snow albedo, averaged into 2 bands (=0 if no sun or no snow) [frc]
flx_abs[nlevsno+1][numrad]    [double] absorbed flux in each layer per unit flux incident [frc]

*/
void SnowAlbedoRadiationFlux(const bool &urbpoi, const int &flg_slr_in, const int &snl_top, const double &coszen,
const double &h2osno, const ArrayD1 albsoi, const ArrayD1 albout_lcl, const ArrayD2 flx_abs_lcl, 
ArrayD1 albout, ArrayD1 flx_abs) {

  if (!Land.urbpoi) {
    if ((coszen > 0.0) && (h2osno > min_snw)) {

      // Incident flux weighting parameters
      //  - sum of all VIS bands must equal 1
      //  - sum of all NIR bands must equal 1
      //
      // Spectral bands (5-band case)
      //  Band 1: 0.3-0.7um (VIS)
      //  Band 2: 0.7-1.0um (NIR)
      //  Band 3: 1.0-1.2um (NIR)
      //  Band 4: 1.2-1.5um (NIR)
      //  Band 5: 1.5-5.0um (NIR)

      // 5-band weights
      // Direct:
      double flx_wgt[numrad_snw];
      if (flg_slr_in == 1) {
         flx_wgt[0] = 1.0;
         flx_wgt[1] = 0.49352158521175;
         flx_wgt[2] = 0.18099494230665;
         flx_wgt[3] = 0.12094898498813;
         flx_wgt[4] = 0.20453448749347;
      // Diffuse:
      } else if (flg_slr_in == 2) {
         flx_wgt[0] = 1.0;
         flx_wgt[1] = 0.58581507618433;
         flx_wgt[2] = 0.20156903770812;
         flx_wgt[3] = 0.10917889346386;
         flx_wgt[4] = 0.10343699264369;
      }

      // Weight output NIR albedo appropriately
      albout[0] = albout_lcl[0];
      double flx_sum = 0.0;
      double flx_wgt_sum = 0.0;
      for (int bnd_idx = nir_bnd_bgn; bnd_idx <= nir_bnd_end; ++bnd_idx) {
        flx_sum += flx_wgt[bnd_idx] * albout_lcl[bnd_idx];
        flx_wgt_sum += flx_wgt[bnd_idx];
      }
      albout[1] = flx_sum / flx_wgt_sum;
  
      // Weight output NIR absorbed layer fluxes (flx_abs) appropriately
      for (int i = 0; i <= nlevsno; ++i) { flx_abs[i][0] = flx_abs_lcl[i][0]; }
  
      for (int i = snl_top; i <= nlevsno; ++i) {
        flx_sum = 0.0;
        for (int bnd_idx = nir_bnd_bgn; bnd_idx <= nir_bnd_end; ++bnd_idx) {
          flx_sum += flx_wgt[bnd_idx] * flx_abs_lcl[i][bnd_idx];
        }
        flx_abs[i][1] = flx_sum / flx_wgt_sum;
      }

      // near-IR direct albedo/absorption adjustment for high solar zenith angles
      // solar zenith angle parameterization
      // calculate the scaling factor for NIR direct albedo if SZA>75 degree
      // coefficients used for SZA parameterization
      if ((mu_not < mu_75) && (flg_slr_in == 1)) {
        const double sza_c1 = sza_a0 + sza_a1 * mu_not + sza_a2 * pow(mu_not, 2.0); // coefficient, SZA parameteirzation
        const double sza_c0 = sza_b0 + sza_b1 * mu_not + sza_b2 * pow(mu_not, 2.0); // coefficient, SZA parameterization
        const double sza_factor = sza_c1 * (log10(snw_rds_lcl[snl_top] * c1) - c6) + sza_c0; // factor used to adjust NIR direct albedo
        const double flx_sza_adjust  = albout[1] * (sza_factor - c1) * flx_wgt_sum; // direct NIR flux adjustment from sza_factor
        albout[1] *= sza_factor;
        flx_abs[snl_top][1] -= flx_sza_adjust;
      }
    // If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
    } else if ((coszen > 0.0) && (h2osno < min_snw) && (h2osno > 0.0)) {
      albout[0] = albsoi[0];
      albout[1] = albsoi[1]
    // There is either zero snow, or no sun
    } else {
      albout[0] = 0.0;
      albout[1] = 0.0;
    } // if ((coszen > 0.0) && (h2osno > min_snw))
  } // if !urbpoi
} // SnowAlbedoRadiationFlux()


} // namespace SNICAR
} // namespace ELM
