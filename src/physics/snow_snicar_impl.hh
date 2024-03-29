// functions derived from subroutine SNICAR_AD_RT() in SnowSNICARMod.F90

#pragma once

namespace ELM::snow_snicar {

template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
ACCELERATE
void init_timestep(const bool& urbpoi, const int& flg_slr_in, const double& coszen, const double& h2osno, const int& snl,
                   const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, const ArrayD1 snw_rds, int& snl_top,
                   int& snl_btm, ArrayD2 flx_abs_lcl, ArrayD2 flx_abs, int& flg_nosnl, ArrayD1 h2osoi_ice_lcl,
                   ArrayD1 h2osoi_liq_lcl, ArrayI1 snw_rds_lcl, double& mu_not, ArrayD1 flx_slrd_lcl,
                   ArrayD1 flx_slri_lcl)
{
  using ELMdims::nlevsno;
  using ELMdims::numrad;
  using ELMdims::numrad_snw;
  using ELMconst::SNW_RDS_MIN;

  if (!urbpoi) {
    // Zero absorbed radiative fluxes:
    for (int i = 0; i <= nlevsno(); ++i) {
      for (int ib = 0; i < numrad(); ++i) {
        flx_abs(i, ib) = 0.0;
      }
    }
    for (int i = 0; i <= nlevsno(); ++i) {
      for (int ib = 0; i < numrad_snw(); ++i) {
        flx_abs_lcl(i, ib) = 0.0;
      }
    }

    // Qualifier for computing snow RT:
    //  1) sunlight from atmosphere model
    //  2) minimum amount of snow on ground.
    //     Otherwise, set snow albedo to zero
    if ((coszen > 0.0) && (h2osno > detail::min_snw)) {

      // If there is snow, but zero snow layers, we must create a layer locally.
      // This layer is presumed to have the fresh snow effective radius.
      int snl_lcl;
      if (snl == 0) {
        flg_nosnl = 1;
        snl_lcl = 1;
        h2osoi_ice_lcl(nlevsno() - 1) = h2osno;
        h2osoi_liq_lcl(nlevsno() - 1) = 0.0;
        snw_rds_lcl(nlevsno() - 1) = round(SNW_RDS_MIN());
      } else {
        flg_nosnl = 0;
        snl_lcl = snl;
        for (int i = 0; i < nlevsno(); ++i) {
          h2osoi_liq_lcl(i) = h2osoi_liq(i);
          h2osoi_ice_lcl(i) = h2osoi_ice(i);
          snw_rds_lcl(i) = round(snw_rds(i));
        }
      }

      snl_btm = nlevsno() - 1;       // index of bottom snow layer
      snl_top = nlevsno() - snl_lcl; // index of top snow layer

      // Assume fixed BC effective radii of 100nm. This is close to
      // the effective radius of 95nm (number median radius of
      // 40nm) assumed for freshly-emitted BC in MAM.  Future
      // implementations may prognose the BC effective radius in
      // snow.

      // put this in header
      // for (int i = 0; i < nlevsno(); ++i) {
      //  rds_bcint_lcl[i] = 100.0;
      //  rds_bcext_lcl[i] = 100.0;
      //}

      // Error check for snow grain size:
      for (int i = snl_top; i <= snl_btm; ++i) {
        if ((snw_rds_lcl(i) < detail::snw_rds_min_tbl) || (snw_rds_lcl(i) > detail::snw_rds_max_tbl)) {
          throw std::runtime_error("ELM ERROR: SNICAR snow grain radius out of bounds.");
        }
      }

      // mu_not is cosine solar zenith angle above the fresnel level; make
      // sure mu_not is large enough for stable and meaningful radiation
      // solution: .01 is like sun just touching horizon with its lower edge
      // equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
      mu_not = std::max(coszen, alg::cp01);

      // Set direct or diffuse incident irradiance to 1
      // (This has to be within the bnd loop because mu_not is adjusted in rare cases)
      if (flg_slr_in == 1) {
        for (int bnd_idx = 0; bnd_idx < numrad_snw(); ++bnd_idx) {
          flx_slrd_lcl(bnd_idx) = 1.0 / (mu_not * ELMconst::ELM_PI()); // this corresponds to incident irradiance of 1.0
          flx_slri_lcl(bnd_idx) = 0.0;
        }
      } else if (flg_slr_in == 2) {
        for (int bnd_idx = 0; bnd_idx < numrad_snw(); ++bnd_idx) {
          flx_slrd_lcl(bnd_idx) = 0.0;
          flx_slri_lcl(bnd_idx) = 1.0;
        }
      } else {
        throw std::runtime_error("ELM ERROR: SNICAR solar waveband flag out of bounds - only 2 bands supported.");
      }
    }
  }
}

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
                             SubviewD2 tau_star)
{
  using ELMdims::sno_nbr_aer;
  using ELMdims::nlevsno;
  using ELMdims::numrad_snw;
  static constexpr double rds_bcint_lcl{100.0}; // effective radius of within-ice BC [nm]
  static constexpr double rds_bcext_lcl{100.0}; // effective radius of external BC [nm]
  static constexpr int DELTA{1};                // always use Delta approximation for snow

  // Define local Mie parameters based on snow grain size and aerosol species,
  //  retrieved from a lookup table.
  if (!urbpoi) {
    if ((coszen > 0.0) && (h2osno > detail::min_snw)) {

      // Set local aerosol array
      double mss_cnc_aer_lcl[nlevsno()][sno_nbr_aer()];
      for (int i = 0; i < nlevsno(); ++i) {
        for (int j = 0; j < sno_nbr_aer(); ++j) {
          mss_cnc_aer_lcl[i][j] = mss_cnc_aer_in(i, j);
        }
      }

      for (int bnd_idx = 0; bnd_idx < numrad_snw(); ++bnd_idx) {

        if ((numrad_snw() == 5) && ((bnd_idx == 4) || (bnd_idx == 3))) {
          for (int i = 0; i < nlevsno(); ++i) {
            for (int j = 0; j < sno_nbr_aer(); ++j) {
              mss_cnc_aer_lcl[i][j] = 0.0;
            }
          }
        }

        double ss_alb_snw_lcl[nlevsno()];
        double asm_prm_snw_lcl[nlevsno()];
        double ext_cff_mss_snw_lcl[nlevsno()];
        if (flg_slr_in == 1) {
          for (int i = snl_top; i <= snl_btm; ++i) {
            const int rds_idx = snw_rds_lcl(i) - detail::snw_rds_min_tbl;
            // snow optical properties (direct radiation)
            ss_alb_snw_lcl[i] = ss_alb_snw_drc(bnd_idx, rds_idx);
            asm_prm_snw_lcl[i] = asm_prm_snw_drc(bnd_idx, rds_idx);
            ext_cff_mss_snw_lcl[i] = ext_cff_mss_snw_drc(bnd_idx, rds_idx);
          }
        } else if (flg_slr_in == 2) {
          for (int i = snl_top; i <= snl_btm; ++i) {
            const int rds_idx = snw_rds_lcl(i) - detail::snw_rds_min_tbl;
            // snow optical properties (diffuse radiation)
            ss_alb_snw_lcl[i] = ss_alb_snw_dfs(bnd_idx, rds_idx);
            asm_prm_snw_lcl[i] = asm_prm_snw_dfs(bnd_idx, rds_idx);
            ext_cff_mss_snw_lcl[i] = ext_cff_mss_snw_dfs(bnd_idx, rds_idx);
          }
        }

        // H. Wang
        //  aerosol species 3 optical properties
        double ss_alb_aer_lcl[sno_nbr_aer()];
        double asm_prm_aer_lcl[sno_nbr_aer()];
        double ext_cff_mss_aer_lcl[sno_nbr_aer()];
        ss_alb_aer_lcl[2] = ss_alb_oc1(bnd_idx);
        asm_prm_aer_lcl[2] = asm_prm_oc1(bnd_idx);
        ext_cff_mss_aer_lcl[2] = ext_cff_mss_oc1(bnd_idx);

        //  aerosol species 4 optical properties
        ss_alb_aer_lcl[3] = ss_alb_oc2(bnd_idx);
        asm_prm_aer_lcl[3] = asm_prm_oc2(bnd_idx);
        ext_cff_mss_aer_lcl[3] = ext_cff_mss_oc2(bnd_idx);

        // aerosol species 5 optical properties
        ss_alb_aer_lcl[4] = ss_alb_dst1(bnd_idx);
        asm_prm_aer_lcl[4] = asm_prm_dst1(bnd_idx);
        ext_cff_mss_aer_lcl[4] = ext_cff_mss_dst1(bnd_idx);

        // aerosol species 6 optical properties
        ss_alb_aer_lcl[5] = ss_alb_dst2(bnd_idx);
        asm_prm_aer_lcl[5] = asm_prm_dst2(bnd_idx);
        ext_cff_mss_aer_lcl[5] = ext_cff_mss_dst2(bnd_idx);

        // aerosol species 7 optical properties
        ss_alb_aer_lcl[6] = ss_alb_dst3(bnd_idx);
        asm_prm_aer_lcl[6] = asm_prm_dst3(bnd_idx);
        ext_cff_mss_aer_lcl[6] = ext_cff_mss_dst3(bnd_idx);

        // aerosol species 8 optical properties
        ss_alb_aer_lcl[7] = ss_alb_dst4(bnd_idx);
        asm_prm_aer_lcl[7] = asm_prm_dst4(bnd_idx);
        ext_cff_mss_aer_lcl[7] = ext_cff_mss_dst4(bnd_idx);

        // 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
        // 2. optical Depths (tau_snw, tau_aer)
        // 3. weighted Mie properties (tau, omega, g)

        // Weighted Mie parameters of each layer
        double tau[nlevsno()];
        double omega[nlevsno()];
        double g[nlevsno()];
        for (int i = snl_top; i <= snl_btm; ++i) {
          // mgf++ within-ice and external BC optical properties
          // Lookup table indices for BC optical properties,
          // dependent on snow grain size and BC particle
          // size.

          // valid for 25 < snw_rds < 1625 um:
          int idx_bcint_icerds;
          if (snw_rds_lcl(i) < 125) {
            double tmp1 = snw_rds_lcl(i) / 50;
            idx_bcint_icerds = round(tmp1) - 1;
          } else if (snw_rds_lcl(i) < 175) {
            idx_bcint_icerds = 1;
          } else {
            double tmp1 = (snw_rds_lcl(i) / 250) + 2;
            idx_bcint_icerds = round(tmp1) - 1;
          }

          // valid for 25 < bc_rds < 525 nm
          int idx_bcint_nclrds = round(rds_bcint_lcl / 50) - 1;
          int idx_bcext_nclrds = round(rds_bcext_lcl / 50) - 1;

          // check bounds:
          if (idx_bcint_icerds < detail::idx_bcint_icerds_min)
            idx_bcint_icerds = detail::idx_bcint_icerds_min;
          if (idx_bcint_icerds > detail::idx_bcint_icerds_max)
            idx_bcint_icerds = detail::idx_bcint_icerds_max;
          if (idx_bcint_nclrds < detail::idx_bc_nclrds_min)
            idx_bcint_nclrds = detail::idx_bc_nclrds_min;
          if (idx_bcint_nclrds > detail::idx_bc_nclrds_max)
            idx_bcint_nclrds = detail::idx_bc_nclrds_max;
          if (idx_bcext_nclrds < detail::idx_bc_nclrds_min)
            idx_bcext_nclrds = detail::idx_bc_nclrds_min;
          if (idx_bcext_nclrds > detail::idx_bc_nclrds_max)
            idx_bcext_nclrds = detail::idx_bc_nclrds_max;

          // retrieve absorption enhancement factor for within-ice BC
          double enh_fct = bcenh(idx_bcint_icerds, idx_bcint_nclrds, bnd_idx);

          // get BC optical properties (moved from above)
          // aerosol species 1 optical properties (within-ice BC)
          ss_alb_aer_lcl[0] = ss_alb_bc1(idx_bcint_nclrds, bnd_idx);
          asm_prm_aer_lcl[0] = asm_prm_bc1(idx_bcint_nclrds, bnd_idx);
          ext_cff_mss_aer_lcl[0] = ext_cff_mss_bc1(idx_bcint_nclrds, bnd_idx) * enh_fct;

          // aerosol species 2 optical properties (external BC)
          ss_alb_aer_lcl[1] = ss_alb_bc2(idx_bcext_nclrds, bnd_idx);
          asm_prm_aer_lcl[1] = asm_prm_bc2(idx_bcext_nclrds, bnd_idx);
          ext_cff_mss_aer_lcl[1] = ext_cff_mss_bc2(idx_bcext_nclrds, bnd_idx);

          double L_snw = h2osoi_ice_lcl(i) + h2osoi_liq_lcl(i);
          double tau_snw = L_snw * ext_cff_mss_snw_lcl[i];

          double tau_aer[sno_nbr_aer()];
          for (int j = 0; j < sno_nbr_aer(); ++j) {
            double L_aer = L_snw * mss_cnc_aer_lcl[i][j];
            tau_aer[j] = L_aer * ext_cff_mss_aer_lcl[j];
          }

          double tau_sum = 0.0;
          double omega_sum = 0.0;
          double g_sum = 0.0;

          for (int j = 0; j < sno_nbr_aer(); ++j) {
            tau_sum += tau_aer[j];
            omega_sum += (tau_aer[j] * ss_alb_aer_lcl[j]);
            g_sum += (tau_aer[j] * ss_alb_aer_lcl[j] * asm_prm_aer_lcl[j]);
          }

          tau[i] = tau_sum + tau_snw;
          omega[i] = (1.0 / tau[i]) * (omega_sum + (ss_alb_snw_lcl[i] * tau_snw));
          g[i] = (1.0 / (tau[i] * omega[i])) * (g_sum + (asm_prm_snw_lcl[i] * ss_alb_snw_lcl[i] * tau_snw));

        } // end Weighted Mie snl loop

        // DELTA transformations, if requested
        if (DELTA == 1) {
          for (int i = snl_top; i <= snl_btm; ++i) {
            g_star(bnd_idx, i) = g[i] / (1.0 + g[i]);
            omega_star(bnd_idx, i) = ((1.0 - pow(g[i], 2.0)) * omega[i]) / (1.0 - (omega[i] * pow(g[i], 2.0)));
            tau_star(bnd_idx, i) = (1.0 - (omega[i] * pow(g[i], 2.0))) * tau[i];
          }
        } else {
          for (int i = snl_top; i <= snl_btm; ++i) {
            g_star(bnd_idx, i) = g[i];
            omega_star(bnd_idx, i) = omega[i];
            tau_star(bnd_idx, i) = tau[i];
          }
        }
      } // end bnd_idx waveband loop
    }   // end if coszen > 0 && h2osno > min_snow
  }     // end !urbpoi
}

template <typename ArrayD1, typename ArrayD2>
ACCELERATE
void snow_radiative_transfer_solver(const bool& urbpoi, const int& flg_slr_in, const int& flg_nosnl, const int& snl_top,
                                    const int& snl_btm, const double& coszen, const double& h2osno,
                                    const double& mu_not, const ArrayD1 flx_slrd_lcl, const ArrayD1 flx_slri_lcl,
                                    const ArrayD1 albsoi, const ArrayD2 g_star, const ArrayD2 omega_star,
                                    const ArrayD2 tau_star, ArrayD1 albout_lcl, ArrayD2 flx_abs_lcl)
{
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
   --    0                  4  ---  [nlevsno() - 1]
  | 0|                         |4|
  -----  1 ground interface 5 ----- [nlevsno()]
  | 1|                         |5|
   --    2                  6  ---  [nlevsno() + 1]  */

  using ELMdims::nlevsno;

  // Gaussian integration angle and coefficients for diffuse radiation
  static constexpr double difgauspt[8]      // gaussian angles (radians)
      = {0.9894009, 0.9445750, 0.8656312, 0.7554044, 0.6178762, 0.4580168, 0.2816036, 0.0950125};
  static constexpr double difgauswt[8]      // gaussian weights
      = {0.0271525, 0.0622535, 0.0951585, 0.1246290, 0.1495960, 0.1691565, 0.1826034, 0.1894506};

  static constexpr int ngmax{8};                 // gaussian integration index
  static constexpr double puny{1.0e-11};         // very small number used in this function
  static constexpr double argmax{10.0};          // maximum argument of exponential
  static constexpr double exp_min{exp(-argmax)}; // minimum exponential value

  if (!urbpoi) {
    if ((coszen > 0.0) && (h2osno > detail::min_snw)) {
      // local interface reflect/transmit vars
      double trndir[nlevsno() + 1]; // solar beam down transmission from top
      double trntdr[nlevsno() + 1]; // total transmission to direct beam for layers above
      double trndif[nlevsno() + 1]; // diffuse transmission to diffuse beam for layers above
      double rupdir[nlevsno() + 1]; // reflectivity to direct radiation for layers below
      double rupdif[nlevsno() + 1]; // reflectivity to diffuse radiation for layers below
      double rdndif[nlevsno() + 1]; // reflectivity to diffuse radiation for layers above
      double dfdir[nlevsno() + 1];  // down-up flux at interface due to direct beam at top surface
      double dfdif[nlevsno() + 1];  // down-up flux at interface due to diffuse beam at top surface
      double dftmp[nlevsno() + 1];  // temporary variable for down-up flux at interface
      // local layer reflect/transmit vars
      double rdir[nlevsno()];   // layer reflectivity to direct radiation
      double rdif_a[nlevsno()]; // layer reflectivity to diffuse radiation from above
      double rdif_b[nlevsno()]; // layer reflectivity to diffuse radiation from below
      double tdir[nlevsno()];   // layer transmission to direct radiation (solar beam + diffuse)
      double tdif_a[nlevsno()]; // layer transmission to diffuse radiation from above
      double tdif_b[nlevsno()]; // layer transmission to diffuse radiation from below
      double trnlay[nlevsno()]; // solar beam transm for layer (direct beam only)
      // net absorbed radiative energy (lyr) [W/m^2]
      double F_abs[nlevsno()];

      const int snl_btm_itf = nlevsno(); // index of ground/snow interface (same as snl_btm + 1)

      for (int bnd_idx = 0; bnd_idx < ELMdims::numrad_snw(); ++bnd_idx) {

        for (int i = snl_top; i <= snl_btm_itf; ++i) {
          trndir[i] = alg::c0;
          trntdr[i] = alg::c0;
          trndif[i] = alg::c0;
          rupdir[i] = alg::c0;
          rupdif[i] = alg::c0;
          rdndif[i] = alg::c0;
        }

        // initialize top interface of top layer
        trndir[snl_top] = alg::c1;
        trntdr[snl_top] = alg::c1;
        trndif[snl_top] = alg::c1;
        rdndif[snl_top] = alg::c0;

        // begin main level loop
        // for layer interfaces except for the very bottom
        for (int i = snl_top; i <= snl_btm; ++i) {
          // initialize all layer apparent optical properties to 0
          rdir[i] = alg::c0;
          rdif_a[i] = alg::c0;
          rdif_b[i] = alg::c0;
          tdir[i] = alg::c0;
          tdif_a[i] = alg::c0;
          tdif_b[i] = alg::c0;
          trnlay[i] = alg::c0;

          // compute next layer Delta-eddington solution only if total transmission
          // of radiation to the interface just above the layer exceeds trmin.

          if (trntdr[i] > alg::trmin) {

            // calculation over layers with penetrating radiation

            // delta-transformed single-scattering properties
            // of this layer
            const double ts{tau_star(bnd_idx, i)};
            const double ws{omega_star(bnd_idx, i)};
            const double gs{g_star(bnd_idx, i)};

            // Delta-Eddington solution expressions
            // n(uu,et)         = ((uu+alg::c1)*(uu+alg::c1)/et ) - ((uu-alg::c1)*(uu-alg::c1)*et)
            // u(w,gg,e)        = alg::c1p5*(alg::c1 - w*gg)/e
            // el(w,gg)         = sqrt(alg::c3*(alg::c1-w)*(alg::c1 - w*gg))
            const double lm{std::sqrt(alg::c3 * (alg::c1 - ws) * (alg::c1 - ws * gs))}; // lm = el(ws,gs)
            const double ue{alg::c1p5 * (alg::c1 - ws * gs) / lm};                 // ue = u(ws,gs,lm)
            const double extins{std::max(exp_min, exp(-lm * ts))};
            const double ne{((ue + alg::c1) * (ue + alg::c1) / extins) - ((ue - alg::c1) * (ue - alg::c1) * extins)}; // ne = n(ue,extins)

            // first calculation of rdif, tdif using Delta-Eddington formulas
            // rdif_a(k) = (ue+alg::c1)*(ue-alg::c1)*(alg::c1/extins - extins)/ne
            rdif_a[i] = (pow(ue, 2.0) - alg::c1) * (alg::c1 / extins - extins) / ne;
            tdif_a[i] = alg::c4 * ue / ne;

            // evaluate rdir,tdir for direct beam
            trnlay[i] = std::max(exp_min, exp(-ts / mu_not));

            // Delta-Eddington solution expressions
            // alpha(w,uu,gg,e) = p75*w*uu*((alg::c1 + gg*(alg::c1-w))/(alg::c1 - e*e*uu*uu))
            // agamm(w,uu,gg,e) = p5*w*((alg::c1 + alg::c3*gg*(alg::c1-w)*uu*uu)/(alg::c1-e*e*uu*uu))
            // alp = alpha(ws,mu_not,gs,lm)
            // gam = agamm(ws,mu_not,gs,lm)
            double alp = alg::cp75 * ws * mu_not * ((alg::c1 + gs * (alg::c1 - ws)) / (alg::c1 - lm * lm * mu_not * mu_not));
            double gam = alg::cp5 * ws * ((alg::c1 + alg::c3 * gs * (alg::c1 - ws) * mu_not * mu_not) / (alg::c1 - lm * lm * mu_not * mu_not));
            double apg = alp + gam;
            double amg = alp - gam;

            rdir[i] = apg * rdif_a[i] + amg * (tdif_a[i] * trnlay[i] - alg::c1);
            tdir[i] = apg * tdif_a[i] + (amg * rdif_a[i] - apg + alg::c1) * trnlay[i];

            // recalculate rdif,tdif using direct angular integration over rdir,tdir,
            // since Delta-Eddington rdif formula is not well-behaved (it is usually
            // biased low and can even be negative); use ngmax angles and gaussian
            // integration for most accuracy:
            const double R1{rdif_a[i]}; // use R1 as temporary
            const double T1{tdif_a[i]}; // use T1 as temporary
            double swt = alg::c0;
            double smr = alg::c0;
            double smt = alg::c0;

            for (int ng = 0; ng < ngmax; ++ng) {
              const double mu{difgauspt[ng]};
              const double gwt{difgauswt[ng]};
              swt = swt + mu * gwt;
              const double trn{std::max(exp_min, exp(-ts / mu))};
              // alp = alpha(ws,mu,gs,lm)
              // gam = agamm(ws,mu,gs,lm)
              alp = alg::cp75 * ws * mu * ((alg::c1 + gs * (alg::c1 - ws)) / (alg::c1 - lm * lm * mu * mu));
              gam = alg::cp5 * ws * ((alg::c1 + alg::c3 * gs * (alg::c1 - ws) * mu * mu) / (alg::c1 - lm * lm * mu * mu));
              apg = alp + gam;
              amg = alp - gam;
              const double rdr = apg * R1 + amg * T1 * trn - amg;
              const double tdr = apg * T1 + amg * R1 * trn - apg * trn + trn;
              smr = smr + mu * rdr * gwt;
              smt = smt + mu * tdr * gwt;
            }

            rdif_a[i] = smr / swt;
            tdif_a[i] = smt / swt;

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

          trndir[i + 1] = trndir[i] * trnlay[i];
          const double refkm1{alg::c1 / (alg::c1 - rdndif[i] * rdif_a[i])};
          const double tdrrdir{trndir[i] * rdir[i]};
          const double tdndif{trntdr[i] - trndir[i]};
          trntdr[i + 1] = trndir[i] * tdir[i] + (tdndif + tdrrdir * rdndif[i]) * refkm1 * tdif_a[i];
          rdndif[i + 1] = rdif_b[i] + (tdif_b[i] * rdndif[i] * refkm1 * tdif_a[i]);
          trndif[i + 1] = trndif[i] * refkm1 * tdif_a[i];
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
        rupdir[snl_btm_itf] = albsoi(1);
        rupdif[snl_btm_itf] = albsoi(1);
        if (bnd_idx == 0) {
          rupdir[snl_btm_itf] = albsoi(0);
          rupdif[snl_btm_itf] = albsoi(0);
        }

        for (int i = snl_btm; i >= snl_top; --i) {
          // interface scattering
          const double refkp1{alg::c1 / (alg::c1 - rdif_b[i] * rupdif[i + 1])};
          // dir from top layer plus exp tran ref from lower layer, interface
          // scattered and tran thru top layer from below, plus diff tran ref
          // from lower layer with interface scattering tran thru top from below
          rupdir[i] =
              rdir[i] + (trnlay[i] * rupdir[i + 1] + (tdir[i] - trnlay[i]) * rupdif[i + 1]) * refkp1 * tdif_b[i];
          // dif from top layer from above, plus dif tran upwards reflected and
          // interface scattered which tran top from below
          rupdif[i] = rdif_a[i] + tdif_a[i] * rupdif[i + 1] * refkp1 * tdif_b[i];
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

        double refk;
        for (int i = snl_top; i <= snl_btm_itf; ++i) {
          // interface scattering
          refk = alg::c1 / (alg::c1 - rdndif[i] * rupdif[i]);
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
          dfdir[i] = trndir[i] + (trntdr[i] - trndir[i]) * (alg::c1 - rupdif[i]) * refk -
                     trndir[i] * rupdir[i] * (alg::c1 - rdndif[i]) * refk;
          if (dfdir[i] < puny)
            dfdir[i] = alg::c0;
          // dfdif = fdifdn - fdifup
          dfdif[i] = trndif[i] * (alg::c1 - rupdif[i]) * refk;
          if (dfdif[i] < puny)
            dfdif[i] = alg::c0;
        }

        // SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
        // direct incident
        double albedo, F_sfc_pls;
        if (flg_slr_in == 1) {
          albedo = rupdir[snl_top];
          for (int i = snl_top; i <= snl_btm_itf; ++i) {
            dftmp[i] = dfdir[i];
          }
          refk = alg::c1 / (alg::c1 - rdndif[snl_top] * rupdif[snl_top]);
          F_sfc_pls =
              (trndir[snl_top] * rupdir[snl_top] + (trntdr[snl_top] - trndir[snl_top]) * rupdif[snl_top]) * refk;
          // diffuse incident
        } else {
          albedo = rupdif[snl_top];
          for (int i = snl_top; i <= snl_btm_itf; ++i) {
            dftmp[i] = dfdif[i];
          }
          refk = alg::c1 / (alg::c1 - rdndif[snl_top] * rupdif[snl_top]);
          F_sfc_pls = trndif[snl_top] * rupdif[snl_top] * refk;
        }

        // Absorbed flux in each layer
        for (int i = snl_top; i <= snl_btm; ++i) {
          F_abs[i] = dftmp[i] - dftmp[i + 1];
          flx_abs_lcl(i, bnd_idx) = F_abs[i];

          // ERROR check: negative absorption
          if (flx_abs_lcl(i, bnd_idx) < -0.00001) {
            throw std::runtime_error("ELM ERROR: SNICAR negative absoption.");
          }
        }

        // absobed flux by the underlying ground
        const double F_btm_net{dftmp[snl_btm_itf]};

        // note here, snl_btm_itf = 1 by snow column set up in ELM
        flx_abs_lcl(nlevsno(), bnd_idx) = F_btm_net;

        if (flg_nosnl == 1) {
          // If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
          // flx_abs_lcl(:,bnd_idx) = 0._r8
          // flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

          // changed on 20070408:
          // OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
          // handles the case of no snow layers. Then, if a snow layer is addded between now and
          // SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
          flx_abs_lcl(nlevsno() - 1, bnd_idx) = F_abs[nlevsno() - 1];
          flx_abs_lcl(nlevsno(), bnd_idx) = F_btm_net;
        }

        // Underflow check (we've already tripped the error condition above)
        for (int i = snl_top; i <= nlevsno(); ++i) {
          if (flx_abs_lcl(i, bnd_idx) < 0.0) {
            flx_abs_lcl(i, bnd_idx) = 0.0;
          }
        }

        double F_abs_sum = 0.0;
        for (int i = snl_top; i <= snl_btm; ++i) {
          F_abs_sum = F_abs_sum + F_abs[i];
        }

        // Energy conservation check:
        // Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
        const double energy_sum{(mu_not * ELMconst::ELM_PI() * flx_slrd_lcl(bnd_idx))
            + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)};
        if (std::abs(energy_sum) > 0.00001) {
          throw std::runtime_error("ELM ERROR: SNICAR Energy conservation error.");
        }

        albout_lcl(bnd_idx) = albedo;
        // Check that albedo is less than 1
        if (albout_lcl(bnd_idx) > 1.0) {
          throw std::runtime_error("ELM ERROR: SNICAR Albedo > 1.0.");
        }
      } // end bnd_idx waveband loop
    }   // end if coszen > 0 && h2osno > min_snow
  }     // end !urbpoi
}

template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
ACCELERATE
void snow_albedo_radiation_factor(const bool& urbpoi, const int& flg_slr_in, const int& snl_top, const double& coszen,
                                  const double& mu_not, const double& h2osno, const ArrayI1 snw_rds_lcl,
                                  const ArrayD1 albsoi, const ArrayD1 albout_lcl, const ArrayD2 flx_abs_lcl,
                                  ArrayD1 albout, ArrayD2 flx_abs)
{
  using ELMdims::nlevsno;
  
  // constant coefficients used for SZA parameterization
  static constexpr double sza_a0{0.085730};
  static constexpr double sza_a1{-0.630883};
  static constexpr double sza_a2{1.303723};
  static constexpr double sza_b0{1.467291};
  static constexpr double sza_b1{-3.338043};
  static constexpr double sza_b2{6.807489};

  static constexpr int nir_bnd_bgn{1};      // first band index in near-IR spectrum [idx]
  static constexpr int nir_bnd_end{4};      // ending near-IR band index [idx]
  static constexpr double mu_75{0.2588};    // cosine of 75 degree

  if (!urbpoi) {
    if ((coszen > 0.0) && (h2osno > detail::min_snw)) {

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
      double flx_wgt[ELMdims::numrad_snw()];
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
      albout(0) = albout_lcl(0);
      double flx_sum = 0.0;
      double flx_wgt_sum = 0.0;
      for (int bnd_idx = nir_bnd_bgn; bnd_idx <= nir_bnd_end; ++bnd_idx) {
        flx_sum += flx_wgt[bnd_idx] * albout_lcl(bnd_idx);
        flx_wgt_sum += flx_wgt[bnd_idx];
      }
      albout(1) = flx_sum / flx_wgt_sum;

      // Weight output NIR absorbed layer fluxes (flx_abs) appropriately
      for (int i = 0; i <= nlevsno(); ++i) {
        flx_abs(i, 0) = flx_abs_lcl(i, 0);
      }

      for (int i = snl_top; i <= nlevsno(); ++i) {
        flx_sum = 0.0;
        for (int bnd_idx = nir_bnd_bgn; bnd_idx <= nir_bnd_end; ++bnd_idx) {
          flx_sum += flx_wgt[bnd_idx] * flx_abs_lcl(i, bnd_idx);
        }
        flx_abs(i, 1) = flx_sum / flx_wgt_sum;
      }

      // near-IR direct albedo/absorption adjustment for high solar zenith angles
      // solar zenith angle parameterization
      // calculate the scaling factor for NIR direct albedo if SZA>75 degree
      // coefficients used for SZA parameterization
      if ((mu_not < mu_75) && (flg_slr_in == 1)) {
        const double sza_c1{sza_a0 + sza_a1 * mu_not + sza_a2 * pow(mu_not, 2.0)}; // coefficient, SZA parameteirzation
        const double sza_c0{sza_b0 + sza_b1 * mu_not + sza_b2 * pow(mu_not, 2.0)}; // coefficient, SZA parameterization
        const double sza_factor{sza_c1 *
            (log10(snw_rds_lcl(snl_top) * alg::c1) - alg::c6) + sza_c0}; // factor used to adjust NIR direct albedo
        const double flx_sza_adjust{albout(1) *
            (sza_factor - alg::c1) * flx_wgt_sum}; // direct NIR flux adjustment from sza_factor
        albout(1) *= sza_factor;
        flx_abs(snl_top, 1) -= flx_sza_adjust;
      }
      // If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
    } else if ((coszen > 0.0) && (h2osno < detail::min_snw) && (h2osno > 0.0)) {
      albout(0) = albsoi(0);
      albout(1) = albsoi(1);
      // There is either zero snow, or no sun
    } else {
      albout(0) = 0.0;
      albout(1) = 0.0;
    } // if ((coszen > 0.0) && (h2osno > detail::min_snw))
  }   // if !urbpoi
} // snow_albedo_radiation_factor()

} // namespace ELM::snow_snicar
