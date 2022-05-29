
#pragma once

#include <cmath>

#include "land_data.h"

namespace ELM::snow {


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
  template <typename ArrayD1, typename ArrayD3>
  ACCELERATE
  void snow_aging(const bool& do_capsnow,
                 const int& snl,
                 const double& frac_sno,
                 const double& dtime,
                 const double& qflx_snwcp_ice,
                 const double& qflx_snow_grnd,
                 const double& h2osno,
                 const ArrayD1 dz,
                 const ArrayD1 h2osoi_liq,
                 const ArrayD1 h2osoi_ice,
                 const ArrayD1 t_soisno,
                 const ArrayD1 qflx_snofrz_lyr,
                 const SnwRdsTable<ArrayD3>& snw_table,
                 double& snot_top,
                 double& dTdz_top,
                 double& snw_rds_top,
                 double& sno_liq_top,
                 ArrayD1 snw_rds)
  {
    using ELMdims::nlevsno;
    using ELMconst::ELM_PI;
    using ELMconst::SNW_RDS_MIN;
    using snow_snicar::detail::idx_T_min;
    using snow_snicar::detail::idx_T_max;
    using snow_snicar::detail::idx_Tgrd_min;
    using snow_snicar::detail::idx_Tgrd_max;
    using snow_snicar::detail::idx_rhos_min;
    using snow_snicar::detail::idx_rhos_max;

    static constexpr double snw_rds_refrz{1000.0}; // effective radius of re-frozen snow [microns]
    static constexpr double C1_liq_Brun89{0.0}; // constant for liquid water grain growth [m3 s-1],
                                                 // from Brun89: zeroed to accomodate dry snow aging
    static constexpr double C2_liq_Brun89{4.22e-13}; // constant for liquid water grain growth [m3 s-1],
                                                       // from Brun89: corrected for LWC in units of percent
    static constexpr bool flg_snoage_scl{false}; // flag for scaling the snow aging rate by some arbitrary factor
    static constexpr double xdrdt{1.0}; // arbitrary factor applied to snow aging rate

   if (snl > 0) {

      const int snl_btm = nlevsno - 1;
      const int snl_top = nlevsno - snl;

      // loop over snow layers
      for (int i = snl_top; i <= snl_btm; ++i) {
        //
        //**********  1. DRY SNOW AGING  ***********
        //
        double h2osno_lyr = h2osoi_liq(i) + h2osoi_ice(i);

        // temperature gradient
        double t_snotop, t_snobtm; // temperature at upper and lower layer boundaries [K]
        if (i == snl_top) {
          // top layer
          t_snotop = t_soisno(snl_top);
          t_snobtm = (t_soisno(i+1)*dz(i) + t_soisno(i)*dz(i+1)) / (dz(i)+dz(i+1));
        } else {
          t_snotop = (t_soisno(i-1)*dz(i) + t_soisno(i)*dz(i-1)) / (dz(i)+dz(i-1));
          t_snobtm = (t_soisno(i+1)*dz(i) + t_soisno(i)*dz(i+1)) / (dz(i)+dz(i+1));
        }

        double cdz = frac_sno * dz(i);
        double dTdz = std::abs((t_snotop - t_snobtm) / cdz);

        // snow density
        double rhos = (h2osoi_liq(i) + h2osoi_ice(i)) / cdz;

        // make sure rhos doesn't drop below 50 (see rhos_idx below)
        rhos = std::max(50.0, rhos);

        // best-fit table indecies
        int T_idx = static_cast<int>(std::round((t_soisno(i) - 223) / 5)); // snow aging lookup table temperature index [idx]
        int Tgrd_idx = static_cast<int>(std::round(dTdz / 10)); // snow aging lookup table temperature gradient index [idx]
        int rhos_idx = static_cast<int>(std::round((rhos-50) / 50)); // snow aging lookup table snow density index [idx]

        // boundary check:
        if (T_idx < idx_T_min) { T_idx = idx_T_min; }
        if (T_idx > idx_T_max) { T_idx = idx_T_max; }
        if (Tgrd_idx < idx_Tgrd_min) { Tgrd_idx = idx_Tgrd_min; }
        if (Tgrd_idx > idx_Tgrd_max) { Tgrd_idx = idx_Tgrd_max; }
        if (rhos_idx < idx_rhos_min) { rhos_idx = idx_rhos_min; }
        if (rhos_idx > idx_rhos_max) { rhos_idx = idx_rhos_max; }

        // best-fit parameters
        // retrieved from lookup table
        // snow aging parameter [hour]
        double bst_tau   = snw_table.snowage_tau(T_idx, Tgrd_idx, rhos_idx);
        // snow aging parameter [unitless]
        double bst_kappa = snw_table.snowage_kappa(T_idx, Tgrd_idx, rhos_idx);
        // snow aging parameter [um hr-1]
        double bst_drdt0 = snw_table.snowage_drdt0(T_idx, Tgrd_idx, rhos_idx);

        // change in snow effective radius, using best-fit parameters
        // added checks suggested by mgf. --HW 10/15/2015
        double dr_fresh = snw_rds(i) - SNW_RDS_MIN; // difference between fresh snow r_e and current r_e [um]
        if (std::abs(dr_fresh) < 1.0e-8) {
           dr_fresh = 0.0;
        } else if (dr_fresh < 0.0) {
          throw std::runtime_error("ELM ERROR: SnowAge dr_fresh < 0.0.");
        }

        // incremental change in snow effective radius [um]
        double dr = (bst_drdt0 * std::pow(bst_tau / (dr_fresh+bst_tau), 1.0 / bst_kappa)) * (dtime/3600.0);

        //
        //**********  2. WET SNOW AGING  ***********
        //
        // We are assuming wet and dry evolution occur simultaneously, and
        // the contributions from both can be summed.
        // This is justified by setting the linear offset constant C1_liq_Brun89 to zero [Brun, 1989]

        // fraction of layer mass that is liquid water
        double frc_liq = std::min(0.1, (h2osoi_liq(i) / (h2osoi_liq(i) + h2osoi_ice(i))));

        // incremental change in snow effective radius from wet growth [um]
        double dr_wet = 1.0e18 * (dtime * (C2_liq_Brun89 * std::pow(frc_liq, 3.0)) /
                        (4.0 * ELM_PI * std::pow(snw_rds(i), 2.0)));
        dr += dr_wet;

        //
        //**********  3. SNOWAGE SCALING (TURNED OFF BY DEFAULT)  *************
        //
        // Multiply rate of change of effective radius by some constant, xdrdt
        if (flg_snoage_scl) { dr *= xdrdt; }

        //
        //**********  4. INCREMENT EFFECTIVE RADIUS, ACCOUNTING FOR:  ***********
        //               DRY AGING
        //               WET AGING
        //               FRESH SNOW
        //               RE-FREEZING
        //
        // new snowfall [kg/m2]
        double newsnow;
        if (do_capsnow) {
          newsnow = std::max(0.0, (qflx_snwcp_ice * dtime));
        } else {
          newsnow = std::max(0.0, (qflx_snow_grnd * dtime));
        }

        // snow that has re-frozen [kg/m2]
        double refrzsnow = std::max(0.0, (qflx_snofrz_lyr(i) * dtime));

        // fraction of layer mass that is re-frozen
        double frc_refrz = refrzsnow / h2osno_lyr;

        // fraction of layer mass that is new snow
        double frc_newsnow;
        if (i == snl_top) {
          frc_newsnow = newsnow / h2osno_lyr;
        } else {
          frc_newsnow = 0.0;
        }

        double frc_oldsnow;
        if ((frc_refrz + frc_newsnow) > 1.0) {
          frc_refrz = frc_refrz / (frc_refrz + frc_newsnow);
          frc_newsnow = 1.0 - frc_refrz;
          frc_oldsnow = 0.0;
        } else {
          frc_oldsnow = 1.0 - frc_refrz - frc_newsnow;
        }

        // mass-weighted mean of fresh snow, old snow, and re-frozen snow effective radius
        snw_rds(i) = (snw_rds(i) + dr) * frc_oldsnow + SNW_RDS_MIN * frc_newsnow + snw_rds_refrz * frc_refrz;
        //
        //**********  5. CHECK BOUNDARIES   ***********
        //
        // boundary check
        if (snw_rds(i) < SNW_RDS_MIN) {
           snw_rds(i) = SNW_RDS_MIN;
        }

        if (snw_rds(i) > SNW_RDS_MIN) {
           snw_rds(i) = SNW_RDS_MIN;
        }

        // set top layer variables for history files
        if (i == snl_top) {
          snot_top = t_soisno(i);
          dTdz_top = dTdz;
          snw_rds_top = snw_rds(i);
          sno_liq_top = h2osoi_liq(i) / (h2osoi_liq(i) + h2osoi_ice(i));
        }
      } // for i = snl_top .. snl_btm
    } // if snl > 0

    // Special case: snow on ground, but not enough to have defined a snow layer:
    //   set snw_rds to fresh snow grain size:
    if (snl == 0) {
      if (h2osno > 0.0) { snw_rds(nlevsno-1) = SNW_RDS_MIN; }
    }
  } // snow_aging



/*
!DESCRIPTION:
Evaluate the change of snow mass and the snow water onto soil.
Water flow within snow is computed by an explicit and non-physical
based scheme, which permits a part of liquid water over the holding
capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
percolate into the underlying layer.  Except for cases where the
porosity of one of the two neighboring layers is less than 0.05, zero
flow is assumed. The water flow out of the bottom of the snow pack will
participate as the input of the soil water and runoff.  This subroutine
uses a filter for columns containing snow which must be constructed prior
to being called.


ncells
mflx_neg_snow_col_1d =>  col_wf%mflx_neg_snow_1d , & ! Output:  [real(r8) (:)   ]  mass flux from top soil layer due to negative water content in snow layers (kg H2O /s)

*/
  template<typename ArrayD1>
  ACCELERATE
  void snow_water(const bool& do_capsnow,
                  const int& snl,
                  const double& dtime,
                  const double& frac_sno_eff,
                  const double& h2osno,
                  const double& qflx_sub_snow,
                  const double& qflx_evap_grnd,
                  const double& qflx_dew_snow,
                  const double& qflx_dew_grnd,
                  const double& qflx_rain_grnd,
                  const double& qflx_snomelt,
                  double& qflx_snow_melt,
                  double& qflx_top_soil,
                  double& int_snow,
                  double& frac_sno,
                  double& mflx_neg_snow,
                  ArrayD1 h2osoi_liq,
                  ArrayD1 h2osoi_ice,
                  ArrayD1 mss_bcphi,
                  ArrayD1 mss_bcpho,
                  ArrayD1 mss_dst1,
                  ArrayD1 mss_dst2,
                  ArrayD1 mss_dst3,
                  ArrayD1 mss_dst4,
                  ArrayD1 dz)
  {
    using ELMconst::DENICE;
    using ELMconst::DENH2O;

    // renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
    // surface snow layer resulting from sublimation (frost) / evaporation (condense)
    mflx_neg_snow  = 0.0;

    const int top = nlevsno - snl;
    if (do_capsnow) {
      const double wgdif = h2osoi_ice(top) - frac_sno_eff * qflx_sub_snow * dtime;
      h2osoi_ice(top) = wgdif;
      if (wgdif < 0.0) {
        h2osoi_ice(top) = 0.9;
        h2osoi_liq(top) = h2osoi_liq(top) + wgdif;
      }
      h2osoi_liq(top) = h2osoi_liq(top) - frac_sno_eff * qflx_evap_grnd * dtime;
    } else {
      const double wgdif = h2osoi_ice(top) + frac_sno_eff * (qflx_dew_snow - qflx_sub_snow) * dtime;
      h2osoi_ice(top) = wgdif;
      if (wgdif < 0.0) {
        h2osoi_ice(top) = 0.9;
        h2osoi_liq(top) = h2osoi_liq(top) + wgdif;
      }
      h2osoi_liq(top) = h2osoi_liq(top) + frac_sno_eff * (qflx_rain_grnd + qflx_dew_grnd - qflx_evap_grnd) * dtime;
    }
    // if negative, reduce deeper layer's liquid water content sequentially
    if (h2osoi_liq(top) < 0.0) {
      for (int i = top; i <= nlevsno; ++i) {
        double wgdif = h2osoi_liq(i);
        if (wgdif >= 0.0) break;
        h2osoi_liq(i) = 0.0;
        mflx_neg_snow = wgdif / dtime;
      }
    }

    // porosity and partial volume
    double vol_ice[nlevsno];
    double vol_liq[nlevsno];
    double eff_porosity[nlevsno];
    for (int i = top; i < nlevsno; ++i) {
      // need to scale dz by frac_sno to convert to grid cell average depth
      vol_ice[i]      = std::min(1.0, h2osoi_ice(i) / (dz(i) * frac_sno_eff * DENICE));
      eff_porosity[i] = 1.0 - vol_ice[i];
      vol_liq[i]      = std::min(eff_porosity[i], h2osoi_liq[i] / (dz(i) * frac_sno_eff * DENH2O));
    }
        
    // Capillary forces within snow are usually two or m
    // less than those of gravity. Only gravity terms ar
    // the genernal expression for water flow is "K * ss
    // no effective parameterization for "K".  Thus, a v
    // (not physically based) is introduced:
    // when the liquid water of layer exceeds the layer'
    // capacity, the excess meltwater adds to the underl
    // Also compute aerosol fluxes through snowpack in t
    // 1) compute aerosol mass in each layer
    // 2) add aerosol mass flux from above layer to mass
    // 3) qout_xxx is mass flux of aerosol species xxx o
    //    layer in water flow, proportional to (current)
    //    of aerosol in layer multiplied by a scavenging
    // 4) update mass of aerosol in top layer, according
    // 5) update mass concentration of aerosol according

    double qin = 0.0;
    double qin_bc_phi = 0.0;
    double qin_bc_pho = 0.0;
    double qin_dst1 = 0.0;
    double qin_dst2 = 0.0;
    double qin_dst3 = 0.0;
    double qin_dst4 = 0.0;
    double qout;

    static constexpr double scvng_fct_mlt_bcphi = 0.20; // scavenging factor for hydrophillic BC inclusion in meltwater
    static constexpr double scvng_fct_mlt_bcpho = 0.03; // scavenging factor for hydrophobic BC inclusion in meltwater
    static constexpr double scvng_fct_mlt_dst1  = 0.02; // scavenging factor for dust species 1 inclusion in meltwater
    static constexpr double scvng_fct_mlt_dst2  = 0.02; // scavenging factor for dust species 2 inclusion in meltwater
    static constexpr double scvng_fct_mlt_dst3  = 0.01; // scavenging factor for dust species 3 inclusion in meltwater
    static constexpr double scvng_fct_mlt_dst4  = 0.01; // scavenging factor for dust species 4 inclusion in meltwater
    static constexpr double wimp = 0.05; // Water impremeable if porosity less than wimp
    static constexpr double ssi = 0.033; // irreducible water saturation of snow

    for (int i = top; i < nlevsno; ++i) {

      h2osoi_liq(i) = h2osoi_liq(i) + qin;
      mss_bcphi(i) = mss_bcphi(i) + qin_bc_phi;
      mss_bcpho(i) = mss_bcpho(i) + qin_bc_pho;
      mss_dst1(i)  = mss_dst1(i) + qin_dst1;
      mss_dst2(i)  = mss_dst2(i) + qin_dst2;
      mss_dst3(i)  = mss_dst3(i) + qin_dst3;
      mss_dst4(i)  = mss_dst4(i) + qin_dst4;

      if (i < nlevsno - 1) {
        // no runoff over snow surface, just ponding on surface
        if (eff_porosity[i] < wimp || eff_porosity[i+1] < wimp) {
          qout = 0.0;
        } else {
          // dz must be scaled by frac_sno to obtain gridcell average value
          qout = std::max(0.0, (vol_liq[i] - ssi * eff_porosity[i]) * dz(i) * frac_sno_eff);
          qout = std::min(qout, (1.0 - vol_ice[i+i] - vol_liq[i+1]) * dz(i+1) * frac_sno_eff);
        }
      } else {
        qout = std::max(0.0, (vol_liq[i] - ssi * eff_porosity[i]) * dz(i) * frac_sno_eff);
      }
      qout *= 1000.0;
      h2osoi_liq(i) -= qout;
      qin = qout;


      // mass of ice+water: in extremely rare circumstances, this can
      // be zero, even though there is a snow layer defined. In
      // this case, set the mass to a very small value to
      // prevent division by zero.
      double mss_liqice = h2osoi_liq(i) + h2osoi_ice(i);
      if (mss_liqice < 1.0e-30) {
        mss_liqice = 1.0e-30;
      }


      // BCPHI:
      // 1. flux with meltwater:
      double qout_bc_phi = qout * scvng_fct_mlt_bcphi * (mss_bcphi(i) / mss_liqice);
      if (qout_bc_phi > mss_bcphi(i)) {
        qout_bc_phi = mss_bcphi(i);
      }
      mss_bcphi(i) = mss_bcphi(i) - qout_bc_phi;
      qin_bc_phi = qout_bc_phi;
      
      // BCPHO:
      // 1. flux with meltwater:
      double qout_bc_pho = qout * scvng_fct_mlt_bcpho * (mss_bcpho(i) / mss_liqice);
      if (qout_bc_pho > mss_bcpho(i)) {
        qout_bc_pho = mss_bcpho(i);
      }
      mss_bcpho(i) = mss_bcpho(i) - qout_bc_pho;
      qin_bc_pho = qout_bc_pho;
      
      // DUST 1:
      // 1. flux with meltwater:
      double qout_dst1 = qout * scvng_fct_mlt_dst1 * (mss_dst1(i) / mss_liqice);
      if (qout_dst1 > mss_dst1(i)) {
        qout_dst1 = mss_dst1(i);
      }
      mss_dst1(i) = mss_dst1(i) - qout_dst1;
      qin_dst1 = qout_dst1;
      
      // DUST 2:
      // 1. flux with meltwater:
      double qout_dst2 = qout * scvng_fct_mlt_dst2 * (mss_dst2(i) / mss_liqice);
      if (qout_dst2 > mss_dst2(i)) {
        qout_dst2 = mss_dst2(i);
      }
      mss_dst2(i) = mss_dst2(i) - qout_dst2;
      qin_dst2 = qout_dst2;
      
      // DUST 3:
      // 1. flux with meltwater:
      double qout_dst3 = qout * scvng_fct_mlt_dst3 * (mss_dst3(i) / mss_liqice);
      if (qout_dst3 > mss_dst3(i)) {
        qout_dst3 = mss_dst3(i);
      }
      mss_dst3(i) = mss_dst3(i) - qout_dst3;
      qin_dst3 = qout_dst3;
      
      // DUST 4:
      // 1. flux with meltwater:
      double qout_dst4 = qout * scvng_fct_mlt_dst4 * (mss_dst4(i) / mss_liqice);
      if (qout_dst4 > mss_dst4(i)) {
        qout_dst4 = mss_dst4(i);
      }
      mss_dst4(i) = mss_dst4(i) - qout_dst4;
      qin_dst4 = qout_dst4;
    }

    // Adjust layer thickness for any water+ice content changes in excess of previous 
    // layer thickness. Strictly speaking, only necessary for top snow layer, but doing
    // it for all snow layers will catch problems with older initial files.
    // Layer interfaces (zi) and node depths (z) do not need adjustment here because they
    // are adjusted in CombineSnowLayers and are not used up to that point.
    for (int i = top; i < nlevsno; ++i) {
      dz(i) = std::max(dz(i), h2osoi_liq(i) / DENH2O + h2osoi_ice(i) / DENICE);
    }

    if (snl > 0) {
      // Qout from snow bottom
      qflx_snow_melt += qout / dtime;
      qflx_top_soil = (qout / dtime) + (1.0 - frac_sno_eff) * qflx_rain_grnd;
      int_snow += frac_sno_eff * (qflx_dew_snow + qflx_dew_grnd + qflx_rain_grnd) * dtime;
    } else {
      qflx_snow_melt = qflx_snomelt;
      qflx_top_soil = qflx_rain_grnd + qflx_snomelt;
      // reset accumulated snow when no snow present
      if (h2osno <= 0.0) { int_snow = 0.0; }
      if (h2osno <= 0.0) { frac_sno = 0.0; }
    }
  } // snow_water

  // call ComputeAerosolDeposition between these functions

  // Transfer BC and OC from the within-ice state to the external
  // state based on snow sublimation and re-freezing of liquid water.
  // Re-freezing effect is inactived by default because of
  // uncertainty in how this process operates.
  template <typename ArrayD1>
  ACCELERATE
  void aerosol_phase_change(const int& snl,
                            const double& dtime,
                            const double& qflx_sub_snow,
                            const ArrayD1 h2osoi_liq,
                            const ArrayD1 h2osoi_ice,
                            ArrayD1 mss_bcphi,
                            ArrayD1 mss_bcpho)
  {
    using ELMdims::nlevsno;

    const int top = nlevsno - snl;

    // snow that has sublimated [kg/m2] (top layer only)
    double subsnow = std::max(0.0, (qflx_sub_snow * dtime));
    // fraction of layer mass that has sublimated:
    double frc_sub;
    if ((h2osoi_liq(top) + h2osoi_ice(top)) > 0.0) {
      frc_sub = subsnow / (h2osoi_liq(top) + h2osoi_ice(top));
    } else {
      frc_sub = 0.0;
    }

    for (int i = top; i < nlevsno; ++i) {
      // // snow that has re-frozen [kg/m2]
      // double refrzsnow = std::max(0.0, (qflx_snofrz_lyr(i) * dtime));
      // // fraction of layer mass that is re-frozen
      // double frc_refrz;
      // if ((h2osoi_liq(i) + h2osoi_ice(i)) > 0.0) {
      //   frc_refrz = refrzsnow / (h2osoi_liq(i) + h2osoi_ice(i));
      // } else {
      //   frc_refrz = 0.0;
      // }

      // prohibit sublimation effect to operate on
      // sub-surface layers:
      if (i != top) { frc_sub = 0.0; }

      // fraction of layer mass transformed (sublimation only)
      //double frc_transfer = frc_refrz + frc_sub
      double frc_transfer = frc_sub;
      // cap the fraction at 1
      if (frc_transfer > 1.0) {
        frc_transfer = 1.0;
      }
      // transfer proportionate mass of BC and OC:
      double dm_int  = mss_bcphi(i) * frc_transfer;
      mss_bcphi(i) -= dm_int;
      mss_bcpho(i) += dm_int;
    }
  } // aerosol_phase_change


  template <typename ArrayI1, typename ArrayD1>
  ACCELERATE
  void snow_compaction(const int& snl,
                       const int& subgridflag,
                       const int& ltype,
                       const double& dtime,
                       const double& int_snow,
                       const double& n_melt,
                       const double frac_sno,
                       const ArrayI1 imelt,
                       const ArrayD1 swe_old,
                       const ArrayD1 h2osoi_liq,
                       const ArrayD1 h2osoi_ice,
                       const ArrayD1 t_soisno,
                       const ArrayD1 frac_iceold,
                       ArrayD1 dz)
  {
    using ELMdims::nlevsno;
    using ELMconst::DENH2O;
    using ELMconst::DENICE;
    using ELMconst::TFRZ;
    using ELMconst::ELM_PI;
    using LND::istsoil;
    using LND::istcrop;

    static constexpr double c2 = 23.e-3;   // [m3/kg]
    static constexpr double c3 = 2.777e-6; // [1/s]
    static constexpr double c4 = 0.04;     // [1/K]
    static constexpr double c5 = 2.0;      //
    static constexpr double dm = 100.0;    // Upper Limit on Destructive Metamorphism Compaction [kg/m3]
    static constexpr double eta0 = 9.0e+5;  // The Viscosity Coefficient Eta0 [kg-s/m2]

    const int top = nlevsno - snl;

    double burden = 0.0;
    for (int i = top; i < nlevsno; ++i) {
      double wx = h2osoi_ice(i) + h2osoi_liq(i);
      double vd = 1.0 - (h2osoi_ice(i) / DENICE + h2osoi_liq(i) / DENH2O) / dz(i);
      wx = (h2osoi_ice(i) + h2osoi_liq(i));
      vd = 1.0 - (h2osoi_ice(i) / DENICE + h2osoi_liq(i) / DENH2O) / (frac_sno * dz(i));

      // Allow compaction only for non-saturated node and higher ice lens node.
      if (vd > 0.001 && h2osoi_ice(i) > 0.1) {

        double bi = h2osoi_ice(i) / (frac_sno * dz(i));
        double fi = h2osoi_ice(i) / wx;
        double td = TFRZ - t_soisno(i);
        double dexpf = std::exp(-c4 * td);

        // Settling as a result of destructive metamorphism
        double ddz1 = -c3 * dexpf;
        if (bi > dm) { ddz1 *= std::exp(-46.0e-3 * (bi-dm)); }

        // Liquid water term
        if (h2osoi_liq(i) > 0.01 * dz(i) * frac_sno) { ddz1 *= c5; }

        // Compaction due to overburden
        double ddz2 = -(burden + wx/2.0) * std::exp(-0.08 * td - c2 * bi) / eta0; 
        double ddz3;
        // Compaction occurring during melt
        if (imelt(i) == 1) {
          if (subgridflag == 1 && (ltype == istsoil || ltype == istcrop)) {
            // first term is delta mass over mass
            ddz3 = std::max(0.0, std::min(1.0, (swe_old(i) - wx) / wx));
            // 2nd term is delta fsno over fsno, allowing for negative values for ddz3
            double wsum = 0.0;
            if ((swe_old(i) - wx) > 0.0) {
              if (i == top) {
                for (int j = top; j < nlevsno; ++j) {
                  wsum += h2osoi_liq(j) + h2osoi_ice(j);
                }
              }
              double fsno_melt = 1.0 - std::pow(std::acos(2.0 * std::min(1.0, wsum / int_snow) - 1.0) / ELM_PI, n_melt);
              ddz3 -= std::max(0.0, (fsno_melt - frac_sno) / frac_sno);
            }
            ddz3 = -1.0 / dtime * ddz3;
          } else {
            ddz3 = -1.0 / dtime * std::max(0.0, (frac_iceold(i) - fi) / frac_iceold(i));
          }
        } else {
          ddz3 = 0.0;
        }
        // Time rate of fractional change in dz (units of s-1)
        double pdzdtc = ddz1 + ddz2 + ddz3;

        // The change in dz due to compaction
        // Limit compaction to be no greater than fully saturated layer thickness
        dz(i) = std::max(dz(i) * (1.0 + pdzdtc * dtime), (h2osoi_ice(i) / DENICE + h2osoi_liq(i) / DENH2O) /frac_sno);
      }
      burden += wx;
    }
  }



/*
!DESCRIPTION:
Combine snow layers that are less than a minimum thickness or mass
If the snow element thickness or mass is less than a prescribed minimum,
then it is combined with a neighboring element.  The subroutine
clm\_combo.f90 then executes the combination of mass and energy.
*/
  template <typename ArrayD1>
  ACCELERATE
  void combine_layers(const bool& urbpoi,
                      const int& ltype,
                      const double& dtime,
                      int& snl,
                      double& h2osno,
                      double& snow_depth,
                      double& frac_sno_eff,
                      double& frac_sno,
                      double& int_snow,
                      double& qflx_sl_top_soil,
                      double& qflx_snow2topsoi,
                      double& mflx_snowlyr_col,
                      ArrayD1 t_soisno,
                      ArrayD1 h2osoi_ice,
                      ArrayD1 h2osoi_liq,
                      ArrayD1 snw_rds,
                      ArrayD1 mss_bcphi,
                      ArrayD1 mss_bcpho,
                      ArrayD1 mss_dst1,
                      ArrayD1 mss_dst2,
                      ArrayD1 mss_dst3,
                      ArrayD1 mss_dst4,
                      ArrayD1 dz,
                      ArrayD1 z,
                      ArrayD1 zi)
  {
    using LND::istsoil;
    using LND::istcrop;
    using LND::istwet;
    using LND::istice;
    using LND::istice_mec;

    static constexpr double dzmin[5] = {0.010, 0.015, 0.025, 0.055, 0.115};

    qflx_sl_top_soil = 0.0;
    qflx_snow2topsoi = 0.0;
    mflx_snowlyr_col = 0.0;

    int top_old = nlevsno - snl;
    for (int i = top_old; i < nlevsno; ++i) {
      // use 0.01 to avoid runaway ice buildup
      if (h2osoi_ice(i) <= .01) {
        if (ltype == istsoil || urbpoi || ltype == istcrop) {
          h2osoi_liq(i+1) += h2osoi_liq(i);
          h2osoi_ice(i+1) += h2osoi_ice(i);
          if (i == nlevsno - 1) {
             qflx_sl_top_soil = (h2osoi_liq(i) + h2osoi_ice(i)) / dtime;
             mflx_snowlyr_col += qflx_sl_top_soil;
          }
          if (i != nlevsno - 1) { 

            dz(i+1) += dz(i);

            // NOTE: Temperature, and similarly snw_rds, of the
            // underlying snow layer are NOT adjusted in this case. 
            // Because the layer being eliminated has a small mass, 
            // this should not make a large difference, but it 
            // would be more thorough to do so.
            mss_bcphi(i+1) += mss_bcphi(i);
            mss_bcpho(i+1) += mss_bcpho(i);
            mss_dst1(i+1) += mss_dst1(i);
            mss_dst2(i+1) += mss_dst2(i);
            mss_dst3(i+1) += mss_dst3(i);
            mss_dst4(i+1) += mss_dst4(i);
          }

        } else if (ltype != istsoil && !urbpoi && ltype != istcrop && i != nlevsno - 1) {

            h2osoi_liq(i+1) += h2osoi_liq(i);
            h2osoi_ice(i+1) += h2osoi_ice(i);
            dz(i+1) += dz(i);
            mss_bcphi(i+1) += mss_bcphi(i);
            mss_bcpho(i+1) += mss_bcpho(i);
            mss_dst1(i+1) += mss_dst1(i);
            mss_dst2(i+1) += mss_dst2(i);
            mss_dst3(i+1) += mss_dst3(i);
            mss_dst4(i+1) += mss_dst4(i);
        }

        // shift all elements above this down one.
        int top = nlevsno - snl;
        if (i > top && snl > 1) {
          for (int ii = i; ii > top; --ii) {
            // If the layer closest to the surface is less than 0.1 mm and the ltype is not
            // urban, soil or crop, the h2osoi_liq and h2osoi_ice associated with this layer is sent 
            // to qflx_qrgwl later on in the code.  To keep track of this for the snow balance
            // error check, we add this to qflx_sl_top_soil here
            if (ltype != istsoil && ltype != istcrop &&  !urbpoi && ii == nlevsno - 1) {
              qflx_sl_top_soil = (h2osoi_liq(ii) + h2osoi_ice(ii)) / dtime;
            }
            t_soisno(ii) = t_soisno(ii-1);
            h2osoi_liq(ii) = h2osoi_liq(ii-1);
            h2osoi_ice(ii) = h2osoi_ice(ii-1);
            mss_bcphi(ii) = mss_bcphi(ii-1);
            mss_bcpho(ii) = mss_bcpho(ii-1);
            mss_dst1(ii) = mss_dst1(ii-1);
            mss_dst2(ii) = mss_dst2(ii-1);
            mss_dst3(ii) = mss_dst3(ii-1);
            mss_dst4(ii) = mss_dst4(ii-1);
            snw_rds(ii) = snw_rds(ii-1);
            dz(ii) = dz(ii-1);
          }
        }
        snl -= 1;
      }
    }

    h2osno = 0.0;
    snow_depth = 0.0;
    double zwice = 0.0;
    double zwliq = 0.0;

    top_old = nlevsno - snl;
    for (int i = top_old; i < nlevsno; ++i) {
      h2osno += h2osoi_ice(i) + h2osoi_liq(i);
      snow_depth += dz(i);
      zwice += h2osoi_ice(i);
      zwliq += h2osoi_liq(i);
    }


    // Check the snow depth - all snow gone
    // The liquid water assumes ponding on soil surface.

    if (snow_depth > 0.0 &&
      ((frac_sno_eff * snow_depth < 0.01) || (h2osno / (frac_sno_eff * snow_depth) < 50.0))) {
      snl = 0;
      h2osno = zwice;
      for (int i = 0; i < nlevsno; ++i) {
        mss_bcphi(i) = 0.0;
        mss_bcpho(i) = 0.0;
        mss_dst1(i)  = 0.0;
        mss_dst2(i)  = 0.0;
        mss_dst3(i)  = 0.0;
        mss_dst4(i)  = 0.0;
      }

      if (h2osno <= 0.0) { snow_depth = 0.0; }

      // this is where water is transfered from layer nlevsno - 1 (snow) to layer nlevsno (soil)
      if (ltype == istsoil || urbpoi || ltype == istcrop) {
        h2osoi_liq(nlevsno - 1) = 0.0;
        h2osoi_liq(nlevsno) += zwliq;
        qflx_snow2topsoi = zwliq / dtime;
        mflx_snowlyr_col += zwliq / dtime;
      }
      if (ltype == istwet || ltype == istice || ltype == istice_mec) {
        h2osoi_liq(nlevsno - 1) = 0.0;
      }
    }
    
    if (h2osno <= 0.0) {
      snow_depth = 0.0;
      frac_sno = 0.0;
      frac_sno_eff = 0.0;
      int_snow = 0.0;
    }

    // Check the snow depth - snow layers combined
    // The following loop IS NOT VECTORIZED

    // Two or more layers
    if (snl > 1) {

      int mssi = 0;
      top_old = nlevsno - snl;

      for (int i = top_old; i < nlevsno; ++i) {

        if ((frac_sno_eff * dz(i) < dzmin[mssi]) ||
          ((h2osoi_ice(i) + h2osoi_liq(i)) / (frac_sno_eff * dz(i)) < 50.0)) {

          int neibor;
          if (i == nlevsno - snl) {
             // If top node is removed, combine with bottom neighbor.
             neibor = i + 1;
          } else if (i == nlevsno - 1) {
             // If the bottom neighbor is not snow, combine with the top neighbor.
             neibor = i - 1;
          } else {
             // If none of the above special cases apply, combine with the thinnest neighbor
             neibor = i + 1;
             if ((dz(i-1) + dz(i)) < (dz(i+1) + dz(i))) { neibor = i - 1; }
          }

          // Node l and j are combined and stored as node j.
          int j, l;
          if (neibor > i) {
             j = neibor;
             l = i;
          } else {
             j = i;
             l = neibor;
          }

          // this should be included in 'Combo' for consistency,
          // but functionally it is the same to do it here
          mss_bcphi(j) += mss_bcphi(l);
          mss_bcpho(j) += mss_bcpho(l);
          mss_dst1(j) += mss_dst1(l);
          mss_dst2(j) += mss_dst2(l);
          mss_dst3(j) += mss_dst3(l);
          mss_dst4(j) += mss_dst4(l);


          // mass-weighted combination of effective grain size:
          snw_rds(j) = (snw_rds(j) * (h2osoi_liq(j) + h2osoi_ice(j)) +
                        snw_rds(l) * (h2osoi_liq(l) + h2osoi_ice(l))) /
                       (h2osoi_liq(j) + h2osoi_ice(j) + h2osoi_liq(l) + h2osoi_ice(l));

          combine(dz(l), h2osoi_liq(l), h2osoi_ice(l), t_soisno(l),
                  dz(j), h2osoi_liq(j), h2osoi_ice(j), t_soisno(j));

          // Now shift all elements above this down one.
          if (j-1 > nlevsno - snl) {
            for (int k = j - 1; k > nlevsno - snl - 1; --k) {
              t_soisno(k) = t_soisno(k-1);
              h2osoi_ice(k) = h2osoi_ice(k-1);
              h2osoi_liq(k) = h2osoi_liq(k-1);
              mss_bcphi(k) = mss_bcphi(k-1);
              mss_bcpho(k) = mss_bcpho(k-1);
              mss_dst1(k) = mss_dst1(k-1);
              mss_dst2(k) = mss_dst2(k-1);
              mss_dst3(k) = mss_dst3(k-1);
              mss_dst4(k) = mss_dst4(k-1);
              snw_rds(k) = snw_rds(k-1);
              dz(k) = dz(k-1);
            }
          }

          // Decrease the number of snow layers
          snl -= 1;
          if (snl <= 1) { break; }

        } else {
          // The layer thickness is greater than the prescribed minimum value
          mssi += 1;
        }
      }
    }

    // Reset the node depth and the depth of layer interface
    for (int i = nlevsno - 1; i >= nlevsno - snl; --i) {
      z(i) = zi(i+1) - 0.5 * dz(i);
      zi(i) = zi(i+1) - dz(i);
    }
  }









  template <typename ArrayD1>
  ACCELERATE
  void divide_layers(const double& frac_sno,
                     int& snl,
                     ArrayD1 h2osoi_ice,
                     ArrayD1 h2osoi_liq,
                     ArrayD1 t_soisno,
                     ArrayD1 snw_rds,
                     ArrayD1 mss_bcphi,
                     ArrayD1 mss_bcpho,
                     ArrayD1 mss_dst1,
                     ArrayD1 mss_dst2,
                     ArrayD1 mss_dst3,
                     ArrayD1 mss_dst4,
                     ArrayD1 dz,
                     ArrayD1 z,
                     ArrayD1 zi)
  {
    using ELMdims::nlevsno;
    using ELMconst::TFRZ;
    using snow_snicar::detail::snw_rds_min_tbl;
    using snow_snicar::detail::snw_rds_max_tbl;


    // this is how ELM does it
    // I'm copying it for now
    // will likely change after everything is verified
    double dzsno[nlevsno];
    double swice[nlevsno];
    double swliq[nlevsno];
    double tsno[nlevsno];
    double mbc_phi[nlevsno];
    double mbc_pho[nlevsno];
    double mdst1[nlevsno];
    double mdst2[nlevsno];
    double mdst3[nlevsno];
    double mdst4[nlevsno];
    double rds[nlevsno];

    int msno = snl;
    int top = nlevsno - snl;
    for (int i = 0; i < snl; ++i) {
      dzsno[i] = frac_sno * dz(i+top);
      swice[i] = h2osoi_ice(i+top);
      swliq[i] = h2osoi_liq(i+top);
      tsno[i] = t_soisno(i+top);
      mbc_phi[i] = mss_bcphi(i+top);
      mbc_pho[i] = mss_bcpho(i+top);
      mdst1[i] = mss_dst1(i+top);
      mdst2[i] = mss_dst2(i+top);
      mdst3[i] = mss_dst3(i+top);
      mdst4[i] = mss_dst4(i+top);
      rds[i] = snw_rds(i+top);
    }

    if (msno == 1) {
      if (dzsno[0] > 0.03) {
        msno = 2;
        dzsno[0] /= 2.0;
        swice[0] /= 2.0;
        swliq[0] /= 2.0;
        dzsno[1] = dzsno[0];
        swice[1] = swice[0];
        swliq[1] = swliq[0];
        tsno[1] = tsno[0];
        mbc_phi[0] /= 2.0;
        mbc_phi[1] = mbc_phi[0];
        mbc_pho[0] /= 2.0;
        mbc_pho[1] = mbc_pho[0];
        mdst1[0] /= 2.0;
        mdst1[1] = mdst1[0];
        mdst2[0] /= 2.0;
        mdst2[1] = mdst2[0];
        mdst3[0] /= 2.0;
        mdst3[1] = mdst3[0];
        mdst4[0] /= 2.0;
        mdst4[1] = mdst4[0];
        rds[1] = rds[0];
      }
    }

    if (msno > 1) {
      if (dzsno[0] > 0.02) {
        
        double drr = dzsno[0] - 0.02;
        double propor = drr / dzsno[0];
        double zwice = propor * swice[0];
        double zwliq = propor * swliq[0];
        double zmbc_phi = propor * mbc_phi[0];
        double zmbc_pho = propor * mbc_pho[0];
        double zmdst1 = propor * mdst1[0];
        double zmdst2 = propor * mdst2[0];
        double zmdst3 = propor * mdst3[0];
        double zmdst4 = propor * mdst4[0];

        propor = 0.02 / dzsno[0];
        swice[0] *= propor;
        swliq[0] *= propor;
        mbc_phi[0] *= propor;
        mbc_pho[0] *= propor;
        mdst1[0] *= propor;
        mdst2[0] *= propor;
        mdst3[0] *= propor;
        mdst4[0] *= propor;
        dzsno[0] = 0.02;
        mbc_phi[1] += zmbc_phi;
        mbc_pho[1] += zmbc_pho;
        mdst1[1] += zmdst1;
        mdst2[1] += zmdst2;
        mdst3[1] += zmdst3;
        mdst4[1] += zmdst4;

        rds[1] = (rds[1] * (swliq[1] + swice[1]) + rds[0] * (zwliq + zwice)) / (swliq[1] + swice[1] + zwliq + zwice);

        if (rds[1] < snw_rds_min_tbl || rds[1] > snw_rds_max_tbl) {
          throw std::runtime_error("ELM ERROR: snow radius out of bounds in snow::divide_layers.");
          // calculate simply?
          // rds[1] = rds[0];
        }

        combine(drr, zwliq, zwice, tsno[0], dzsno[1], swliq[1], swice[1], tsno[1]);

        // Subdivide a new layer
        if (msno <= 2 && dzsno[1] > 0.07) {
          msno = 3;
          double dtdz = (tsno[0] - tsno[1]) / ((dzsno[0] + dzsno[1]) / 2.0);
          dzsno[1] /= 2.0;
          swice[1] /= 2.0;
          swliq[1] /= 2.0;
          dzsno[2] = dzsno[1];
          swice[2] = swice[1];
          swliq[2] = swliq[1];
          tsno[2] = tsno[1] - dtdz * dzsno[1] / 2.0;

          if (tsno[2] >= TFRZ) {
            tsno[2] = tsno[1];
          } else {
            tsno[1] += dtdz * dzsno[1] / 2.0;
          }

          mbc_phi[1] /= 2.0;
          mbc_phi[2] = mbc_phi[1];
          mbc_pho[1] /= 2.0;
          mbc_pho[2] = mbc_pho[1];
          mdst1[1] /= 2.0;
          mdst1[2] = mdst1[1];
          mdst2[1] /= 2.0;
          mdst2[2] = mdst2[1];
          mdst3[1] /= 2.0;
          mdst3[2] = mdst3[1];
          mdst4[1] /= 2.0;
          mdst4[2] = mdst4[1];
          rds[2] = rds[1];
        }
      }
    } // if (msno > 1)

    if (msno > 2) {
      if (dzsno[1] > 0.05) {
        double drr = dzsno[1] - 0.05;
        double propor = drr / dzsno[1];
        double zwice = propor * swice[1];
        double zwliq = propor * swliq[1];
        double zmbc_phi = propor * mbc_phi[1];
        double zmbc_pho = propor * mbc_pho[1];
        double zmdst1 = propor * mdst1[1];
        double zmdst2 = propor * mdst2[1];
        double zmdst3 = propor * mdst3[1];
        double zmdst4 = propor * mdst4[1];

        propor = 0.05 / dzsno[1];

        swice[1] *= propor;
        swliq[1] *= propor;
        mbc_phi[1] *= propor;
        mbc_pho[1] *= propor;
        mdst1[1] *= propor;
        mdst2[1] *= propor;
        mdst3[1] *= propor;
        mdst4[1] *= propor;

        dzsno[1] = 0.05;

        mbc_phi[2] += zmbc_phi;
        mbc_pho[2] += zmbc_pho;
        mdst1[2] += zmdst1;
        mdst2[2] += zmdst2;
        mdst3[2] += zmdst3;
        mdst4[2] += zmdst4;

        rds[2] = (rds[2] * (swliq[2] + swice[2]) + rds[1] * (zwliq + zwice)) / (swliq[2] + swice[2] + zwliq + zwice);

        if (rds[2] < snw_rds_min_tbl || rds[2] > snw_rds_max_tbl) {
          throw std::runtime_error("ELM ERROR: snow radius out of bounds in snow::divide_layers.");
          // calculate simply?
          // rds[2] = rds[1];
        }

        combine(drr, zwliq, zwice, tsno[1], dzsno[2], swliq[2], swice[2], tsno[2]);

        // Subdivided a new layer
        if (msno <= 3 && dzsno[2] > 0.18) {
          msno = 4;
          double dtdz = (tsno[1] - tsno[2]) / ((dzsno[1] + dzsno[2]) / 2.0) ;
          dzsno[2] /= 2.0;
          swice[2] /= 2.0;
          swliq[2] /= 2.0;
          dzsno[3] = dzsno[2];
          swice[3] = swice[2];
          swliq[3] = swliq[2];
          tsno[3] = tsno[2] - dtdz * dzsno[2] / 2.0;
          if (tsno[2] >= TFRZ) {
            tsno[3] = tsno[2];
          } else {
            tsno[2] += dtdz * dzsno[2] / 2.0;
          }
          mbc_phi[2] /= 2.0;
          mbc_phi[3] = mbc_phi[2];
          mbc_pho[2] /= 2.0;
          mbc_pho[3] = mbc_pho[2];
          mdst1[2] /= 2.0;
          mdst1[3] = mdst1[2];
          mdst2[2] /= 2.0;
          mdst2[3] = mdst2[2];
          mdst3[2] /= 2.0;
          mdst3[3] = mdst3[2];
          mdst4[2] /= 2.0;
          mdst4[3] = mdst4[2];
          rds[3] = rds[2];
        }
      }
    } // if (msno > 2)

    if (msno > 3) {
      if (dzsno[2] > 0.11) {
        double drr = dzsno[2] - 0.11;
        double propor = drr / dzsno[2];
        double zwice = propor * swice[2];
        double zwliq = propor * swliq[2];
        double zmbc_phi = propor * mbc_phi[2];
        double zmbc_pho = propor * mbc_pho[2];
        double zmdst1 = propor * mdst1[2];
        double zmdst2 = propor * mdst2[2];
        double zmdst3 = propor * mdst3[2];
        double zmdst4 = propor * mdst4[2];

        propor = 0.11 / dzsno[2];

        swice[2] *= propor;
        swliq[2] *= propor;
        mbc_phi[2] *= propor;
        mbc_pho[2] *= propor;
        mdst1[2] *= propor;
        mdst2[2] *= propor;
        mdst3[2] *= propor;
        mdst4[2] *= propor;

        dzsno[2] = 0.11;

        mbc_phi[3] += zmbc_phi;
        mbc_pho[3] += zmbc_pho;
        mdst1[3] += zmdst1;
        mdst2[3] += zmdst2;
        mdst3[3] += zmdst3;
        mdst4[3] += zmdst4;

        rds[3] = (rds[3] * (swliq[3] + swice[3]) + rds[2] * (zwliq + zwice)) / (swliq[3] + swice[3] + zwliq + zwice);

        if (rds[3] < snw_rds_min_tbl || rds[3] > snw_rds_max_tbl) {
          throw std::runtime_error("ELM ERROR: snow radius out of bounds in snow::divide_layers.");
          // calculate simply?
          // rds[3] = rds[2];
        }

        combine(drr, zwliq, zwice, tsno[2], dzsno[3], swliq[3], swice[3], tsno[3]);

        // Subdivided a new layer
        if (msno <= 4 && dzsno[3] > 0.41) {
          msno = 5;
          double dtdz = (tsno[2] - tsno[3])/((dzsno[2]+dzsno[3]) / 2.0);
          dzsno[3] /= 2.0;
          swice[3] /= 2.0;
          swliq[3] /= 2.0;
          dzsno[4] = dzsno[3];
          swice[4] = swice[3];
          swliq[4] = swliq[3];
          tsno[4] = tsno[3] - dtdz * dzsno[3] / 2.0 ;
          if (tsno[4] >= TFRZ) {
            tsno[4] = tsno[3];
          } else {
            tsno[3] += dtdz * dzsno[3] / 2.0;
          }
          mbc_phi[3] /= 2.0;
          mbc_phi[4] = mbc_phi[3];
          mbc_pho[3] /= 2.0;
          mbc_pho[4] = mbc_pho[3];
          mdst1[3] /= 2.0;
          mdst1[4] = mdst1[3];
          mdst2[3] /= 2.0;
          mdst2[4] = mdst2[3];
          mdst3[3] /= 2.0;
          mdst3[4] = mdst3[3];
          mdst4[3] /= 2.0;
          mdst4[4] = mdst4[3];
          rds[4] = rds[3];
        }
      }
    } // if (msno > 3)

    if (msno > 4) {
      if (dzsno[3] > 0.23) {
        double drr = dzsno[3] - 0.23;
        double propor = drr/dzsno[3];
        double zwice = propor * swice[3];
        double zwliq = propor * swliq[3];

        double zmbc_phi = propor * mbc_phi[3];
        double zmbc_pho = propor * mbc_pho[3];
        double zmdst1 = propor * mdst1[3];
        double zmdst2 = propor * mdst2[3];
        double zmdst3 = propor * mdst3[3];
        double zmdst4 = propor * mdst4[3];
        
        propor = 0.23 / dzsno[3];
        
        swice[3] *= propor;
        swliq[3] *= propor;
        mbc_phi[3] *= propor;
        mbc_pho[3] *= propor;
        mdst1[3] *= propor;
        mdst2[3] *= propor;
        mdst3[3] *= propor;
        mdst4[3] *= propor;

        dzsno[3] = 0.23;

        mbc_phi[4] += zmbc_phi;
        mbc_pho[4] += zmbc_pho;
        mdst1[4] += zmdst1;
        mdst2[4] += zmdst2;
        mdst3[4] += zmdst3;
        mdst4[4] += zmdst4;

        rds[4] = (rds[4] * (swliq[4] + swice[4]) + rds[3] * (zwliq + zwice)) / (swliq[4] + swice[4] + zwliq + zwice);

        if (rds[3] < snw_rds_min_tbl || rds[3] > snw_rds_max_tbl) {
          throw std::runtime_error("ELM ERROR: snow radius out of bounds in snow::divide_layers.");
          // calculate simply?
          // rds[4] = rds[3];
        }

        combine(drr, zwliq, zwice, tsno[3], dzsno[4], swliq[4], swice[4], tsno[4]);

      }
    } // if (msno > 4)

    // place into state variables
    snl = msno;
    top = nlevsno - snl;
    for (int i = top; i < nlevsno; ++i) {
      dz(i) = dzsno[i-top] / frac_sno;
      h2osoi_ice(i) = swice[i-top];
      h2osoi_liq(i) = swliq[i-top];
      t_soisno(i) = tsno[i-top];
      mss_bcphi(i) = mbc_phi[i-top];
      mss_bcpho(i) = mbc_pho[i-top];
      mss_dst1(i) = mdst1[i-top];
      mss_dst2(i) = mdst2[i-top];
      mss_dst3(i) = mdst3[i-top];
      mss_dst4(i) = mdst4[i-top];
      snw_rds(i) = rds[i-top];
    }

    // adjust node depths and interfaces
    for (int i = nlevsno - 1; i >= top; --i) {
      z(i) = zi(i+1) - 0.5 * dz(i);
      zi(i) = zi(i+1) - dz(i);
    }
  } // snow::divide_layers










// dz2    nodal thickness of element being absorbed [m]
// wliq2  liquid water of element 2 [kg/m2]
// wice2  ice of element 2 [kg/m2]
// t2     nodal temperature of element 2 [K]
// dz     nodal thickness of absorbing element [m]
// wliq   liquid water of element 1
// wice   ice of element 1 [kg/m2]
// t      nodel temperature of elment 1 [K]
  ACCELERATE
  void combine(const double& dz2,
               const double& wliq2,
               const double& wice2,
               const double& t2,
               double& dz,
               double& wliq,
               double& wice,
               double& t)
  {
    using ELMconst::CPICE;
    using ELMconst::CPWAT;
    using ELMconst::TFRZ;
    using ELMconst::HFUS;

    double h = (CPICE * wice + CPWAT * wliq) * (t - TFRZ) + HFUS * wliq;
    double h2 = (CPICE * wice2 + CPWAT * wliq2) * (t2 - TFRZ) + HFUS * wliq2;
    wice += wice2;
    wliq += wliq2;
    double tc = TFRZ + (h + h2 - HFUS * wliq) / (CPICE * wice + CPWAT * wliq);
    dz += dz2;
    t = tc;
  }


  // Set empty snow layers to zero
  template <typename ArrayD1>
  ACCELERATE
  void prune_snow_layers(const int& snl,
                         ArrayD1 h2osoi_ice,
                         ArrayD1 h2osoi_liq,
                         ArrayD1 t_soisno,
                         ArrayD1 dz,
                         ArrayD1 z,
                         ArrayD1 zi)
  {
    const int top = nlevsno - snl;
    for (int i = 0; i < top; ++i) {
      h2osoi_ice(i) = 0.0;
      h2osoi_liq(i) = 0.0;
      t_soisno(i)  = 0.0;
      dz(i)    = 0.0;
      z(i)     = 0.0;
      zi(i)  = 0.0;
    }
  }

} // namespace ELM::snow
