
#pragma once

namespace ELM::snow {



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

      if (i <= nlevsno + 2) {
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

} // namespace ELM::snow
