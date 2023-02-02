
#pragma once

#include <cmath>

namespace ELM::soil_temp {

// h2osoi_ice_sl1, t_soisno_sl1, and fact_sl1 are all from the
// first snow layer (closest to ground) (nlevsno() - 1)
// they are single elements of larger column arrays and the
// driver code is responsible for passing in the correct values
ACCELERATE
void phase_change_h2osfc(const int& snl,
                         const double& dtime,
                         const double& frac_sno,
                         const double& frac_h2osfc,
                         const double& dhsdT,
                         const double& c_h2osfc,
                         const double& fact_sl1,
                         double& t_h2osfc,
                         double& h2osfc,
                         double& xmf_h2osfc,
                         double& qflx_h2osfc_to_ice,
                         double& eflx_h2osfc_to_snow,
                         double& h2osno,
                         double& int_snow,
                         double& snow_depth,
                         double& h2osoi_ice_sl1,
                         double& t_soisno_sl1)
{
  using ELMconst::TFRZ;
  using ELMconst::HFUS;
  using ELMconst::DENICE;
  using ELMconst::CPWAT;

  qflx_h2osfc_to_ice = 0.0;
  eflx_h2osfc_to_snow = 0.0;
  xmf_h2osfc = 0.0;

  // if liquid exists below melt point, freeze some to ice.
  if (frac_h2osfc > 0.0 && t_h2osfc <= TFRZ()) {

    const double tinc = TFRZ() - t_h2osfc;
    t_h2osfc = TFRZ();
    // energy absorbed beyond freezing temperature
    double hm = frac_h2osfc * (dhsdT * tinc - tinc * c_h2osfc / dtime);
    // mass of water converted from liquid to ice
    double xm = hm * dtime / HFUS();
    const double temp1 = h2osfc + xm;
    const double z_avg = frac_sno * snow_depth;
    double rho_avg;
    if (z_avg > 0.0) {
      rho_avg = std::min(800.0, h2osno / z_avg);
    } else {
      rho_avg = 200.0;
    }

    //=====================  xm < h2osfc  ====================================
    if (temp1 >= 0.0) { //add some frozen water to snow column
      // add ice to snow column
      h2osno -= xm;
      int_snow -= xm;
      if (snl > 0) h2osoi_ice_sl1 -= xm;
      // remove ice from h2osfc
      h2osfc += xm;
      xmf_h2osfc = hm;
      qflx_h2osfc_to_ice = -xm / dtime;
      // update snow depth
      if (frac_sno > 0 && snl > 0) {
        snow_depth = h2osno / (rho_avg * frac_sno);
      } else {
        snow_depth = h2osno / DENICE();
      }
      // adjust temperature of lowest snow layer to account for addition of ice
      if (snl == 0) {
        //initialize for next time step
        t_soisno_sl1 = t_h2osfc;
        eflx_h2osfc_to_snow = 0.0;
      } else {
        double c1, c2;
        if (snl == 1) {
          c1 = frac_sno * (dtime / fact_sl1 - dhsdT * dtime);
        } else {
          c1 = frac_sno / fact_sl1 * dtime;
        }
        if (frac_h2osfc != 0.0) {
          c2 = (-CPWAT() * xm - frac_h2osfc * dhsdT * dtime);
        } else {
          c2 = 0.0;
        }
        t_soisno_sl1 = (c1 * t_soisno_sl1 + c2 * t_h2osfc) / (c1 + c2);
        eflx_h2osfc_to_snow = (t_h2osfc - t_soisno_sl1) * c2 / dtime;
      }
       //=========================  xm > h2osfc  =============================
    } else { //all h2osfc converted to ice

      rho_avg = (h2osno * rho_avg + h2osfc * DENICE()) / (h2osno + h2osfc);
      h2osno += h2osfc;
      int_snow += h2osfc;

      qflx_h2osfc_to_ice = h2osfc / dtime;

      // excess energy is used to cool ice layer
      if (snl > 0) h2osoi_ice_sl1 = h2osoi_ice_sl1 + h2osfc;

      // NOTE: should compute and then use the heat capacity of frozen h2osfc layer
      //       rather than using heat capacity of the liquid layer. But this causes
      //       balance check errors as it doesn't know about it.
      // compute heat capacity of frozen h2osfc layer

      // cool frozen h2osfc layer with extra heat
      t_h2osfc = t_h2osfc - temp1 * HFUS() / (dtime * dhsdT - c_h2osfc);

      xmf_h2osfc = hm - frac_h2osfc * temp1 * HFUS() / dtime;

      // next, determine equilibrium temperature of combined ice/snow layer
      double c1, c2;
      if (snl == 0) {
        // initialize for next time step
        t_soisno_sl1 = t_h2osfc;
      } else if (snl == 1) {
        c1 = frac_sno * (dtime / fact_sl1 - dhsdT * dtime);
        if (frac_h2osfc != 0.0) {
          c2 = frac_h2osfc * (c_h2osfc - dtime * dhsdT);
        } else {
          c2 = 0.0;
        }
        // account for the change in t_soisno(c,0) via xmf_h2osfc(c)
        t_soisno_sl1 = (c1 * t_soisno_sl1 + c2 * t_h2osfc) / (c1 + c2);
        t_h2osfc = t_soisno_sl1;
      } else {
         c1 = frac_sno / fact_sl1 * dtime;
         if (frac_h2osfc != 0.0) {
            c2 = frac_h2osfc * (c_h2osfc - dtime * dhsdT);
         } else {
            c2=0.0;
         }
         t_soisno_sl1 = (c1 * t_soisno_sl1 + c2 * t_h2osfc) / (c1 + c2);
         t_h2osfc = t_soisno_sl1;
      }

      // set h2osfc to zero (all liquid converted to ice)
      h2osfc = 0.0;

      // update snow depth
      if (frac_sno > 0.0 && snl > 0) { 
        snow_depth = h2osno / (rho_avg * frac_sno);
      } else {
        snow_depth = h2osno / DENICE();
      }
    }
  }
}





/*
!DESCRIPTION:
Calculation of the phase change within snow and soil layers:
(1) Check the conditions for which the phase change may take place,
    i.e., the layer temperature is great than the freezing point
    and the ice mass is not equal to zero (i.e. melting),
    or the layer temperature is less than the freezing point
    and the liquid water mass is greater than the allowable supercooled 
    liquid water calculated from freezing point depression (i.e. freezing).
(2) Assess the rate of phase change from the energy excess (or deficit)
    after setting the layer temperature to freezing point.
(3) Re-adjust the ice and liquid mass, and the layer temperature



qflx_snomelt                snow melt (mm H2O /s)
qflx_snow_melt              snow melt (net)
qflx_snofrz_lyr[nlevsno()]    snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1] 
qflx_snofrz                 column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]

xmf                     (:)   => null() ! total latent heat of phase change of ground water
xmf_h2osfc              (:)   => null() ! latent heat of phase change of surface water
imelt                   (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno()+1:nlevgrnd()) 
*/

template<typename ArrayI1, typename ArrayD1>
ACCELERATE
void phase_change_soisno(const int& snl,
                         const int& ltype,
                         const double& dtime,
                         const double& dhsdT,
                         const double& frac_h2osfc,
                         const double& frac_sno_eff,
                         const ArrayD1 fact,
                         const ArrayD1 watsat,
                         const ArrayD1 sucsat,
                         const ArrayD1 bsw,
                         const ArrayD1 dz,
                         double& h2osno,
                         double& snow_depth,
                         double& xmf,
                         double& qflx_snofrz,
                         double& qflx_snow_melt,
                         double& qflx_snomelt,
                         double& eflx_snomelt,
                         ArrayI1 imelt,
                         ArrayD1 qflx_snofrz_lyr,
                         ArrayD1 h2osoi_ice,
                         ArrayD1 h2osoi_liq,
                         ArrayD1 t_soisno)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevgrnd;
  using ELMconst::TFRZ;
  using ELMconst::HFUS;
  using ELMconst::GRAV;

  using LND::istsoil;
  using LND::istcrop;
  using LND::icol_road_perv;

  // initialization
  xmf = 0.0;
  qflx_snofrz = 0.0;
  qflx_snow_melt = 0.0;
  qflx_snomelt = 0.0;
  for (int i = 0; i < nlevsno(); ++i) { qflx_snofrz_lyr(i) = 0.0; }

  const int top = nlevsno() - snl;
  for (int i = top; i < nlevsno() + nlevgrnd(); ++i) { imelt(i) = 0; }

  // local arrays
  double tinc[nlevgrnd()+nlevsno()];
  double supercool[nlevgrnd()]; 

  //--  snow layers  --------------------------------------------------- 
  for (int i = top; i < nlevsno(); ++i) {
    // melting identification
    // if ice exists above melt point, melt some to liquid.
    if (h2osoi_ice(i) > 0.0 && t_soisno(i) > TFRZ()) {
      imelt(i) = 1;
      tinc[i] = TFRZ() - t_soisno(i); 
      t_soisno(i) = TFRZ();
    }

    // Freezing identification
    // If liquid exists below melt point, freeze some to ice.
    if (h2osoi_liq(i) > 0.0 && t_soisno(i) < TFRZ()) {
      imelt(i) = 2;
      tinc[i] = TFRZ() - t_soisno(i);
      t_soisno(i) = TFRZ();
    }
  }

  //-- soil layers   ---------------------------------------------------
  for (int i = nlevsno(); i < nlevsno() + nlevgrnd(); ++i) {

    if (h2osoi_ice(i) > 0.0 && t_soisno(i) > TFRZ()) {
      imelt(i) = 1;
      tinc[i] = TFRZ() - t_soisno(i);
      t_soisno(i) = TFRZ();
    }

    // from Zhao (1997) and Koren (1999)
    supercool[i-nlevsno()] = 0.0;
    if (ltype == istsoil || ltype == istcrop || ltype == icol_road_perv) {
      if(t_soisno(i) < TFRZ()) {
        double smp = HFUS() * (TFRZ() - t_soisno(i)) / (GRAV() * t_soisno(i)) * 1000.0; // (mm)
        supercool[i-nlevsno()] = watsat(i-nlevsno()) * pow(smp / sucsat(i-nlevsno()), -1.0 / bsw(i-nlevsno()));
        supercool[i-nlevsno()] *= dz(i) * 1000.0; // (mm)
      }
    }

    if (h2osoi_liq(i) > supercool[i-nlevsno()] && t_soisno(i) < TFRZ()) {
      imelt(i) = 2;
      tinc[i] = TFRZ() - t_soisno(i);
      t_soisno(i) = TFRZ();
    }

    // if snow exists, but its thickness is less than the critical value (0.01 m)
    if (snl == 0 && h2osno > 0.0 && i == nlevsno()) {
      if (t_soisno(i) > TFRZ()) {
        imelt(i) = 1;
        tinc[i] = TFRZ() - t_soisno(i);
        t_soisno(i) = TFRZ();
      }
    }
  }


  //-- all layers   ---------------------------------------------------
  for (int i = top; i < nlevsno() + nlevgrnd(); ++i) {
    
    double hm = 0.0;
    // Calculate the energy surplus and loss for melting and freezing
    if (imelt(i) > 0) {
      // this nest of logic defines the following scenarios:
      // this snow layer is the top active layer
      // this soil layer is the top active layer
      // this soil with surface water layer is the top active layer
      // this soil layer shares an interface with snow
      // interior soil layer
      // interior snow layer
      //=================================================================
      if (i == top) { // top active layer
        if (i < nlevsno()) { // top layer is snow
           hm = frac_sno_eff * (dhsdT * tinc[i] - tinc[i] / fact(i));
        } else {
          double temp_hm = dhsdT * tinc[i] - tinc[i] / fact(i); // top layer is soil or soil with surface water
          hm = (frac_h2osfc != 0.0) ? temp_hm - frac_h2osfc * (dhsdT * tinc[i]) : temp_hm;
        }
      } else if (i == nlevsno()) { // top soil layer under snow
        hm = (1.0 - frac_sno_eff - frac_h2osfc) * dhsdT * tinc[i] - tinc[i] / fact(i);
      } else { // non-interfacial layers
        if (i < nlevsno()) { // snow
          hm = -frac_sno_eff * (tinc[i] / fact(i));
        } else { // soil
          hm = -tinc[i] / fact(i);
        }
      }
    }

    // These two errors were checked carefully (Y. Dai).  They result from the
    // computed error of "Tridiagonal-Matrix" in subroutine "thermal".
    if (imelt(i) == 1 && hm < 0.0) {
      hm = 0.0;
      imelt(i) = 0;
    }
    if (imelt(i) == 2 && hm > 0.0) {
      hm = 0.0;
      imelt(i) = 0;
    }

    // The rate of melting and freezing
    if (imelt(i) > 0 && std::abs(hm) > 0.0) {
      double xm = hm * dtime / HFUS(); // kg/m2
      // If snow exists, but its thickness is less than the critical value
      // (1 cm). Note: more work is needed to determine how to tune the
      // snow depth for this case
      if (i == nlevsno()) {
        if (snl == 0 && h2osno > 0.0 && xm > 0.0) {
          double temp1 = h2osno; // kg/m2
          h2osno = std::max(0.0, temp1 - xm);
          double propor = h2osno / temp1;
          snow_depth *= propor;
          double heatr = hm - HFUS() * (temp1 - h2osno) / dtime; // W/m2
          if (heatr > 0.0) {
            xm = heatr * dtime / HFUS(); // kg/m2
            hm = heatr; // W/m2
          } else {
            xm = 0.0;
            hm = 0.0;
          }
          qflx_snomelt = std::max(0.0, temp1 - h2osno) / dtime; // kg/(m2 s)
          xmf = HFUS() * qflx_snomelt;
          qflx_snow_melt = qflx_snomelt;
        }
      }

      double heatr = 0.0;
      double wmass0 = h2osoi_ice(i) + h2osoi_liq(i);
      double wice0 = h2osoi_ice(i);
      if (xm > 0.0) {
        h2osoi_ice(i) = std::max(0.0, wice0 - xm);
        heatr = hm - HFUS() * (wice0 - h2osoi_ice(i)) / dtime;
      } else if (xm < 0.0) {
        if (i < nlevsno()) {
          h2osoi_ice(i) = std::min(wmass0, wice0 - xm); // snow
        } else {
          if (wmass0 < supercool[i-nlevsno()]) {
            h2osoi_ice(i) = 0.0;
          } else {
            h2osoi_ice(i) = std::min(wmass0 - supercool[i-nlevsno()], wice0 - xm);
          }
        }
        heatr = hm - HFUS() * (wice0 - h2osoi_ice(i)) / dtime;
      }

      h2osoi_liq(i) = std::max(0.0, wmass0 - h2osoi_ice(i));
      if (std::abs(heatr) > 0.0) {
        if (i == top) {
          if (snl == 0) {
            t_soisno(i) += fact(i) * heatr / (1.0 - (1.0 - frac_h2osfc) * fact(i) * dhsdT);
          } else {
            t_soisno(i) += (fact(i) / frac_sno_eff) * heatr / (1.0 - fact(i) * dhsdT);
          }
        } else if (i == nlevsno()) {
          t_soisno(i) += fact(i) * heatr / (1.0 - (1.0 - frac_sno_eff - frac_h2osfc) * fact(i) * dhsdT);
        } else {
          if (i >= nlevsno()) {
            t_soisno(i) += fact(i) * heatr;
          } else {
            if (frac_sno_eff > 0.0) t_soisno(i) += (fact(i) / frac_sno_eff) * heatr;
          }
        }
        if (i < nlevsno()) { // snow
          if (h2osoi_liq(i) * h2osoi_ice(i) > 0.0) { t_soisno(i) = TFRZ(); }
        }
      } // end of heatr > 0 if-block

      xmf += HFUS() * (wice0 - h2osoi_ice(i)) / dtime;
      if (imelt(i) == 1 && i < nlevsno()) {
        qflx_snomelt += std::max(0.0, (wice0 - h2osoi_ice(i))) / dtime;
      }
      // layer freezing mass flux (positive):
      if (imelt(i) == 2 && i < nlevsno()) {
        qflx_snofrz_lyr(i) = std::max(0.0, (h2osoi_ice(i) - wice0)) / dtime;
      }
    } // if (imelt(i) > 0 && abs(hm) > 0.0)
  } // all layers loop


  eflx_snomelt = qflx_snomelt * HFUS();
  for (int i = 0; i < nlevsno(); ++i) {
    if (imelt(i) == 2 && i < nlevsno()) {
      qflx_snofrz += qflx_snofrz_lyr(i);
    }
  }
}


template<typename ArrayD1>
ACCELERATE
void phase_change_correction(const int& snl,
                             const ArrayD1 tk,
                             const ArrayD1 t_soisno,
                             const ArrayD1 z,
                             ArrayD1 fn1)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevgrnd;

  const int top = nlevsno() - snl;
  for (int i = top; i < nlevgrnd() + nlevsno() - 1; ++i) {
    fn1(i) = tk(i) * (t_soisno(i+1) - t_soisno(i)) / (z(i+1)-z(i));
  }
  fn1(nlevgrnd() + nlevsno() - 1) = 0.0;
}


} // namespace ELM::soil_temp
