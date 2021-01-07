/* functions derived from CanopyFluxesMod.F90

*/
#include <algorithm>
#include <cmath>
#include <assert>
#include "clm_constants.h"


class CanopyFluxes {

private:
  double btran0 = 0.0;
  bool found = false;





public:


void InitializeFlux()
// need to decide where to place photosyns_vars%TimeStepInit() calc_effective_soilporosity() calc_volumetric_h2oliq()
// calc_root_moist_stress()
// need to add irrigation function, includes 592-654, 446, maybe more
{
  irrig_nsteps_per_day = ((irrig_length + (dtime - 1))/dtime);  // round up

  // -----------------------------------------------------------------
  // Time step initialization of photosynthesis variables
  // -----------------------------------------------------------------
  //call photosyns_vars%TimeStepInit(bounds)

  if (!Land.lakpoi && !Land.urbpoi) {
    if (frac_veg_nosno == 0) {
      btran = 0.0;
      t_veg = forc_t;
      cf_bare  = forc_pbot / (ELM_RGAS * 0.001 * thm) * 1.e06;
      rssun = 1.0 / 1.e15 * cf_bare;
      rssha = 1.0 / 1.e15 * cf_bare;
      lbl_rsc_h2o = 0.0;
      for (int i = 0; i < nlevgrnd; i++) {
        rootr[i]  = 0.0;
        rresis[i] = 0.0;
      }
    } else {
      del(p)    = 0.0;  // change in leaf temperature from previous iteration
      efeb(p)   = 0.0;  // latent head flux from leaf for previous iteration
      wtlq0(p)  = 0.0;
      wtalq(p)  = 0.0;
      wtgq(p)   = 0.0;
      wtaq0(p)  = 0.0;
      obuold(p) = 0.0;
      btran(p)  = btran0;

      // calculate dayl_factor as the ratio of (current:max dayl)^2
      // set a minimum of 0.01 (1%) for the dayl_factor
      dayl_factor = std::min(1.0, std::max(0.01, (dayl * dayl) / (max_dayl * max_dayl)));

      // compute effective soil porosity
      calc_effective_soilporosity(watsat, h2osoi_ice, dz, eff_porosity);
      // compute volumetric liquid water content
      calc_volumetric_h2oliq(eff_porosity, h2osoi_liq, dz, h2osoi_liqvol);
      // calculate root moisture stress
      calc_root_moist_stress(vtype, h2osoi_liqvol, rootfr, rootfr_unf, t_soisno, 
        tc_stress, sucsat, watsat, h2osoi_vol, bsw, smpso, smpsc, eff_porosity, rootr, rresis, btran);

      // Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
      lt = std::min(elai + esai, tlsai_crit);
      egvf = (1.0 - alpha_aero * exp(-lt)) / (1.0 - alpha_aero * exp(-tlsai_crit));
      displa = egvf * displa;
      z0mv   = exp(egvf * std::log(z0mv) + (1.0 - egvf) * std::log(z0mg));
      z0hv   = z0mv;
      z0qv   = z0mv;

      // Net absorbed longwave radiation by canopy and ground
      air =   emv * (1.0 + (1.0 - emv) * (1.0 - emg)) * forc_lwrad;
      bir = - (2.0 - emv * (1.0 - emg)) * emv * sb;
      cir =   emv * emg * sb;

      // Saturated vapor pressure, specific humidity, and their derivatives at the leaf surface
      QSat(t_veg, forc_pbot, el, deldT, qsatl, qsatldT);
      // Determine atmospheric co2 and o2
      co2 = forc_pco2;
      o2  = forc_po2;

      //Initialize flux profile
      nmozsgn = 0;
      taf = (t_grnd + thm) / 2.0;
      qaf = (forc_q + qg) / 2.0;
      ur = std::max(1.0, std::sqrt(forc_u * forc_u + forc_v * forc_v));
      dth = thm - taf;
      dqh = forc_q - qaf;
      delq = qg - qaf;
      dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
      zldis = forc_hgt_u_patch - displa;

      // Check to see if the forcing height is below the canopy height
      assert(zldis >= 0.0);
      // Initialize Monin-Obukhov length and wind speed
      MoninObukIni(ur, tvh, dtvh, zldis, z0mv, um, obu);
    }
  }
}

void StabilityIteration()
// clm uses a decreasing index filter to act as a stop criteria
// we operate on single cells, so we will replace the filter criteria with something else
{
  if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno != 0) {
    bool stop = false;
    while (itlef <= itmax && !stop) {
      // Determine friction velocity, and potential temperature and humidity profiles of the surface boundary layer
      FrictionVelocityWind(forc_hgt_u_patch, displa, um, obu, z0mg, ustar);
      FrictionVelocityTemperature(forc_hgt_t_patch, displa, obu, z0hg, temp1);
      FrictionVelocityHumidity(forc_hgt_q_patch, forc_hgt_t_patch, displa, obu, z0hg, z0qg, temp1, temp2);

      // save leaf temp and leaf temp delta from previous iteration
      tlbef(p) = t_veg(p)
      del2(p) = del(p)
      // Determine aerodynamic resistances
      ram1   = 1.0 / (ustar * ustar / um);
      rah[0] = 1.0 / (temp1 * ustar);
      raw[0] = 1.0 / (temp2 * ustar);
      // Bulk boundary layer resistance of leaves
      uaf = um * std::sqrt(1.0 / (ram1 * um));
      // Use pft parameter for leaf characteristic width dleaf_patch
      dleaf_patch = dleaf[vtype];
      cf  = 0.01 / (std::sqrt(uaf) * std::sqrt(dleaf_patch));
      rb(p)  = 1.0 / (cf * uaf);

      //Parameterization for variation of csoilc with canopy density from X. Zeng, University of Arizona

      w = exp(-(elai + esai));

      // changed by K.Sakaguchi from here
      // transfer coefficient over bare soil is changed to a local variable
      // just for readability of the code (from line 680)
      csoilb = (vkc / (0.13 * pow((z0mg * uaf / 1.5e-5), 0.45)));
      //compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)
      ri = (grav * htop * (taf - t_grnd)) / pow((taf * uaf), 2.0);

      // modify csoilc value (0.004) if the under-canopy is in stable condition
      if ((taf - t_grnd) > 0.0) {
        // decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
        // ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
        ricsoilc = csoilc / (1.00 + ria * std::min(ri, 10.0));
        csoilcn = csoilb * w + ricsoilc * (1.0 - w);
      } else {
        csoilcn = csoilb * w + csoilc * (1.0 - w);
      }

      // Sakaguchi changes for stability formulation ends here
      rah[1] = 1.0 / (csoilcn * uaf);
      raw[1] = rah[1];
      // Stomatal resistances for sunlit and shaded fractions of canopy.
      // Done each iteration to account for differences in eah, tv.
      svpts = el;                      // pa
      eah = forc_pbot * qaf / 0.622;   // pa
      rhaf = eah / svpts;

      if (vtype == nsoybean || vtype == nsoybeanirrig) { btran = std::min(1.0, btran * 1.25); }

      // call photosynthesis
      // Photosynthesis(phase=sun); need to implement

      if (vtype == nsoybean || vtype == nsoybeanirrig) { btran = std::min(1.0, btran * 1.25); }

      // call photosynthesis
      // Photosynthesis(phase=shade); need to implement

      // Sensible heat conductance for air, leaf and ground
      wta   = 1.0 / rah[0];             // air
      wtl   = (elai + esai) / rb;       // leaf
      wtg   = 1.0 / rah[1];             // ground
      wtshi = 1.0 / (wta + wtl + wtg);
      wtl0  = wtl * wtshi;              // leaf
      wtg0  = wtg * wtshi;              // ground
      wta0  = wta * wtshi;              // air
      wtga  = wta0 + wtg0;              // ground + air
      wtal  = wta0 + wtl0;              // air + leaf

      // Fraction of potential evaporation from leaf
       if (fdry > 0.0) {
        rppdry = fdry * rb * (laisun / (rb + rssun) + laisha / (rb + rssha)) / elai;
       } else { rppdry = 0.0; }

       efpot = forc_rho * wtl * (qsatl - qaf);
        if (efpot > 0.0) {
          if (btran > btran0) {
            qflx_tran_veg = efpot * rppdry;
            rpp = rppdry + fwet;
          } else {
            // No transpiration if btran below 1.e-10
            rpp = fwet;
            qflx_tran_veg = 0.0;
          }
          // Check total evapotranspiration from leaves
          rpp = std::min(rpp, (qflx_tran_veg + h2ocan / dtime) / efpot);
        } else {
          // No transpiration if potential evaporation less than zero
          rpp = 1.0;
          qflx_tran_veg = 0.0;
        }

        // Update conductances for changes in rpp
        // Latent heat conductances for ground and leaf.
        // Air has same conductance for both sensible and latent heat.
        wtaq    = frac_veg_nosno / raw[0];                   // air
        wtlq    = frac_veg_nosno * (elai + esai) / rb * rpp; // leaf
        // Litter layer resistance. Added by K.Sakaguchi
        snow_depth_c = z_dl;                                 // critical depth for 100% litter burial by snow (=litter thickness)
        fsno_dl = snow_depth / snow_depth_c;                 // effective snow cover for (dry)plant litter
        elai_dl = lai_dl * (1.0 - std::min(fsno_dl, 1.0));   // exposed (dry)litter area index
        rdl = ( 1.0 - exp(-elai_dl)) / (0.004 * uaf);        // dry litter layer resistance
        // add litter resistance and Lee and Pielke 1992 beta
        if (delq < 0.0) { // dew. Do not apply beta for negative flux (follow old rsoil)
          wtgq = frac_veg_nosno / (raw[1] + rdl);
        } else {
          wtgq = soilbeta * frac_veg_nosno / (raw[1] + rdl);
        }

        wtsqi = 1.0 / (wtaq + wtlq + wtgq);
        wtgq0 = wtgq * wtsqi;           // ground
        wtlq0 = wtlq * wtsqi;           // leaf
        wtaq0 = wtaq * wtsqi;           // air
        wtgaq = wtaq0 + wtgq0;          // air + ground
        wtalq = wtaq0 + wtlq0;          // air + leaf
        dc1   = forc_rho * cpair * wtl;
        dc2   = hvap * forc_rho * wtlq;
        efsh  = dc1 * (wtga * t_veg - wtg0 * t_grnd - wta0 * thm);

        // Evaporation flux from foliage
        erre = 0._r8
        if (efe * efeb < 0.0) {
          efeold = efe;
          efe(p)  = 0.1 * efeold;
          erre = efe - efeold;
         }

        // fractionate ground emitted longwave
        lw_grnd = (frac_sno * pow(t_soisno[nlevsno-snl], 4.0) + 
          (1.0 - frac_sno - frac_h2osfc) * pow(t_soisno[nlevsno], 4.0) + frac_h2osfc * pow(t_h2osfc, 4.0));
        dt_veg = (sabv + air + bir * pow(t_veg, 4.0) + cir * lw_grnd - efsh - efe) / 
          (-4.0 * bir * pow(t_veg, 3.0) + dc1 * wtga + dc2 * wtgaq * qsatldT);
        t_veg = tlbef + dt_veg;
        dels = dt_veg;
        del  = std::abs(dels);
        err = 0.0;
        if (del > delmax) {
          dt_veg = delmax * dels / del;
          t_veg = tlbef + dt_veg;
          err = sabv + air + bir * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg) + cir * lw_grnd - 
            (efsh + dc1 * wtga * dt_veg) - (efe + dc2 * wtgaq * qsatldT * dt_veg);
        }

        // Fluxes from leaves to canopy space
        // "efe" was limited as its sign changes frequently.  This limit may
        // result in an imbalance in "hvap*qflx_evap_veg" and
        // "efe + dc2*wtgaq*qsatdt_veg"
        efpot = forc_rho * wtl * (wtgaq * (qsatl + qsatldT * dt_veg) - wtgq0 * qg - wtaq0 * forc_q);
        qflx_evap_veg = rpp * efpot;
        // Calculation of evaporative potentials (efpot) and interception losses; flux in kg m**-2 s-1.  
        // ecidif holds the excess energy if all intercepted water is evaporated during the timestep.  
        // This energy is later added to the sensible heat flux.
        ecidif = 0.0;
        if (efpot > 0.0 && btran > btran0) {
          qflx_tran_veg = efpot * rppdry;
        } else { qflx_tran_veg = 0.0; }
        ecidif = std::max(0.0, qflx_evap_veg - qflx_tran_veg - h2ocan / dtime);
        qflx_evap_veg = std::min(qflx_evap_veg, qflx_tran_veg + h2ocan / dtime);
        // The energy loss due to above two limits is added to the sensible heat flux.
        eflx_sh_veg = efsh + dc1 * wtga * dt_veg + err + erre + hvap * ecidif;
        //Re-calculate saturated vapor pressure, specific humidity, and their derivatives at the leaf surface
        QSat(t_veg, forc_pbot, el, deldT, qsatl, qsatldT);

        // Update vegetation/ground surface temperature, canopy air
        // temperature, canopy vapor pressure, aerodynamic temperature, and
        // Monin-Obukhov stability parameter for next iteration.
        taf = wtg0 * t_grnd + wta0 * thm + wtl0 * t_veg;
        qaf = wtlq0 * qsatl + wtgq0 * qg + forc_q * wtaq0;
        // Update Monin-Obukhov length and wind speed including the stability effect
        dth = thm - taf;
        dqh = forc_q - qaf;
        delq = wtalq * qg - wtlq0 * qsatl - wtaq0 * forc_q;
        tstar = temp1 * dth;
        qstar = temp2 * dqh;
        thvstar = tstar * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * qstar;
        zeta = zldis * vkc * grav * thvstar / (pow(ustar, 2.0) * thv);
        if (zeta >= 0._r8) {      //stable
          zeta = std::min(2.0, std::max(zeta, 0.01));
          um = std::max(ur, 0.1);
        } else {
          zeta = std::max(-100.0, std::min(zeta, -0.01));
          wc = beta * pow((-grav * ustar * thvstar * zii / thv), 0.333);
          um = std::sqrt(ur * ur + wc * wc);
        }
        obu = zldis / zeta;
        if (obuold * obu < 0.) { nmozsgn += 1; }
        if (nmozsgn >= 4) { obu = zldis / (-0.01); }
        obuold = obu;
        // laminar boundary resistance for h2o over leaf?
        lbl_rsc_h2o = getlblcef(forc_rho, t_veg) * uaf / (pow(uaf, 2.0) + 1.e-10);

        // Test for convergence
        itlef += 1;
        if (itlef > itmin) {
          dele = std::abs(efe - efeb);
          efeb = efe;
          det  = std::max(del, del2);
          if ((det < dtmin) && (dele < dlemin)) { stop = true; }
        }
    } // stability iteration
  } // land type
} //CanopyFluxes::StabilityIteration()



void ComputeFlux()
{

  if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno != 0) {

    // Energy balance check in canopy

    lw_grnd = (frac_sno * pow(t_soisno[nlevsno-snl], 4.0) + (1.0 - frac_sno - frac_h2osfc) * 
      pow(t_soisno[nlevsno], 4.0) + frac_h2osfc * pow(t_h2osfc, 4.0));
    err = sabv + air + bir * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg) + cir * lw_grnd - 
      eflx_sh_veg - hvap * qflx_evap_veg;
    // Fluxes from ground to canopy space
    delt = wtal * t_grnd - wtl0 * t_veg - wta0 * thm;
    taux = -forc_rho * forc_u / ram1;
    tauy = -forc_rho * forc_v / ram1;
    eflx_sh_grnd = cpair * forc_rho * wtg * delt;
    // compute individual sensible heat fluxes
    delt_snow = wtal * t_soisno[nlevsno-snl] - wtl0 * t_veg - wta0 * thm;
    eflx_sh_snow = cpair * forc_rho * wtg * delt_snow;
    delt_soil = wtal * t_soisno[nlevsno] - wtl0 * t_veg - wta0 * thm;
    eflx_sh_soil = cpair * forc_rho * wtg * delt_soil;
    delt_h2osfc  = wtal * t_h2osfc - wtl0 * t_veg - wta0 * thm;
    eflx_sh_h2osfc = cpair * forc_rho * wtg * delt_h2osfc;
    qflx_evap_soi = forc_rho * wtgq * delq;
    // compute individual latent heat fluxes
    delq_snow = wtalq * qg_snow - wtlq0 * qsatl - wtaq0 * forc_q;
    qflx_ev_snow = forc_rho * wtgq * delq_snow;
    delq_soil = wtalq * qg_soil - wtlq0 * qsatl - wtaq0 * forc_q;
    qflx_ev_soil = forc_rho * wtgq * delq_soil;
    delq_h2osfc = wtalq * qg_h2osfc - wtlq0 * qsatl - wtaq0 * forc_q;
    qflx_ev_h2osfc = forc_rho * wtgq * delq_h2osfc;
    // Downward longwave radiation below the canopy
    dlrad = (1.0 - emv) * emg * forc_lwrad + emv * emg * sb * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg);
    // Upward longwave radiation above the canopy
    ulrad = ((1.0 - emg) * (1.0 - emv) * (1.0 - emv) * forc_lwrad + emv * (1.0 + (1.0 - emg) * (1.0 - emv)) * 
      sb * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg) + emg * (1.0 - emv) * sb * lw_grnd);
    //Derivative of soil energy flux with respect to soil temperature
    cgrnds = cgrnds + cpair * forc_rho * wtg * wtal;
    cgrndl = cgrndl + forc_rho * wtgq * wtalq * dqgdT;
    cgrnd  = cgrnds + cgrndl * htvp;
    // Update dew accumulation (kg/m2)
    h2ocan = std::max(0.0, h2ocan + (qflx_tran_veg - qflx_evap_veg) * dtime);

    // Determine total photosynthesis -- need to implement
    // PhotosynthesisTotal(fn, filterp, atm2lnd_vars, cnstate_vars, canopystate_vars, photosyns_vars)

    // evaluate error and write??
    //  if (abs(err(p)) > 0.1_r8)




  }



}

};

