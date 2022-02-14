// derived from PhotosynthesisMod.F90

#pragma once

namespace ELM::photosynthesis {

template <class ArrayD1>
void photosynthesis(const PSNVegData& psnveg, const int& nrad, const double& forc_pbot, const double& t_veg,
                    const double& t10, const double& esat_tv, const double& eair, const double& oair,
                    const double& cair, const double& rb, const double& btran, const double& dayl_factor,
                    const double& thm, const ArrayD1 tlai_z, const double& vcmaxcint, const ArrayD1 par_z,
                    const ArrayD1 lai_z, double ci_z[nlevcan], double& rs) {

  // photosynthesis and stomatal conductance parameters, from: Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
  const double fnps = 0.15;      // fraction of light absorbed by non-photosynthetic pigments
  const double theta_psii = 0.7; // empirical curvature parameter for electron transport rate
  // C3 or C4 photosynthesis logical variable
  bool c3flag{false};
  if (round(psnveg.c3psn) == 1) {
    c3flag = true;
  } else if (round(psnveg.c3psn) == 0) {
    c3flag = false;
  }

  // Multi-layer parameters scaled by leaf nitrogen profile. Loop through each canopy layer to calculate
  // nitrogen profile using cumulative lai at the midpoint of the layer. Only including code for nu_com == RD && !use_cn

  // Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
  double lnc = 1.0 / (psnveg.slatop * psnveg.leafcn); // leaf N concentration (gN leaf/m^2)
  double act25 = psnveg.act25 * 1000.0 / 60.0;        // (umol/gRubisco/s) Rubisco activity at 25 C
  // vcmax25 at canopy top, as in CN but using lnc at top of the canopy
  double vcmax25top = lnc * psnveg.flnr * psnveg.fnr * act25 *
                      dayl_factor; // canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
  vcmax25top *= psnveg.fnitr;
  // Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
  // canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
  double jmax25top = (2.59 - 0.035 * std::min(std::max((t10 - tfrz), 11.0), 35.0)) * vcmax25top;
  double tpu25top = 0.167 * vcmax25top;  // canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
  double kp25top = 20000.0 * vcmax25top; // canopy top: initial slope of CO2 response curve (C4 plants) at 25C

  // Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used kn = 0.11. Here, derive kn from
  // vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859. Remove daylength factor from vcmax25 so that kn is
  // based on maximum vcmax25 But not used as defined here if using sun/shade big leaf code. Instead, will use canopy
  // integrated scaling factors from SurfaceAlbedo.
  double kn; // leaf nitrogen decay coefficient
  if (dayl_factor == 0.0) {
    kn = 0.0;
  } else {
    kn = exp(0.00963 * vcmax25top / dayl_factor - 2.43);
  }
  // Leaf maintenance respiration in proportion to vcmax25top
  double lmr25top; // canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
  if (c3flag) {
    lmr25top = vcmax25top * 0.015;
  } else {
    lmr25top = vcmax25top * 0.025;
  }

  // Loop through canopy layers (above snow). Respiration needs to be calculated every timestep. Others are calculated
  // only if daytime
  double laican = 0.0;     // canopy sum of lai_z
  double lmr_z[nlevcan];   // canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
  double vcmax_z[nlevcan]; // maximum rate of carboxylation (umol co2/m**2/s)
  double tpu_z[nlevcan];   // patch triose phosphate utilization rate (umol CO2/m**2/s)
  double kp_z[nlevcan];    // patch initial slope of CO2 response curve (C4 plants)
  double jmax_z[nlevcan];  // maximum electron transport rate (umol electrons/m**2/s)
  for (int iv = 0; iv < nrad; iv++) {
    // Cumulative lai at middle of layer
    if (iv == 0) {
      laican = 0.5 * tlai_z[iv];
    } else {
      laican += 0.5 * (tlai_z[iv - 1] + tlai_z[iv]);
    }

    // Scale for leaf nitrogen profile. If multi-layer code, use explicit profile. If sun/shade big leaf code, use
    // canopy integrated factor.
    double nscaler; // leaf nitrogen scaling coefficient
    if (nlevcan == 1) {
      nscaler = vcmaxcint;
    } else if (nlevcan > 1) {
      nscaler = exp(-kn * laican);
    }

    // Maintenance respiration
    double lmr25 = lmr25top * nscaler; // leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    if (c3flag) {
      double lmrc = fth25(psnveg.lmrhd, psnveg.lmrse); // scaling factor for high temperature inhibition (25 C = 1.0)
      lmr_z[iv] = lmr25 * ft(t_veg, psnveg.lmrha) * fth(t_veg, psnveg.lmrhd, psnveg.lmrse, lmrc);
    } else {
      lmr_z[iv] = lmr25 * pow(2.0, ((t_veg - (tfrz + 25.0)) / 10.0));
      lmr_z[iv] /= (1.0 + exp(1.3 * (t_veg - (tfrz + 55.0))));
    }

    if (par_z[iv] <= 0.0) { // night time
      vcmax_z[iv] = 0.0;
      jmax_z[iv] = 0.0;
      tpu_z[iv] = 0.0;
      kp_z[iv] = 0.0;
    } else {                                 // day time
      double vcmax25 = vcmax25top * nscaler; // leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
      double jmax25 = jmax25top * nscaler;   // leaf layer: maximum electron transport rate at 25C (umol electrons/m**
      double tpu25 = tpu25top * nscaler;     // leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
      double kp25 = kp25top * nscaler;       // leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
      // Adjust for temperature
      double vcmaxse = 668.39 - 1.07 * std::min(std::max((t10 - tfrz), 11.0), 35.0); // entropy term for vcmax (J/mol/K)
      double jmaxse = 659.70 - 0.75 * std::min(std::max((t10 - tfrz), 11.0), 35.0);  // entropy term for jmax (J/mol/K)
      double tpuse = vcmaxse;                                                        // entropy term for tpu (J/mol/K)
      double vcmaxc = fth25(psnveg.vcmaxhd, vcmaxse); // scaling factor for high temperature inhibition (25 C = 1.0)
      double jmaxc = fth25(psnveg.jmaxhd, jmaxse);    // scaling factor for high temperature inhibition (25 C = 1.0)
      double tpuc = fth25(psnveg.tpuhd, tpuse);       // scaling factor for high temperature inhibition (25 C = 1.0)
      vcmax_z[iv] = vcmax25 * ft(t_veg, psnveg.vcmaxha) * fth(t_veg, psnveg.vcmaxhd, vcmaxse, vcmaxc);
      jmax_z[iv] = jmax25 * ft(t_veg, psnveg.jmaxha) * fth(t_veg, psnveg.jmaxhd, jmaxse, jmaxc);
      tpu_z[iv] = tpu25 * ft(t_veg, psnveg.tpuha) * fth(t_veg, psnveg.tpuhd, tpuse, tpuc);

      if (!c3flag) {
        vcmax_z[iv] = vcmax25 * pow(2.0, ((t_veg - (tfrz + 25.0)) / 10.0));
        vcmax_z[iv] /= (1.0 + exp(0.2 * ((tfrz + 15.0) - t_veg)));
        vcmax_z[iv] /= (1.0 + exp(0.3 * (t_veg - (tfrz + 40.0))));
      }
      kp_z[iv] = kp25 * pow(2.0, ((t_veg - (tfrz + 25.0)) / 10.0));
    }

    // Adjust for soil water
    vcmax_z[iv] *= btran;
    lmr_z[iv] *= btran;
  } // nrad canopy layers loop

  // Leaf-level photosynthesis and stomatal conductance

  // Leaf boundary layer conductance, umol/m**2/s
  double cf = forc_pbot / (ELM_RGAS * 1.0e-3 * thm) * 1.e06; // s m**2/umol -> s/m
  double gb = 1.0 / rb;                                      // leaf boundary layer conductance (m/s)
  double gb_mol = gb * cf;                                   // leaf boundary layer conductance (umol H2O/m**2/s)

  double psn;               // foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_z[nlevcan];    // canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wc;            // Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wj;            // RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wp;            // product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wc_z[nlevcan]; // Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
  double psn_wj_z[nlevcan]; // RuBP-limited contribution to psn_z (umol CO2/m**2/s)
  double psn_wp_z[nlevcan]; // product-limited contribution to psn_z (umol CO2/m**2/s)
  double rs_z[nlevcan];     // canopy layer: leaf stomatal resistance (s/m)
  double gs_mol[nlevcan];   // leaf stomatal conductance (umol H2O/m**2/s)
  double bbb = std::max(psnveg.bbbopt * btran, 1.0); // Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
  double rsmax0 = 2.0e4;                             // maximum stomatal resistance [s/m]

  double kc25 = (404.9 / 1.e06) * forc_pbot; // Michaelis-Menten constant for CO2 at 25C (Pa)
  double ko25 = (278.4 / 1.e03) * forc_pbot; // Michaelis-Menten constant for O2 at 25C (Pa)
  double cp25 = 0.5 * oair / sco;            // CO2 compensation point at 25C (Pa)
  // account for temperature
  double kc = kc25 * ft(t_veg, psnveg.kcha); // patch Michaelis-Menten constant for CO2 (Pa)
  double ko = ko25 * ft(t_veg, psnveg.koha); // patch Michaelis-Menten constant for O2 (Pa)
  double cp = cp25 * ft(t_veg, psnveg.cpha); // patch CO2 compensation point (Pa)
  // Loop through canopy layers (above snow). Only do calculations if daytime
  for (int iv = 0; iv < nrad; iv++) {
    if (par_z[iv] <= 0.0) { // night time
      // ac[iv] = 0.0;
      // aj[iv] = 0.0;
      // ap[iv] = 0.0;
      // ag[iv] = 0.0;
      // an[iv] = ag[iv] - lmr_z[iv];
      // psn_z[iv] = 0.0;
      // psn_wc_z[iv] = 0.0;
      // psn_wj_z[iv] = 0.0;
      // psn_wp_z[iv] = 0.0;
      ci_z[iv] = 0.0;
      rs_z[iv] = std::min(rsmax0, 1.0 / bbb * cf);
    } else { // day time
      // now the constraint is no longer needed, Jinyun Tang
      double ceair = std::min(eair, esat_tv); // vapor pressure of air, constrained (Pa)
      double rh_can = ceair / esat_tv;        // //  canopy air relative humidity
      // Electron transport rate for C3 plants. Convert par from W/m2 to
      // umol photons/m**2/s using the factor 4.6
      double qabs = 0.5 * (1.0 - fnps) * par_z[iv] * 4.6; // PAR absorbed by PS II (umol photons/m**2/s)
      double aquad = theta_psii;                          // terms for quadratic equations
      double bquad = -(qabs + jmax_z[iv]);                // terms for quadratic equations
      double cquad = qabs * jmax_z[iv];                   // terms for quadratic equations
      double r1, r2;                                      // roots of quadratic equation
      quadratic(aquad, bquad, cquad, r1, r2);
      double je = std::min(r1, r2); // electron transport rate (umol electrons/m**2/s)

      // Iterative loop for ci beginning with initial guess
      if (c3flag) {
        ci_z[iv] = 0.7 * cair;
      } else {
        ci_z[iv] = 0.4 * cair;
      }

      double ciold = ci_z[iv]; // previous value of Ci for convergence check

      // find ci and stomatal conductance
      double ac; // patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
      double aj; // patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
      double ap; // patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
      double ag; // patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
      double an; // patch net leaf photosynthesis (umol CO2/m**2/s)
      hybrid(ciold, gb_mol, je, cair, oair, lmr_z[iv], par_z[iv], rh_can, gs_mol[iv], vcmax_z[iv], forc_pbot, c3flag,
             ac, aj, ap, ag, an, cp, kc, ko, psnveg.qe, tpu_z[iv], kp_z[iv], psnveg.theta_cj, bbb, psnveg.mbbopt);

      // End of ci iteration.  Check for an < 0, in which case gs_mol = bbb
      if (an < 0.0) {
        gs_mol[iv] = bbb;
      }

      //
      double cs = cair - 1.4 / gb_mol * an * forc_pbot; // CO2 partial pressure at leaf surface (Pa)
      cs = std::max(cs, 1.0e-6);
      ci_z[iv] = cair - an * forc_pbot * (1.4 * gs_mol[iv] + 1.6 * gb_mol) / (gb_mol * gs_mol[iv]);
      double gs = gs_mol[iv] / cf; // leaf stomatal conductance (m/s)
      rs_z[iv] = std::min(1.0 / gs, rsmax0);
      // photosynthesis. Save rate-limiting photosynthesis
      psn_z[iv] = ag;
      psn_wc_z[iv] = 0.0;
      psn_wj_z[iv] = 0.0;
      psn_wp_z[iv] = 0.0;

      if (ac <= aj && ac <= ap) {
        psn_wc_z[iv] = psn_z[iv];
      } else if (aj < ac && aj <= ap) {
        psn_wj_z[iv] = psn_z[iv];
      } else if (ap < ac && ap < aj) {
        psn_wp_z[iv] = psn_z[iv];
      }

      // Make sure iterative solution is correct
      if (gs_mol[iv] < 0.0) {
        throw std::runtime_error("ELM ERROR: Negative stomatal conductance");
      }

      // Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b
      // fractional humidity at leaf surface (dimensionless)
      double hs = (gb_mol * ceair + gs_mol[iv] * esat_tv) / ((gb_mol + gs_mol[iv]) * esat_tv);
      double gs_mol_err = psnveg.mbbopt * std::max(an, 0.0) * hs / cs * forc_pbot + bbb; // gs_mol for error check
      if (std::abs(gs_mol[iv] - gs_mol_err) > 1.0e-01) {
        std::cout << "Ball-Berry error check - stomatal conductance error:\n"
                  << gs_mol[iv] << " " << gs_mol_err << "\n";
      }
    } // night/day
  }   // nrad canopy layer loop

  // Canopy photosynthesis and stomatal conductance
  // Sum canopy layer fluxes and then derive effective leaf-level fluxes (per unit leaf area), which are used in other
  // parts of the model. Here, laican sums to either laisun or laisha.

  // double psncan = 0.0; // canopy sum of psn_z
  // double psncan_wc = 0.0; // canopy sum of psn_wc_z
  // double psncan_wj = 0.0; // canopy sum of psn_wj_z
  // double psncan_wp = 0.0; // canopy sum of psn_wp_z
  // double lmrcan = 0.0;
  laican = 0.0;
  double gscan = 0.0; // canopy sum of leaf conductance
  for (int iv = 0; iv < nrad; iv++) {
    // psncan += psn_z[iv] * lai_z[iv];
    // psncan_wc += psn_wc_z[iv] * lai_z[iv];
    // psncan_wj += psn_wj_z[iv] * lai_z[iv];
    // psncan_wp += psn_wp_z[iv] * lai_z[iv];
    // lmrcan += lmr_z[iv] * lai_z[iv];
    gscan += lai_z[iv] / (rb + rs_z[iv]);
    laican += lai_z[iv];
  }
  if (laican > 0.0) { // these variables get used in CN mode, but not in physics (except for rs) -- pass out??
                      // diagnostics maybe??
    // psn = psncan / laican;
    // psn_wc = psncan_wc / laican;
    // psn_wj = psncan_wj / laican;
    // psn_wp = psncan_wp / laican;
    //  lmr = lmrcan / laican;
    rs = laican / gscan - rb;
  } else {
    // psn = 0.0;
    // psn_wc = 0.0;
    // psn_wj = 0.0;
    // psn_wp = 0.0;
    //  lmr = 0.0;
    rs = 0.0;
  }
} // void PhotoSynthesis

} // namespace ELM::photosynthesis
