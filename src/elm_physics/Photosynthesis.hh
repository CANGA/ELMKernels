
#pragma once

#include "clm_constants.h"
#include "vegproperties.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ELM {

// DESCRIPTION: evaluate the function f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
void quadratic(const double &a, const double &b, const double &c, double &r1, double &r2) {
  double q;
  if (a == 0.0) {
    throw std::runtime_error("ELM ERROR: quadratic solution a == 0.0");
  } // error! end run
  if (b >= 0.0) {
    q = -0.5 * (b + std::sqrt(b * b - 4.0 * a * c));
  } else {
    q = -0.5 * (b - std::sqrt(b * b - 4.0 * a * c));
  }
  r1 = q / a;
  if (q != 0.0) {
    r2 = c / q;
  } else {
    r2 = 1.0e36;
  }
}

void TimeStepInit() {}

// DESCRIPTION: evaluate the function f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
void ci_func(const double &ci, // intracellular leaf CO2 (Pa)
             double &fval,     // return function of the value f(ci)
             // const int& iv, // canopy index
             const double &gb_mol, // leaf boundary layer conductance (umol H2O/m**2/s)
             const double &je,     // electron transport rate (umol electrons/m**2/s)
             const double &cair,   // Atmospheric CO2 partial pressure (Pa)
             const double &oair,   // Atmospheric O2 partial pressure (Pa)
             const double &lmr_z,  // (single canopy layer) leaf maintenance respiration rate (umol CO2/m**2/s)
             const double &par_z,  // (single canopy layer) par absorbed per unit lai for canopy layer (w/m**2)
             const double &rh_can, // canopy air relative humidity
             double &gs_mol,       // leaf stomatal conductance (umol H2O/m**2/s)

             const double &vcmax_z,   // (single canopy layer)  maximum rate of carboxylation (umol co2/m**2/s)
             const double &forc_pbot, // atmospheric pressure (Pa)
             const bool &c3flag,      // true if C3 and false if C4
             double &ac,              //  (single canopy layer) Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
             double &aj,              //  (single canopy layer) RuBP-limited gross photosynthesis (umol CO2/m**2/s)
             double &ap, //  (single canopy layer) product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol
                         //  CO2/m**2/s)
             double &ag, //  (single canopy layer) co-limited gross leaf photosynthesis (umol CO2/m**2/s)
             double &an, //  (single canopy layer) net leaf photosynthesis (umol CO2/m**2/s)
             const double &cp,       //  CO2 compensation point (Pa)
             const double &kc,       //  Michaelis-Menten constant for CO2 (Pa)
             const double &ko,       //  Michaelis-Menten constant for O2 (Pa)
             const double &qe,       //  quantum efficiency, used only for C4 (mol CO2 / mol photons)
             const double &tpu_z,    //  (single canopy layer) triose phosphate utilization rate (umol CO2/m**2/s)
             const double &kp_z,     //  (single canopy layer) initial slope of CO2 response curve (C4 plants)
             const double &theta_cj, //  empirical curvature parameter for ac, aj photosynthesis co-limitation
             const double &bbb,      //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
             const double &mbb)      //  Ball-Berry slope of conductance-photosynthesis relationship
{
  static const double theta_ip = 0.95;
  double r1, r2; // roots of quadratic equation

  if (c3flag) {
    // C3: Rubisco-limited photosynthesis
    ac = vcmax_z * std::max(ci - cp, 0.0) / (ci + kc * (1.0 + oair / ko));
    // C3: RuBP-limited photosynthesis
    aj = je * std::max(ci - cp, 0.0) / (4.0 * ci + 8.0 * cp);
    // C3: Product-limited photosynthesis
    ap = 3.0 * tpu_z;
  } else {
    // C4: Rubisco-limited photosynthesis
    ac = vcmax_z;
    // C4: RuBP-limited photosynthesis
    aj = qe * par_z * 4.6;
    // C4: PEP carboxylase-limited (CO2-limited)
    ap = kp_z * std::max(ci, 0.0) / forc_pbot;
  }

  // Gross photosynthesis. First co-limit ac and aj. Then co-limit ap
  double aquad = theta_cj; // terms for quadratic equations
  double bquad = -(ac + aj);
  double cquad = ac * aj;
  quadratic(aquad, bquad, cquad, r1, r2);
  double ai = std::min(r1, r2); // intermediate co-limited photosynthesis (umol CO2/m**2/s)

  aquad = theta_ip;
  bquad = -(ai + ap);
  cquad = ai * ap;
  quadratic(aquad, bquad, cquad, r1, r2);
  ag = std::min(r1, r2);

  // Net photosynthesis. Exit iteration if an < 0
  an = ag - lmr_z;
  if (an < 0.0) {
    fval = 0.0;
    return;
  }

  // Quadratic gs_mol calculation with an known. Valid for an >= 0.
  // With an <= 0, then gs_mol = bbb
  double cs = cair - 1.4 / gb_mol * an * forc_pbot; // CO2 partial pressure at leaf surface (Pa)
  cs = std::max(cs, 1.e-6);
  aquad = cs;
  bquad = cs * (gb_mol - bbb) - mbb * an * forc_pbot;
  cquad = -gb_mol * (cs * bbb + mbb * an * forc_pbot * rh_can);
  quadratic(aquad, bquad, cquad, r1, r2);
  gs_mol = std::max(r1, r2);

  // Derive new estimate for ci
  fval = ci - cair + an * forc_pbot * (1.4 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol);
} // ci_func

// DESCRIPTION: Use Brent's method to find the root of a single variable function ci_func, which is known to exist
// between x1 and x2. The found root will be updated until its accuracy is tol. modified from numerical recipes in F90
// by press et al. 1188-1189
void brent(
    double &x, // indepedent variable of the single value function ci_func(x)
    const double
        &x1, // minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    const double &x2, const double &f1, const double &f2,
    const double &tol, // the error tolerance
    // const int& iv, // canopy index
    const double &gb_mol, // leaf boundary layer conductance (umol H2O/m**2/s)
    const double &je,     // electron transport rate (umol electrons/m**2/s)
    const double &cair,   // Atmospheric CO2 partial pressure (Pa)
    const double &oair,   // Atmospheric O2 partial pressure (Pa)
    const double &lmr_z,  // (single canopy layer) canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    const double &par_z,  // (single canopy layer) par absorbed per unit lai for canopy layer (w/m**2)
    const double &rh_can, // canopy air realtive humidity
    double &gs_mol,       // leaf stomatal conductance (umol H2O/m**2/s)

    const double &vcmax_z,   // (single canopy layer)  maximum rate of carboxylation (umol co2/m**2/s)
    const double &forc_pbot, // atmospheric pressure (Pa)
    const bool &c3flag,      // true if C3 and false if C4
    double &ac,              //  (single canopy layer) Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    double &aj,              //  (single canopy layer) RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    double
        &ap, //  (single canopy layer) product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
    double &ag,             //  (single canopy layer) co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    double &an,             //  (single canopy layer) net leaf photosynthesis (umol CO2/m**2/s)
    const double &cp,       //  CO2 compensation point (Pa)
    const double &kc,       //  Michaelis-Menten constant for CO2 (Pa)
    const double &ko,       //  Michaelis-Menten constant for O2 (Pa)
    const double &qe,       //  quantum efficiency, used only for C4 (mol CO2 / mol photons)
    const double &tpu_z,    //  (single canopy layer) triose phosphate utilization rate (umol CO2/m**2/s)
    const double &kp_z,     //  (single canopy layer) initial slope of CO2 response curve (C4 plants)
    const double &theta_cj, //  empirical curvature parameter for ac, aj photosynthesis co-limitation
    const double &bbb,      //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    const double &mbb)      //  Ball-Berry slope of conductance-photosynthesis relationship
{
  static const int ITMAX = 20;      // maximum number of iterations
  static const double EPS = 1.0e-2; // relative error tolerance
  double d, e, p, q, r, s, tol1, xm;
  double a = x1;
  double b = x2;
  double fa = f1;
  double fb = f2;
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    throw std::runtime_error("ELM ERROR: root must be bracketed for brent");
  } // error! end run
  double c = b;
  double fc = fb;

  int iter = 0;
  while (iter != ITMAX) {
    iter += 1;
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a; // Rename a, b, c and adjust bounding interval d.
      fc = fa;
      d = b - a;
      e = d;
    }
    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * EPS * std::abs(b) + 0.5 * tol; // Convergence check.
    xm = 0.5 * (c - b);
    if (std::abs(xm) <= tol1 || fb == 0.0) {
      x = b;
      return;
    }
    if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
      s = fb / fa; // Attempt inverse quadratic interpolation.
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) {
        q *= -1.0;
      } // Check whether in bounds.
      p = std::abs(p);
      if (2.0 * p < std::min(3.0 * xm * q - std::abs(tol1 * q), std::abs(e * q))) {
        e = d; // Accept interpolation.
        d = p / q;
      } else {
        d = xm; // Interpolation failed, use bisection.
        e = d;
      }
    } else { // Bounds decreasing too slowly, use bisection.
      d = xm;
      e = d;
    }
    a = b; // Move last best guess to a.
    fa = fb;
    if (std::abs(d) > tol1) { // Evaluate new trial root.
      b = b + d;
    } else {
      b = b + copysign(tol1, xm);
    }

    ci_func(b, fb, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag, ac, aj, ap, ag, an,
            cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

    if (fb == 0.0) {
      break;
    }
  }

  // if (iter == ITMAX)write(iulog,*) 'brent exceeding maximum iterations', b, fb -- include message???
  x = b;
} // brent

// DESCRIPTION: use a hybrid solver to find the root of equation f(x) = x- h(x), we want to find x, s.t. f(x) = 0.
// the hybrid approach combines the strength of the newton secant approach (find the solution domain)
// and the bisection approach implemented with the Brent's method to guarrantee convergence.
void hybrid(
    double &x0, // initial guess and final value of the solution
    // const int& iv, // canopy index
    const double &gb_mol, // leaf boundary layer conductance (umol H2O/m**2/s)
    const double &je,     // electron transport rate (umol electrons/m**2/s)
    const double &cair,   // Atmospheric CO2 partial pressure (Pa)
    const double &oair,   // Atmospheric O2 partial pressure (Pa)
    const double &lmr_z,  // (single canopy layer) canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    const double &par_z,  // (single canopy layer) par absorbed per unit lai for canopy layer (w/m**2)
    const double &rh_can, // canopy air realtive humidity
    double &gs_mol,       // leaf stomatal conductance (umol H2O/m**2/s)

    const double &vcmax_z,   // (single canopy layer)  maximum rate of carboxylation (umol co2/m**2/s)
    const double &forc_pbot, // atmospheric pressure (Pa)
    const bool &c3flag,      // true if C3 and false if C4
    double &ac,              //  (single canopy layer) Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    double &aj,              //  (single canopy layer) RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    double
        &ap, //  (single canopy layer) product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
    double &ag,             //  (single canopy layer) co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    double &an,             //  (single canopy layer) net leaf photosynthesis (umol CO2/m**2/s)
    const double &cp,       //  CO2 compensation point (Pa)
    const double &kc,       //  Michaelis-Menten constant for CO2 (Pa)
    const double &ko,       //  Michaelis-Menten constant for O2 (Pa)
    const double &qe,       //  quantum efficiency, used only for C4 (mol CO2 / mol photons)
    const double &tpu_z,    //  (single canopy layer) triose phosphate utilization rate (umol CO2/m**2/s)
    const double &kp_z,     //  (single canopy layer) initial slope of CO2 response curve (C4 plants)
    const double &theta_cj, //  empirical curvature parameter for ac, aj photosynthesis co-limitation
    const double &bbb,      //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    const double &mbb)      //  Ball-Berry slope of conductance-photosynthesis relationship
{
  static const double eps = 1.0e-2; // relative accuracy
  static const double eps1 = 1.0e-4;
  static const int itmax = 40; // maximum number of iterations
  double x1, f0, f1, x, dx, tol, minx, minf;

  ci_func(x0, f0, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag, ac, aj, ap, ag, an,
          cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

  if (f0 == 0.0) {
    return;
  }
  minx = x0;
  minf = f0;
  x1 = x0 * 0.99;

  ci_func(x1, f1, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag, ac, aj, ap, ag, an,
          cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

  if (f1 == 0.0) {
    x0 = x1;
    return;
  }
  if (f1 < minf) {
    minx = x1;
    minf = f1;
  }

  // first use the secant approach, then use the brent approach as a backup
  int iter = 0;
  while (true) { // infinite loop until break
    iter += 1;
    dx = -f1 * (x1 - x0) / (f1 - f0);
    x = x1 + dx;
    tol = std::abs(x) * eps;
    if (std::abs(dx) < tol) {
      x0 = x;
      break;
    }
    x0 = x1;
    f0 = f1;
    x1 = x;

    ci_func(x1, f1, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag, ac, aj, ap, ag,
            an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

    if (f1 < minf) {
      minx = x1;
      minf = f1;
    }
    if (std::abs(f1) <= eps1) {
      x0 = x1;
      break;
    }

    // if a root zone is found, use the brent method for a robust backup strategy
    if (f1 * f0 < 0.0) {
      brent(x, x0, x1, f0, f1, tol, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag,
            ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);
      x0 = x;
      break;
    }

    if (iter > itmax) {
      // in case of failing to converge within itmax iterations
      // stop at the minimum function
      // this happens because of some other issues besides the stomatal conductance calculation
      // and it happens usually in very dry places and more likely with c4 plants.
      ci_func(minx, f1, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag, ac, aj, ap,
              ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);
      break;
    }
  } // while loop
} // hybrid

double ft(const double &tl, const double &ha) {
  double ans = exp(ha / (ELM_RGAS * 1.0e-3 * (tfrz + 25.0)) * (1.0 - (tfrz + 25.0) / tl));
  return ans;
}

double fth(const double &tl, const double &hd, const double &se, const double &scaleFactor) {
  double ans = scaleFactor / (1.0 + exp((-hd + se * tl) / (ELM_RGAS * 1.0e-3 * tl)));
  return ans;
}

double fth25(const double &hd, const double &se) {
  double ans = 1.0 + exp((-hd + se * (tfrz + 25.0)) / (ELM_RGAS * 1.0e-3 * (tfrz + 25.0)));
  return ans;
}

/*

t10           => temperature_vars%t_a10_patch             , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2
m temperature (K)


*/
template <class ArrayD1>
void Photosynthesis(const VegProperties &veg, const int &vtype, const int &nrad, const double &forc_pbot,
                    const double &t_veg, const double &t10, const double &esat_tv, const double &eair,
                    const double &oair, const double &cair, const double &rb, const double &btran,
                    const double &dayl_factor, const double &thm, const ArrayD1 tlai_z, const double &vcmaxcint,
                    const ArrayD1 par_z, const ArrayD1 lai_z, double &rs) {
  // vars from PhotosynthesisType - not used anywhere else, so keep local
  bool c3flag;          // true if C3 and false if C4
  double ci_z[nlevcan]; // intracellular leaf CO2 (Pa)
  double rs_z[nlevcan]; // canopy layer: leaf stomatal resistance (s/m)
  // double lmr;              // leaf maintenance respiration rate (umol CO2/m**2/s)
  double lmr_z[nlevcan];   // canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
  double psn;              // foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_z[nlevcan];   // canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wc;           // Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wj;           // RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double psn_wp;           // product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
  double ac[nlevcan];      // patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
  double aj[nlevcan];      // patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
  double ap[nlevcan];      // patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
  double ag[nlevcan];      // patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
  double an[nlevcan];      // patch net leaf photosynthesis (umol CO2/m**2/s)
  double vcmax_z[nlevcan]; // maximum rate of carboxylation (umol co2/m**2/s)
  double gb_mol;           // leaf boundary layer conductance (umol H2O/m**2/s)
  double gs_mol[nlevcan];  // leaf stomatal conductance (umol H2O/m**2/s)
  double cp;               // patch CO2 compensation point (Pa)
  double kc;               // patch Michaelis-Menten constant for CO2 (Pa)
  double ko;               // patch Michaelis-Menten constant for O2 (Pa)
  double qe;               // patch quantum efficiency, used only for C4 (mol CO2 / mol photons)
  double tpu_z[nlevcan];   // patch triose phosphate utilization rate (umol CO2/m**2/s)
  double kp_z[nlevcan];    // patch initial slope of CO2 response curve (C4 plants)
  double theta_cj;         // patch empirical curvature parameter for ac, aj photosynthesis co-limitation
  double bbb;              // patch Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
  double mbb;              // patch Ball-Berry slope of conductance-photosynthesis relationship
  double rh_leaf;          // fractional humidity at leaf surface (dimensionless)

  // Leaf photosynthesis parameters
  double jmax_z[nlevcan]; // maximum electron transport rate (umol electrons/m**2/s)
  double lnc;             // leaf N concentration (gN leaf/m^2)
  double bbbopt;          // Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
  double mbbopt;          // Ball-Berry slope of conductance-photosynthesis relationship, unstressed
  double kn;              // leaf nitrogen decay coefficient
  double vcmax25top;      // canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
  double jmax25top;       // canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
  double tpu25top;        // canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
  double lmr25top;        // canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
  double kp25top;         // canopy top: initial slope of CO2 response curve (C4 plants) at 25C
  double vcmax25;         // leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
  double jmax25;          // leaf layer: maximum electron transport rate at 25C (umol electrons/m**
  double tpu25;           // leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
  double lmr25;           // leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
  double kp25;            // leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
  double kc25;            // Michaelis-Menten constant for CO2 at 25C (Pa)
  double ko25;            // Michaelis-Menten constant for O2 at 25C (Pa)
  double cp25;            // CO2 compensation point at 25C (Pa)
  double vcmaxha;         // activation energy for vcmax (J/mol)
  double jmaxha;          // activation energy for jmax (J/mol)
  double tpuha;           // activation energy for tpu (J/mol)
  double lmrha;           // activation energy for lmr (J/mol)
  double kcha;            // activation energy for kc (J/mol)
  double koha;            // activation energy for ko (J/mol)
  double cpha;            // activation energy for cp (J/mol)
  double vcmaxhd;         // deactivation energy for vcmax (J/mol)
  double jmaxhd;          // deactivation energy for jmax (J/mol)
  double tpuhd;           // deactivation energy for tpu (J/mol)
  double lmrhd;           // deactivation energy for lmr (J/mol)
  double vcmaxse;         // entropy term for vcmax (J/mol/K)
  double jmaxse;          // entropy term for jmax (J/mol/K)
  double tpuse;           // entropy term for tpu (J/mol/K)
  double lmrse;           // entropy term for lmr (J/mol/K)
  double vcmaxc;          // scaling factor for high temperature inhibition (25 C = 1.0)
  double jmaxc;           // scaling factor for high temperature inhibition (25 C = 1.0)
  double tpuc;            // scaling factor for high temperature inhibition (25 C = 1.0)
  double lmrc;            // scaling factor for high temperature inhibition (25 C = 1.0)

  // Other
  double cf;                  //  s m**2/umol -> s/m
  double rsmax0;              //  maximum stomatal resistance [s/m]
  double gb;                  //  leaf boundary layer conductance (m/s)
  double cs;                  //  CO2 partial pressure at leaf surface (Pa)
  double gs;                  //  leaf stomatal conductance (m/s)
  double hs;                  //  fractional humidity at leaf surface (dimensionless)
  double sco;                 //  relative specificity of rubisco
  double ciold;               //  previous value of Ci for convergence check
  double gs_mol_err;          //  gs_mol for error check
  double je;                  //  electron transport rate (umol electrons/m**2/s)
  double qabs;                //  PAR absorbed by PS II (umol photons/m**2/s)
  double aquad, bquad, cquad; //  terms for quadratic equations
  double r1, r2;              //  roots of quadratic equation
  double ceair;               //  vapor pressure of air, constrained (Pa)
  double fnr;                 //  (gRubisco/gN in Rubisco)
  double act25;               //  (umol/mgRubisco/min) Rubisco activity at 25 C
  double nscaler;             //  leaf nitrogen scaling coefficient
  double psn_wc_z[nlevcan];   // Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
  double psn_wj_z[nlevcan];   // RuBP-limited contribution to psn_z (umol CO2/m**2/s)
  double psn_wp_z[nlevcan];   // product-limited contribution to psn_z (umol CO2/m**2/s)
  double psncan;              //  canopy sum of psn_z
  double psncan_wc;           //  canopy sum of psn_wc_z
  double psncan_wj;           //  canopy sum of psn_wj_z
  double psncan_wp;           //  canopy sum of psn_wp_z
  double lmrcan;              //  canopy sum of lmr_z
  double gscan;               //  canopy sum of leaf conductance
  double laican;              //  canopy sum of lai_z
  double rh_can;              //  canopy air relative humidity

  // Photosynthesis and stomatal conductance parameters, from: Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
  const double fnps = 0.15;      // fraction of light absorbed by non-photosynthetic pigments
  const double theta_psii = 0.7; // empirical curvature parameter for electron transport rate

  // vcmax25 parameters, from CN
  fnr = veg.fnr[vtype];
  act25 = veg.act25[vtype]; // umol/mgRubisco/min
  // Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
  act25 *= 1000.0 / 60.0;

  // Activation energy, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  // Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
  // except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282
  kcha = veg.kcha[vtype];
  koha = veg.koha[vtype];
  cpha = veg.cpha[vtype];
  vcmaxha = veg.vcmaxha[vtype];
  jmaxha = veg.jmaxha[vtype];
  tpuha = veg.tpuha[vtype];
  lmrha = veg.lmrha[vtype];

  // High temperature deactivation, from:
  // Leuning (2002) Plant, Cell and Environment 25:1205-1210
  // The factor "c" scales the deactivation to a value of 1.0 at 25C

  vcmaxhd = veg.vcmaxhd[vtype];
  jmaxhd = veg.jmaxhd[vtype];
  tpuhd = veg.tpuhd[vtype];
  lmrhd = veg.lmrhd[vtype];
  lmrse = veg.lmrse[vtype];
  lmrc = fth25(lmrhd, lmrse);

  // C3 or C4 photosynthesis logical variable
  if (round(veg.c3psn[vtype]) == 1) {
    c3flag = true;
  } else if (round(veg.c3psn[vtype]) == 0) {
    c3flag = false;
  }

  // C3 and C4 dependent parameters
  if (c3flag) {
    qe = veg.qe[vtype];
    theta_cj = veg.theta_cj[vtype];
    bbbopt = veg.bbbopt[vtype];
    mbbopt = veg.mbbopt[vtype];
  } else {
    qe = veg.qe[vtype];
    theta_cj = veg.theta_cj[vtype];
    bbbopt = veg.bbbopt[vtype];
    mbbopt = veg.mbbopt[vtype];
  }

  // Soil water stress applied to Ball-Berry parameters
  bbb = std::max(bbbopt * btran, 1.0);
  mbb = mbbopt;

  // Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate cp to account for
  // variation in O2 using cp = 0.5 O2 / sco
  kc25 = (404.9 / 1.e06) * forc_pbot;
  ko25 = (278.4 / 1.e03) * forc_pbot;
  sco = 0.5 * 0.209 / (42.75 / 1.e06);
  cp25 = 0.5 * oair / sco;
  kc = kc25 * ft(t_veg, kcha);
  ko = ko25 * ft(t_veg, koha);
  cp = cp25 * ft(t_veg, cpha);

  // Multi-layer parameters scaled by leaf nitrogen profile. Loop through each canopy layer to calculate
  // nitrogen profile using cumulative lai at the midpoint of the layer. Only including code for nu_com == RD && !use_cn

  // Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
  lnc = 1.0 / (veg.slatop[vtype] * veg.leafcn[vtype]);
  // vcmax25 at canopy top, as in CN but using lnc at top of the canopy
  vcmax25top = lnc * veg.flnr[vtype] * fnr * act25 * dayl_factor;
  vcmax25top *= veg.fnitr[vtype];
  // Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
  jmax25top = (2.59 - 0.035 * std::min(std::max((t10 - tfrz), 11.0), 35.0)) * vcmax25top;
  tpu25top = 0.167 * vcmax25top;
  kp25top = 20000.0 * vcmax25top;

  // Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used kn = 0.11. Here, derive kn from
  // vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859. Remove daylength factor from vcmax25 so that kn is
  // based on maximum vcmax25 But not used as defined here if using sun/shade big leaf code. Instead, will use canopy
  // integrated scaling factors from SurfaceAlbedo.
  if (dayl_factor == 0.0) {
    kn = 0.0;
  } else {
    kn = exp(0.00963 * vcmax25top / dayl_factor - 2.43);
  }
  // Leaf maintenance respiration in proportion to vcmax25top
  if (c3flag) {
    lmr25top = vcmax25top * 0.015;
  } else {
    lmr25top = vcmax25top * 0.025;
  }

  // Loop through canopy layers (above snow). Respiration needs to be calculated every timestep. Others are calculated
  // only if daytime
  laican = 0.0;
  for (int iv = 0; iv < nrad; iv++) {
    // Cumulative lai at middle of layer
    if (iv == 0) {
      laican = 0.5 * tlai_z[iv];
    } else {
      laican = laican + 0.5 * (tlai_z[iv - 1] + tlai_z[iv]);
    }

    // Scale for leaf nitrogen profile. If multi-layer code, use explicit profile. If sun/shade big leaf code, use
    // canopy integrated factor.
    if (nlevcan == 1) {
      nscaler = vcmaxcint;
    } else if (nlevcan > 1) {
      nscaler = exp(-kn * laican);
    }

    // Maintenance respiration
    lmr25 = lmr25top * nscaler;
    if (c3flag) {
      lmr_z[iv] = lmr25 * ft(t_veg, lmrha) * fth(t_veg, lmrhd, lmrse, lmrc);
    } else {
      lmr_z[iv] = lmr25 * pow(2.0, ((t_veg - (tfrz + 25.0)) / 10.0));
      lmr_z[iv] /= (1.0 + exp(1.3 * (t_veg - (tfrz + 55.0))));
    }

    if (par_z[iv] <= 0.0) { // night time
      vcmax_z[iv] = 0.0;
      jmax_z[iv] = 0.0;
      tpu_z[iv] = 0.0;
      kp_z[iv] = 0.0;
    } else { // day time
      vcmax25 = vcmax25top * nscaler;
      jmax25 = jmax25top * nscaler;
      tpu25 = tpu25top * nscaler;
      kp25 = kp25top * nscaler;
      // Adjust for temperature
      vcmaxse = 668.39 - 1.07 * std::min(std::max((t10 - tfrz), 11.0), 35.0);
      jmaxse = 659.70 - 0.75 * std::min(std::max((t10 - tfrz), 11.0), 35.0);
      tpuse = vcmaxse;
      vcmaxc = fth25(vcmaxhd, vcmaxse);
      jmaxc = fth25(jmaxhd, jmaxse);
      tpuc = fth25(tpuhd, tpuse);
      vcmax_z[iv] = vcmax25 * ft(t_veg, vcmaxha) * fth(t_veg, vcmaxhd, vcmaxse, vcmaxc);
      jmax_z[iv] = jmax25 * ft(t_veg, jmaxha) * fth(t_veg, jmaxhd, jmaxse, jmaxc);
      tpu_z[iv] = tpu25 * ft(t_veg, tpuha) * fth(t_veg, tpuhd, tpuse, tpuc);

      if (!c3flag) {
        vcmax_z[iv] = vcmax25 * pow(2.0, ((t_veg - (tfrz + 25.0)) / 10.0));
        vcmax_z[iv] = vcmax_z[iv] / (1.0 + exp(0.2 * ((tfrz + 15.0) - t_veg)));
        vcmax_z[iv] = vcmax_z[iv] / (1.0 + exp(0.3 * (t_veg - (tfrz + 40.0))));
      }
      kp_z[iv] = kp25 * pow(2.0, ((t_veg - (tfrz + 25.0)) / 10.0));
    }

    // Adjust for soil water
    vcmax_z[iv] = vcmax_z[iv] * btran;
    lmr_z[iv] = lmr_z[iv] * btran;
  } // nrad canopy layers loop

  // Leaf-level photosynthesis and stomatal conductance
  rsmax0 = 2.0e4;
  // Leaf boundary layer conductance, umol/m**2/s
  cf = forc_pbot / (ELM_RGAS * 1.0e-3 * thm) * 1.e06;
  gb = 1.0 / rb;
  gb_mol = gb * cf;

  // Loop through canopy layers (above snow). Only do calculations if daytime
  for (int iv = 0; iv < nrad; iv++) {
    if (par_z[iv] <= 0.0) { // night time
      ac[iv] = 0.0;
      aj[iv] = 0.0;
      ap[iv] = 0.0;
      ag[iv] = 0.0;
      an[iv] = ag[iv] - lmr_z[iv];
      psn_z[iv] = 0.0;
      psn_wc_z[iv] = 0.0;
      psn_wj_z[iv] = 0.0;
      psn_wp_z[iv] = 0.0;
      rs_z[iv] = std::min(rsmax0, 1.0 / bbb * cf);
      ci_z[iv] = 0.0;
      rh_leaf = 0.0;
    } else { // day time
      // now the constraint is no longer needed, Jinyun Tang
      ceair = std::min(eair, esat_tv);
      rh_can = ceair / esat_tv;
      // Electron transport rate for C3 plants. Convert par from W/m2 to
      // umol photons/m**2/s using the factor 4.6
      qabs = 0.5 * (1.0 - fnps) * par_z[iv] * 4.6;
      aquad = theta_psii;
      bquad = -(qabs + jmax_z[iv]);
      cquad = qabs * jmax_z[iv];
      quadratic(aquad, bquad, cquad, r1, r2);
      je = std::min(r1, r2);

      // Iterative loop for ci beginning with initial guess
      if (c3flag) {
        ci_z[iv] = 0.7 * cair;
      } else {
        ci_z[iv] = 0.4 * cair;
      }

      ciold = ci_z[iv]; // Save old ci

      // find ci and stomatal conductance
      hybrid(ciold, gb_mol, je, cair, oair, lmr_z[iv], par_z[iv], rh_can, gs_mol[iv], vcmax_z[iv], forc_pbot, c3flag,
             ac[iv], aj[iv], ap[iv], ag[iv], an[iv], cp, kc, ko, qe, tpu_z[iv], kp_z[iv], theta_cj, bbb, mbb);

      // End of ci iteration.  Check for an < 0, in which case gs_mol = bbb
      if (an[iv] < 0.0) {
        gs_mol[iv] = bbb;
      }
      // Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)
      cs = cair - 1.4 / gb_mol * an[iv] * forc_pbot;
      cs = std::max(cs, 1.0e-6);
      ci_z[iv] = cair - an[iv] * forc_pbot * (1.4 * gs_mol[iv] + 1.6 * gb_mol) / (gb_mol * gs_mol[iv]);
      // Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)
      gs = gs_mol[iv] / cf;
      rs_z[iv] = std::min(1.0 / gs, rsmax0);
      // Photosynthesis. Save rate-limiting photosynthesis
      psn_z[iv] = ag[iv];
      psn_wc_z[iv] = 0.0;
      psn_wj_z[iv] = 0.0;
      psn_wp_z[iv] = 0.0;

      if (ac[iv] <= aj[iv] && ac[iv] <= ap[iv]) {
        psn_wc_z[iv] = psn_z[iv];
      } else if (aj[iv] < ac[iv] && aj[iv] <= ap[iv]) {
        psn_wj_z[iv] = psn_z[iv];
      } else if (ap[iv] < ac[iv] && ap[iv] < aj[iv]) {
        psn_wp_z[iv] = psn_z[iv];
      }

      // Make sure iterative solution is correct
      if (gs_mol[iv] < 0.0) {
        throw std::runtime_error("ELM ERROR: Negative stomatal conductance");
      } // error! end run -- message??

      // Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b
      hs = (gb_mol * ceair + gs_mol[iv] * esat_tv) / ((gb_mol + gs_mol[iv]) * esat_tv);
      rh_leaf = hs;
      gs_mol_err = mbb * std::max(an[iv], 0.0) * hs / cs * forc_pbot + bbb;
      if (std::abs(gs_mol[iv] - gs_mol_err) > 1.0e-01) {
      } // write output message??
    }   // night/day
  }     // nrad canopy layer loop

  // Canopy photosynthesis and stomatal conductance
  // Sum canopy layer fluxes and then derive effective leaf-level fluxes (per unit leaf area), which are used in other
  // parts of the model. Here, laican sums to either laisun or laisha.
  psncan = 0.0;
  psncan_wc = 0.0;
  psncan_wj = 0.0;
  psncan_wp = 0.0;
  lmrcan = 0.0;
  gscan = 0.0;
  laican = 0.0;
  for (int iv = 0; iv < nrad; iv++) {
    psncan += psn_z[iv] * lai_z[iv];
    psncan_wc += psn_wc_z[iv] * lai_z[iv];
    psncan_wj += psn_wj_z[iv] * lai_z[iv];
    psncan_wp += psn_wp_z[iv] * lai_z[iv];
    lmrcan += lmr_z[iv] * lai_z[iv];
    gscan += lai_z[iv] / (rb + rs_z[iv]);
    laican += +lai_z[iv];
  }
  if (laican > 0.0) { // these variables get used in CN mode, but not in physics (except for rs) -- pass out??
                      // diagnostics maybe??
    psn = psncan / laican;
    psn_wc = psncan_wc / laican;
    psn_wj = psncan_wj / laican;
    psn_wp = psncan_wp / laican;
    // lmr = lmrcan / laican;
    rs = laican / gscan - rb;
  } else {
    psn = 0.0;
    psn_wc = 0.0;
    psn_wj = 0.0;
    psn_wp = 0.0;
    // lmr = 0.0;
    rs = 0.0;
  }
} // void PhotoSynthesis

/* none of these variables do anything - diagnostics maybe??



*/
void PhotosynthesisTotal(const double &psnsun, const double &psnsun_wc, const double &psnsun_wj,
                         const double &psnsun_wp, const double &laisun, const double &psnsha, const double &psnsha_wc,
                         const double &psnsha_wj, const double &psnsha_wp, const double &laisha, double &fpsn,
                         double &fpsn_wc, double &fpsn_wj, double &fpsn_wp) {
  fpsn = psnsun * laisun + psnsha * laisha;
  fpsn_wc = psnsun_wc * laisun + psnsha_wc * laisha;
  fpsn_wj = psnsun_wj * laisun + psnsha_wj * laisha;
  fpsn_wp = psnsun_wp * laisun + psnsha_wp * laisha;
}

} // namespace ELM
