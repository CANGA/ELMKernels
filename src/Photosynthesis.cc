#include <stdexcept>
#include <cmath>

namespace Photosynthesis {

// DESCRIPTION: evaluate the function f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
void quadratic(const double& a, const double& b, const double& c, double& r1, double& r2) {
  double q;
  if (a == 0.0) { throw std::runtime_error("ELM ERROR: quadratic solution a == 0.0");  } //error! end run
  if (b >= 0.0) {
    q = -0.5 * (b + std::sqrt(b * b - 4.0 * a * c));
  } else {
    q = -0.5 * (b - std::sqrt(b * b - 4.0 * a * c));
  }
  r1 = q / a;
  if (q /= 0.0) {
    r2 = c / q;
  } else {
    r2 = 1.0e36;
  }
}

void TimeStepInit()
{
}

// DESCRIPTION: evaluate the function f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
void ci_func(
  const double& ci,  // intracellular leaf CO2 (Pa)
  double& fval, // return function of the value f(ci)
  const int& iv, // canopy index
  const double& gb_mol,  // leaf boundary layer conductance (umol H2O/m**2/s)
  const double& je,  // electron transport rate (umol electrons/m**2/s)
  const double& cair,  // Atmospheric CO2 partial pressure (Pa)
  const double& oair,  // Atmospheric O2 partial pressure (Pa)
  const double& lmr_z,  // canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
  const double& par_z,  // par absorbed per unit lai for canopy layer (w/m**2)
  const double& rh_can,  // canopy air realtive humidity
  double& gs_mol, // leaf stomatal conductance (umol H2O/m**2/s)

  const double *vcmax_z, // [length of nlevcan & nrad)]  maximum rate of carboxylation (umol co2/m**2/s)
  const double& forc_pbot, // atmospheric pressure (Pa)
  const bool& c3flag, // true if C3 and false if C4
  double *ac,        //  [length of nlevcan & nrad)] Rubisco-limited gross photosynthesis (umol CO2/m**2/s)              
  double *aj,        //  [length of nlevcan & nrad)] RuBP-limited gross photosynthesis (umol CO2/m**2/s)                 
  double *ap,        //  [length of nlevcan & nrad)] product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
  double *ag,        //  [length of nlevcan & nrad)] co-limited gross leaf photosynthesis (umol CO2/m**2/s)              
  double *an,        //  [length of nlevcan & nrad)] net leaf photosynthesis (umol CO2/m**2/s)                                               
  const double& cp,        //  CO2 compensation point (Pa)                                           
  const double& kc,        //  Michaelis-Menten constant for CO2 (Pa)                                
  const double& ko,        //  Michaelis-Menten constant for O2 (Pa)                                 
  const double& qe,        //  quantum efficiency, used only for C4 (mol CO2 / mol photons)          
  const double *tpu_z,     //  triose phosphate utilization rate (umol CO2/m**2/s)                 
  const double *kp_z,      //  initial slope of CO2 response curve (C4 plants)                     
  const double& theta_cj,  //  empirical curvature parameter for ac, aj photosynthesis co-limitation 
  const double& bbb,       //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)                 
  const double& mbb)       //  Ball-Berry slope of conductance-photosynthesis relationship           
{
  static const double theta_ip = 0.95;
  double r1,r2; // roots of quadratic equation

  if (c3flag) {
    // C3: Rubisco-limited photosynthesis
    ac[iv] = vcmax_z[iv] * std::max(ci - cp, 0.0) / (ci + kc * (1.0 + oair / ko));
    // C3: RuBP-limited photosynthesis
    aj[iv] = je * std::max(ci - cp, 0.0) / (4.0 * ci + 8.0 * cp);
    // C3: Product-limited photosynthesis 
    ap[iv] = 3.0 * tpu_z[iv];
  } else {
    // C4: Rubisco-limited photosynthesis
    ac[iv] = vcmax_z[iv];
    // C4: RuBP-limited photosynthesis
    aj[iv] = qe * par_z * 4.6;
    // C4: PEP carboxylase-limited (CO2-limited)
    ap[iv] = kp_z[iv] * std::max(ci, 0.0) / forc_pbot;
  }

  // Gross photosynthesis. First co-limit ac and aj. Then co-limit ap
  double aquad = theta_cj; // terms for quadratic equations
  double bquad = -(ac[iv] + aj[iv]);
  double cquad = ac[iv] * aj[iv];
  quadratic(aquad, bquad, cquad, r1, r2);
  double ai = std::min(r1,r2); // intermediate co-limited photosynthesis (umol CO2/m**2/s)

  aquad = theta_ip;
  bquad = -(ai + ap[iv]);
  cquad = ai * ap[iv];
  quadratic(aquad, bquad, cquad, r1, r2);
  ag[iv] = std::min(r1,r2);

  // Net photosynthesis. Exit iteration if an < 0
  an[iv] = ag[iv] - lmr_z;
  if (an[iv] < 0.0) {
     fval = 0.0;
     return;
  }

  // Quadratic gs_mol calculation with an known. Valid for an >= 0.
  // With an <= 0, then gs_mol = bbb
  double cs = cair - 1.4 / gb_mol * an[iv] * forc_pbot; // CO2 partial pressure at leaf surface (Pa)
  cs = std::max(cs, 1.e-6);
  aquad = cs;
  bquad = cs * (gb_mol - bbb) - mbb * an[iv] * forc_pbot;
  cquad = -gb_mol * (cs * bbb + mbb * an[iv] * forc_pbot * rh_can);
  quadratic(aquad, bquad, cquad, r1, r2);
  gs_mol = std::max(r1, r2);

  // Derive new estimate for ci
  fval = ci - cair + an[iv] * forc_pbot * (1.4 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol);
} // ci_func


// DESCRIPTION: Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
// The found root will be updated until its accuracy is tol.
// modified from numerical recipes in F90 by press et al. 1188-1189
void brent (
  double& x, // indepedent variable of the single value function ci_func(x)
  const double& x1, // minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
  const double& x2, 
  const double& f1,
  const double& f2,
  const double& tol, // the error tolerance
  const int& iv, // canopy index
  const double& gb_mol,  // leaf boundary layer conductance (umol H2O/m**2/s)
  const double& je,  // electron transport rate (umol electrons/m**2/s)
  const double& cair,  // Atmospheric CO2 partial pressure (Pa)
  const double& oair,  // Atmospheric O2 partial pressure (Pa)
  const double& lmr_z,  // canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
  const double& par_z,  // par absorbed per unit lai for canopy layer (w/m**2)
  const double& rh_can,  // canopy air realtive humidity
  double& gs_mol, // leaf stomatal conductance (umol H2O/m**2/s)

  const double *vcmax_z, // [length of nlevcan & nrad)]  maximum rate of carboxylation (umol co2/m**2/s)
  const double& forc_pbot, // atmospheric pressure (Pa)
  const bool& c3flag, // true if C3 and false if C4
  double *ac,        //  [length of nlevcan & nrad)] Rubisco-limited gross photosynthesis (umol CO2/m**2/s)              
  double *aj,        //  [length of nlevcan & nrad)] RuBP-limited gross photosynthesis (umol CO2/m**2/s)                 
  double *ap,        //  [length of nlevcan & nrad)] product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
  double *ag,        //  [length of nlevcan & nrad)] co-limited gross leaf photosynthesis (umol CO2/m**2/s)              
  double *an,        //  [length of nlevcan & nrad)] net leaf photosynthesis (umol CO2/m**2/s)                                               
  const double& cp,        //  CO2 compensation point (Pa)                                           
  const double& kc,        //  Michaelis-Menten constant for CO2 (Pa)                                
  const double& ko,        //  Michaelis-Menten constant for O2 (Pa)                                 
  const double& qe,        //  quantum efficiency, used only for C4 (mol CO2 / mol photons)          
  const double *tpu_z,     //  triose phosphate utilization rate (umol CO2/m**2/s)                 
  const double *kp_z,      //  initial slope of CO2 response curve (C4 plants)                     
  const double& theta_cj,  //  empirical curvature parameter for ac, aj photosynthesis co-limitation 
  const double& bbb,       //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)                 
  const double& mbb)       //  Ball-Berry slope of conductance-photosynthesis relationship
{
  static const int ITMAX = 20;  // maximum number of iterations
  static const double EPS = 1.0e-2; // relative error tolerance
  double d,e,p,q,r,s,tol1,xm;
  double a = x1;
  double b = x2;
  double fa = f1;
  double fb = f2;
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {throw std::runtime_error("ELM ERROR: root must be bracketed for brent"); } //error! end run
  double c = b;
  double fc = fb;
  
  int iter = 0;
  while (iter != ITMAX) {
    iter += 1;
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;   // Rename a, b, c and adjust bounding interval d.
      fc = fa;
      d = b-a;
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
    tol1 = 2.0 * EPS * std::abs(b) + 0.5 * tol;  // Convergence check.
    xm = 0.5 * (c-b);
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
      if (p > 0.0) { q *= -1.0; } // Check whether in bounds.
      p = std::abs(p);
      if (2.0 * p < std::min(3.0 * xm * q - std::abs(tol1*q), std::abs(e*q))) {
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
    a = b; //Move last best guess to a.
    fa = fb;
    if (std::abs(d) > tol1) { // Evaluate new trial root.
      b = b + d;
    } else {
      b = b + copysign(tol1,xm);
    }

    ci_func(b, fb, iv, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag,
      ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

    if (fb == 0.0) { break; }
  }

  // if (iter == ITMAX)write(iulog,*) 'brent exceeding maximum iterations', b, fb -- include message???
  x = b;
} // brent

// DESCRIPTION: use a hybrid solver to find the root of equation f(x) = x- h(x), we want to find x, s.t. f(x) = 0.
// the hybrid approach combines the strength of the newton secant approach (find the solution domain)
// and the bisection approach implemented with the Brent's method to guarrantee convergence.
void hybrid(
  double& x0,  // initial guess and final value of the solution
  const int& iv, // canopy index
  const double& gb_mol,  // leaf boundary layer conductance (umol H2O/m**2/s)
  const double& je,  // electron transport rate (umol electrons/m**2/s)
  const double& cair,  // Atmospheric CO2 partial pressure (Pa)
  const double& oair,  // Atmospheric O2 partial pressure (Pa)
  const double& lmr_z,  // canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
  const double& par_z,  // par absorbed per unit lai for canopy layer (w/m**2)
  const double& rh_can,  // canopy air realtive humidity
  double& gs_mol, // leaf stomatal conductance (umol H2O/m**2/s)

  const double *vcmax_z, // [length of nlevcan & nrad)]  maximum rate of carboxylation (umol co2/m**2/s)
  const double& forc_pbot, // atmospheric pressure (Pa)
  const bool& c3flag, // true if C3 and false if C4
  double *ac,        //  [length of nlevcan & nrad)] Rubisco-limited gross photosynthesis (umol CO2/m**2/s)              
  double *aj,        //  [length of nlevcan & nrad)] RuBP-limited gross photosynthesis (umol CO2/m**2/s)                 
  double *ap,        //  [length of nlevcan & nrad)] product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
  double *ag,        //  [length of nlevcan & nrad)] co-limited gross leaf photosynthesis (umol CO2/m**2/s)              
  double *an,        //  [length of nlevcan & nrad)] net leaf photosynthesis (umol CO2/m**2/s)                                               
  const double& cp,        //  CO2 compensation point (Pa)                                           
  const double& kc,        //  Michaelis-Menten constant for CO2 (Pa)                                
  const double& ko,        //  Michaelis-Menten constant for O2 (Pa)                                 
  const double& qe,        //  quantum efficiency, used only for C4 (mol CO2 / mol photons)          
  const double *tpu_z,     //  triose phosphate utilization rate (umol CO2/m**2/s)                 
  const double *kp_z,      //  initial slope of CO2 response curve (C4 plants)                     
  const double& theta_cj,  //  empirical curvature parameter for ac, aj photosynthesis co-limitation 
  const double& bbb,       //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)                 
  const double& mbb)       //  Ball-Berry slope of conductance-photosynthesis relationship
{
  static const double eps = 1.0e-2;      // relative accuracy
  static const double eps1= 1.0e-4;
  static const int itmax = 40;          // maximum number of iterations
  double x1, f0, f1, x, dx, tol, minx, minf;

  ci_func(x0, f0, iv, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag,
      ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

  if (f0 == 0.0) { return; }
  minx = x0;
  minf = f0;
  x1 = x0 * 0.99;

  ci_func(x1, f1, iv, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag,
      ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

  if(f1 == 0.0) {
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
    dx = -f1 * (x1-x0) / (f1-f0);
    x = x1 + dx;
    tol = std::abs(x) * eps;
    if (std::abs(dx) < tol) {
      x0 = x;
      break;
    }
    x0 = x1;
    f0 = f1;
    x1 = x;

    ci_func(x1, f1, iv, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag,
      ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);

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
      brent(x, x0, x1, f0, f1, tol, iv, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, 
        forc_pbot, c3flag, ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);
      x0 = x;
      break;
    }

    if (iter > itmax) { 
      // in case of failing to converge within itmax iterations
      // stop at the minimum function
      // this happens because of some other issues besides the stomatal conductance calculation
      // and it happens usually in very dry places and more likely with c4 plants.
      ci_func(minx, f1, iv, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, vcmax_z, forc_pbot, c3flag,
        ac, aj, ap, ag, an, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, bbb, mbb);
      break;
    }
  } // while loop
} // hybrid






} // namespace Photosynthesis