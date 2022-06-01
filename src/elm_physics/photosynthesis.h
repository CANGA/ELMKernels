/*! \file photosynthesis.h
\brief Internal functions derived from PhotosynthesisMod.F90
*/
#pragma once

#include "elm_constants.h"
#include "pft_data.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "kokkos_includes.hh"

namespace ELM::photosynthesis {

ACCELERATE
double ft(const double& tl, const double& ha);

ACCELERATE
double fth(const double& tl, const double& hd, const double& se, const double& scaleFactor);

ACCELERATE
double fth25(const double& hd, const double& se);

/*! Evaluate the function f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an. (internal) */
ACCELERATE
void quadratic(const double& a, const double& b, const double& c, double& r1, double& r2);

/*! Evaluate the function f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an. (internal) */
ACCELERATE
void ci_func(const double& ci, // intracellular leaf CO2 (Pa)
             double& fval,     // return function of the value f(ci)
             // const int& iv, // canopy index
             const double& gb_mol, // leaf boundary layer conductance (umol H2O/m**2/s)
             const double& je,     // electron transport rate (umol electrons/m**2/s)
             const double& cair,   // Atmospheric CO2 partial pressure (Pa)
             const double& oair,   // Atmospheric O2 partial pressure (Pa)
             const double& lmr_z,  // (single canopy layer) leaf maintenance respiration rate (umol CO2/m**2/s)
             const double& par_z,  // (single canopy layer) par absorbed per unit lai for canopy layer (w/m**2)
             const double& rh_can, // canopy air relative humidity
             double& gs_mol,       // leaf stomatal conductance (umol H2O/m**2/s)

             const double& vcmax_z,   // (single canopy layer)  maximum rate of carboxylation (umol co2/m**2/s)
             const double& forc_pbot, // atmospheric pressure (Pa)
             const bool& c3flag,      // true if C3 and false if C4
             double& ac,              //  (single canopy layer) Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
             double& aj,              //  (single canopy layer) RuBP-limited gross photosynthesis (umol CO2/m**2/s)
             double& ap, //  (single canopy layer) product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol
                         //  CO2/m**2/s)
             double& ag, //  (single canopy layer) co-limited gross leaf photosynthesis (umol CO2/m**2/s)
             double& an, //  (single canopy layer) net leaf photosynthesis (umol CO2/m**2/s)
             const double& cp,       //  CO2 compensation point (Pa)
             const double& kc,       //  Michaelis-Menten constant for CO2 (Pa)
             const double& ko,       //  Michaelis-Menten constant for O2 (Pa)
             const double& qe,       //  quantum efficiency, used only for C4 (mol CO2 / mol photons)
             const double& tpu_z,    //  (single canopy layer) triose phosphate utilization rate (umol CO2/m**2/s)
             const double& kp_z,     //  (single canopy layer) initial slope of CO2 response curve (C4 plants)
             const double& theta_cj, //  empirical curvature parameter for ac, aj photosynthesis co-limitation
             const double& bbb,      //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
             const double& mbb);     //  Ball-Berry slope of conductance-photosynthesis relationship

/*! Use Brent's method to find the root of a single variable function ci_func, which is known to exist
 between x1 and x2. The found root will be updated until its accuracy is tol. modified from numerical recipes in F90
 by press et al. 1188-1189. (internal) */
ACCELERATE
void brent(
    double& x, // indepedent variable of the single value function ci_func(x)
    const double&
        x1, // minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    const double& x2, const double& f1, const double& f2,
    const double& tol, // the error tolerance
    // const int& iv, // canopy index
    const double& gb_mol, // leaf boundary layer conductance (umol H2O/m**2/s)
    const double& je,     // electron transport rate (umol electrons/m**2/s)
    const double& cair,   // Atmospheric CO2 partial pressure (Pa)
    const double& oair,   // Atmospheric O2 partial pressure (Pa)
    const double& lmr_z,  // (single canopy layer) canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    const double& par_z,  // (single canopy layer) par absorbed per unit lai for canopy layer (w/m**2)
    const double& rh_can, // canopy air realtive humidity
    double& gs_mol,       // leaf stomatal conductance (umol H2O/m**2/s)

    const double& vcmax_z,   // (single canopy layer)  maximum rate of carboxylation (umol co2/m**2/s)
    const double& forc_pbot, // atmospheric pressure (Pa)
    const bool& c3flag,      // true if C3 and false if C4
    double& ac,              //  (single canopy layer) Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    double& aj,              //  (single canopy layer) RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    double&
        ap, //  (single canopy layer) product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
    double& ag,             //  (single canopy layer) co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    double& an,             //  (single canopy layer) net leaf photosynthesis (umol CO2/m**2/s)
    const double& cp,       //  CO2 compensation point (Pa)
    const double& kc,       //  Michaelis-Menten constant for CO2 (Pa)
    const double& ko,       //  Michaelis-Menten constant for O2 (Pa)
    const double& qe,       //  quantum efficiency, used only for C4 (mol CO2 / mol photons)
    const double& tpu_z,    //  (single canopy layer) triose phosphate utilization rate (umol CO2/m**2/s)
    const double& kp_z,     //  (single canopy layer) initial slope of CO2 response curve (C4 plants)
    const double& theta_cj, //  empirical curvature parameter for ac, aj photosynthesis co-limitation
    const double& bbb,      //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    const double& mbb);     //  Ball-Berry slope of conductance-photosynthesis relationship

/*! Use a hybrid solver to find the root of equation f(x) = x- h(x), we want to find x, s.t. f(x) = 0.
 the hybrid approach combines the strength of the newton secant approach (find the solution domain)
 and the bisection approach implemented with the Brent's method to guarantee convergence. (internal) */
ACCELERATE
void hybrid(
    double& x0, // initial guess and final value of the solution
    // const int& iv, // canopy index
    const double& gb_mol, // leaf boundary layer conductance (umol H2O/m**2/s)
    const double& je,     // electron transport rate (umol electrons/m**2/s)
    const double& cair,   // Atmospheric CO2 partial pressure (Pa)
    const double& oair,   // Atmospheric O2 partial pressure (Pa)
    const double& lmr_z,  // (single canopy layer) canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    const double& par_z,  // (single canopy layer) par absorbed per unit lai for canopy layer (w/m**2)
    const double& rh_can, // canopy air realtive humidity
    double& gs_mol,       // leaf stomatal conductance (umol H2O/m**2/s)

    const double& vcmax_z,   // (single canopy layer)  maximum rate of carboxylation (umol co2/m**2/s)
    const double& forc_pbot, // atmospheric pressure (Pa)
    const bool& c3flag,      // true if C3 and false if C4
    double& ac,              //  (single canopy layer) Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    double& aj,              //  (single canopy layer) RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    double&
        ap, //  (single canopy layer) product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
    double& ag,             //  (single canopy layer) co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    double& an,             //  (single canopy layer) net leaf photosynthesis (umol CO2/m**2/s)
    const double& cp,       //  CO2 compensation point (Pa)
    const double& kc,       //  Michaelis-Menten constant for CO2 (Pa)
    const double& ko,       //  Michaelis-Menten constant for O2 (Pa)
    const double& qe,       //  quantum efficiency, used only for C4 (mol CO2 / mol photons)
    const double& tpu_z,    //  (single canopy layer) triose phosphate utilization rate (umol CO2/m**2/s)
    const double& kp_z,     //  (single canopy layer) initial slope of CO2 response curve (C4 plants)
    const double& theta_cj, //  empirical curvature parameter for ac, aj photosynthesis co-limitation
    const double& bbb,      //  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    const double& mbb);     //  Ball-Berry slope of conductance-photosynthesis relationship

/*! Compute photosynthesis with iterative solution for vegetation in both sun and shade. (internal) */
template <typename ArrayD1>
ACCELERATE
void photosynthesis(const PFTDataPSN& psn_pft, const int& nrad, const double& forc_pbot, const double& t_veg,
                    const double& t10, const double& esat_tv, const double& eair, const double& oair,
                    const double& cair, const double& rb, const double& btran, const double& dayl_factor,
                    const double& thm, const ArrayD1 tlai_z, const double& vcmaxcint, const ArrayD1 par_z,
                    const ArrayD1 lai_z, ArrayD1 ci_z, double& rs);

/*! Compute photosynthesis totals. (internal)
note: none of these variables do anything - diagnostics maybe??
*/
ACCELERATE
void photosynthesis_total(const double& psnsun, const double& psnsun_wc, const double& psnsun_wj,
                          const double& psnsun_wp, const double& laisun, const double& psnsha, const double& psnsha_wc,
                          const double& psnsha_wj, const double& psnsha_wp, const double& laisha, double& fpsn,
                          double& fpsn_wc, double& fpsn_wj, double& fpsn_wp);

} // namespace ELM::photosynthesis

#include <photosynthesis_impl.hh>
