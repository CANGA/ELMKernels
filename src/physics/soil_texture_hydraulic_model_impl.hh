
#pragma once

namespace ELM {

ACCELERATE
void pedotransfer(const double& pct_sand, const double& pct_clay, double& watsat,
                  double& bsw, double& sucsat, double& xksat)
{
  // compute hydraulic properties based on functions derived from Table 5 in cosby et al, 1984
  // Cosby et al. Table 5
  watsat = 0.489 - 0.00126 * pct_sand;
  bsw = 2.91 + 0.159 * pct_clay;
  sucsat = 10.0 * pow(10.0, (1.88 - 0.0131 * pct_sand));
  xksat = 0.0070556 * pow(10.0, (-0.884 + 0.0153 * pct_sand)); // mm/s, from table 5
}

ACCELERATE
void soil_hydraulic_params(const double& pct_sand, const double& pct_clay,
                           const double& zsoi, const double& om_frac,
                           double& watsat, double& bsw, double& sucsat,
                           double& watdry, double& watopt, double& watfc,
                           double& tkmg, double& tkdry, double& csol)
{
  static constexpr double zsapric{0.5};  // depth (m) that organic matter takes on characteristics of sapric peat
  static constexpr double pcalpha{0.5};  // percolation threshold
  static constexpr double pcbeta{0.139}; // percolation exponent
  static constexpr double om_tkd{0.05};  // thermal conductivity of dry organic soil (Farouki, 1981)
  static constexpr double om_tkm{0.25};  // thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
  static constexpr double csol_bedrock{2.0e6}; // vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
  static constexpr double om_csol{2.5};  // heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)

  double xksat;
  pedotransfer(pct_sand, pct_clay, watsat, bsw, sucsat, xksat);
  const double om_watsat{std::max(0.93 - 0.1 * (zsoi / zsapric), 0.83)};
  const double om_b{std::min(2.7 + 9.3 * (zsoi / zsapric), 12.0)};
  const double om_sucsat{std::min(10.3 - 0.2 * (zsoi / zsapric), 10.1)};
  const double om_hksat{std::max(0.28 - 0.2799 * (zsoi / zsapric), 0.0001)};

  const double bulk_den{(1.0 - watsat) * 2.7e3};
  const double tkm = (1.0 - om_frac) * (8.8 * pct_sand + 2.92 * pct_clay) /
        (pct_sand + pct_clay) + om_tkm * om_frac; // W/(m K)
  watsat = (1.0 - om_frac) * watsat + om_watsat * om_frac;
  bsw = (1.0 - om_frac) * (2.91 + 0.159 * pct_clay) + om_frac * om_b;
  sucsat = (1.0 - om_frac) * sucsat + om_sucsat * om_frac;
  // hksat_min(i) = xksat;

  // perc_frac is zero unless perf_frac greater than percolation threshold
  double perc_frac;
  if (om_frac > pcalpha) {
    double perc_norm = pow((1.0 - pcalpha), -pcbeta);
    perc_frac = perc_norm * pow((om_frac - pcalpha), pcbeta);
  } else {
    perc_frac = 0.0;
  }

  // uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
  double uncon_frac = (1.0 - om_frac) + (1.0 - perc_frac) * om_frac;

  // uncon_hksat is series addition of mineral/organic conductivites
  double uncon_hksat;
  if (om_frac < 1.0) {
    uncon_hksat = uncon_frac / ((1.0 - om_frac) / xksat + ((1.0 - perc_frac) * om_frac) / om_hksat);
  } else {
    uncon_hksat = 0.0;
  }

  double hksat = uncon_frac * uncon_hksat + (perc_frac * om_frac) * om_hksat;

  tkmg = pow(tkm, (1.0 - watsat));

  // this%tksatu_col(c,lev) = this%tkmg_col(c,lev)*0.57_r8**this%watsat_col(c,lev)

  tkdry = ((0.135 * bulk_den + 64.7) / (2.7e3 - 0.947 * bulk_den)) * (1.0 - om_frac) + om_tkd * om_frac;

  csol = ((1.0 - om_frac) * (2.128 * pct_sand + 2.385 * pct_clay) /
         (pct_sand + pct_clay) + om_csol * om_frac) * 1.0e6; // J/(m3 K)

  watdry = watsat * pow((316230.0 / sucsat), (-1.0 / bsw));
  watopt = watsat * pow((158490.0 / sucsat), (-1.0 / bsw));

  // added by K.Sakaguchi for beta from Lee and Pielke, 1992
  // water content at field capacity, defined as hk = 0.1 mm/day
  // used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / secspday (day/sec)
  watfc = watsat * pow((0.1 / (hksat * ELMconst::SECSPDAY())), (1.0 / (2.0 * bsw + 3.0)));

  // this%sucmin_col(c,lev) = min_liquid_pressure

  // this%watmin_col(c,lev) = &
  //      this%watsat_col(c,lev)*(-min_liquid_pressure/this%sucsat_col(c,lev))**(-1._r8/this%bsw_col(c,lev))
}

template <typename ArrayD1>
ACCELERATE
void init_soil_hydraulics(const double& organic_max, const ArrayD1 pct_sand,
                          const ArrayD1 pct_clay, const ArrayD1 organic,
                          const ArrayD1 zsoi, ArrayD1 watsat, ArrayD1 bsw,
                          ArrayD1 sucsat, ArrayD1 watdry, ArrayD1 watopt,
                          ArrayD1 watfc, ArrayD1 tkmg, ArrayD1 tkdry,
                          ArrayD1 csol)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsoi;

  static constexpr double csol_bedrock{2.0e6}; // copied here for now

  double om_frac;
  for (int i = 0; i < nlevsoi(); ++i) {
    om_frac = pow((organic(i) / organic_max), 2.0);
    soil_hydraulic_params(pct_sand(i), pct_clay(i), zsoi(i + nlevsno()), om_frac, watsat(i), bsw(i), sucsat(i),
                          watdry(i), watopt(i), watfc(i),
                          tkmg(i), tkdry(i), csol(i));
  }

  for (int i = nlevsoi(); i < nlevgrnd(); ++i) {
    om_frac = 0.0;
    soil_hydraulic_params(pct_sand(nlevsoi() - 1), pct_clay(nlevsoi() - 1), zsoi(i + nlevsno()), om_frac,
                          watsat(i), bsw(i), sucsat(i), watdry(i), watopt(i), watfc(i), tkmg(i), tkdry(i), csol(i));

    csol(i) = csol_bedrock;
  }
}

} // namespace ELM
