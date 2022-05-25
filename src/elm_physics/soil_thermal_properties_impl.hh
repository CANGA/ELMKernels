
#pragma once

namespace ELM::soil_thermal {

/*
const int& ltype,
const ArrayD1 h2osoi_liq  [nlevgrnd + nlevsno]
const ArrayD1 h2osoi_ice [nlevgrnd + nlevsno]
const ArrayD1 t_soisno [nlevgrnd + nlevsno]
const ArrayD1 dz  [nlevgrnd + nlevsno]
const ArrayD1 watsat, [nlevgrnd]
tkmg [nlevgrnd]
thk  [nlevgrnd + nlevsno]
tkdry [nlevgrnd]


*/

template <typename ArrayD2>
ACCELERATE
void calc_soil_tk(const int& c,
                  const int& ltype,
                  const ArrayD2 h2osoi_liq,
                  const ArrayD2 h2osoi_ice,
                  const ArrayD2 t_soisno,
                  const ArrayD2 dz,
                  const ArrayD2 watsat,
                  const ArrayD2 tkmg,
                  const ArrayD2 tkdry,
                  ArrayD2 thk)
{
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;

  using detail::TKICE;
  using detail::TKWAT;
  using detail::TKBDRK;

  using ELMconst::DENICE;
  using ELMconst::DENH2O;
  using ELMconst::TFRZ;

  using LND::istwet;
  using LND::istice;
  using LND::istice_mec;

  for (int i = nlevsno; i < nlevgrnd + nlevsno; ++i) {

    if (ltype != istwet && ltype != istice && ltype != istice_mec) {
      double satw = (h2osoi_liq(c, i) / DENH2O + h2osoi_ice(c, i) / DENICE) / (dz(c, i) * watsat(c, i-nlevsno));
      satw = std::min(1.0, satw);
      
      if (satw > 1.0e-6) {
        double dke;
        if (t_soisno(c, i) >= TFRZ) {                 // Unfrozen soil
          dke = std::max(0.0, std::log(satw) + 1.0);
        } else {                                   // Frozen soil
          dke = satw;
        }

        double fl = (h2osoi_liq(c, i) / (DENH2O * dz(c, i))) / (h2osoi_liq(c, i) / (DENH2O * dz(c, i)) +
             h2osoi_ice(c, i) / (DENICE * dz(c, i)));
        double dksat = tkmg(c, i-nlevsno) * pow(TKWAT, fl * watsat(c, i-nlevsno)) * pow(TKICE, (1.0 - fl) * watsat(c, i-nlevsno));
        thk(c, i) = dke * dksat + (1.0 - dke) * tkdry(c, i-nlevsno);
      } else {
         thk(c, i) = tkdry(c, i-nlevsno);
      }

      if (i >= nlevsno + nlevbed) { thk(c, i) = TKBDRK; }

    } else if (ltype == istice || ltype == istice_mec) {

      thk(c, i) = TKWAT;
      if (t_soisno(c, i) < TFRZ) { thk(c, i) = TKICE; }

    } else if (ltype == istwet) {
      
      if (i >= nlevsno + nlevbed) { 
        thk(c, i) = TKBDRK;
      } else {
        thk(c, i) = TKWAT;
        if (t_soisno(c, i) < TFRZ) { thk(c, i) = TKICE; }
      }
    }
  }
}


template <typename ArrayD2>
ACCELERATE
void calc_snow_tk(const int& c,
                  const int& snl,
                  const double& frac_sno,
                  const ArrayD2 h2osoi_liq,
                  const ArrayD2 h2osoi_ice,
                  const ArrayD2 dz,
                  ArrayD2 thk)
{
  static constexpr double TKAIR{0.023};  // thermal conductivity of air   [W/m/K]
  using detail::TKICE;
  using ELMconst::TFRZ;
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;

  // zero out inactive layers
  const int top = nlevsno - snl;
  for (int i = 0; i < top; ++i) {
    thk(c, i) = 0.0;
  }

  // calculate for all active snow layers
  for (int i = top; i < nlevsno; ++i) {
    // partial density of water in the snow pack (ice + liquid) [kg/m3]
    double bw = (h2osoi_ice(c, i) + h2osoi_liq(c, i)) / (frac_sno * dz(c, i));
    thk(c, i) = TKAIR + (7.75e-5 * bw + 1.105e-6 * bw * bw) * (TKICE - TKAIR);
  }
}


// Thermal conductivity at the layer interface
// tk[nlevgrnd+nlevsno]
// tk(i) is the interface between cells i and i+1
// this is different than zi, where zi(i) is between cells i-1 and i
template <typename ArrayD2>
ACCELERATE
void calc_face_tk(const int& c,
                  const int& snl,
                  const ArrayD2 thk,
                  const ArrayD2 z,
                  const ArrayD2 zi,
                  ArrayD2 tk)
{
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;

  // zero out inactive interfaces
  const int top = nlevsno - snl;
  for (int i = 0; i < top; ++i) {
    tk(c, i) = 0.0;
  }

  // active interfaces above bottom interface
  const int bot = nlevgrnd + nlevsno - 1;
  for (int i = top; i < bot - 1; ++i) {
    tk(c, i) = thk(c, i) * thk(c, i+1) * (z(c, i+1) - z(c, i)) /
      (thk(c, i) * (z(c, i+1) - zi(c, i+1)) + thk(c, i+1) * (zi(c, i+1) - z(c, i)));
  }

  // bottom interface
  tk(c, bot) = 0.0;
}


// Soil heat capacity
template <typename ArrayD2>
ACCELERATE
void calc_soil_heat_capacity(const int& c,
                             const int& ltype,
                             const int& snl,
                             const double& h2osno,
                             const ArrayD2 watsat,
                             const ArrayD2 h2osoi_ice,
                             const ArrayD2 h2osoi_liq,
                             const ArrayD2 dz,
                             const ArrayD2 csol,
                             ArrayD2 cv)
{
  using ELMdims::nlevsno;
  using ELMdims::nlevgrnd;

  using ELMconst::CPICE;
  using ELMconst::CPWAT;

  using LND::istwet;
  using LND::istice;
  using LND::istice_mec;

  for (int i = nlevsno; i < nlevgrnd + nlevsno; ++i) {

    if (ltype != istwet && ltype != istice && ltype != istice_mec) {
      cv(c, i) = csol(c, i) * (1.0 - watsat(c, i-nlevsno)) * dz(c, i) + (h2osoi_ice(c, i) * CPICE + h2osoi_liq(c, i) * CPWAT);
    } else if (ltype == istwet) {
      cv(c, i) = (h2osoi_ice(c, i) * CPICE + h2osoi_liq(c, i) * CPWAT);
      if (i >= nlevsno + nlevbed) cv(c, i) = csol(c, i) * dz(c, i);
    } else if (ltype == istice || ltype == istice_mec) {
      cv(c, i) = (h2osoi_ice(c, i) * CPICE + h2osoi_liq(c, i) * CPWAT);
    }

    if (i == nlevsno && snl == 0 && h2osno > 0.0)
    {  cv(c, i) += CPICE * h2osno; }
  }
}


// Snow heat capacity
template <typename ArrayD2>
ACCELERATE
void calc_snow_heat_capacity(const int& c,
                             const int& snl,
                             const double& frac_sno,
                             const ArrayD2 h2osoi_ice,
                             const ArrayD2 h2osoi_liq,
                             ArrayD2 cv)
{
  using ELMdims::nlevsno;
  using ELMconst::CPICE;
  using ELMconst::CPWAT;
  using detail::THIN_SFCLAYER;

   // zero out inactive layers
  const int top = nlevsno - snl;
  for (int i = 0; i < top; ++i) {
    cv(c, i) = 0.0;
  }

  // calculate for all active snow layers
  for (int i = top; i < nlevsno; ++i) {
    if (frac_sno < 0.0) {
      cv(c, i) = std::max(THIN_SFCLAYER, (CPWAT * h2osoi_liq(c, i) + CPICE * h2osoi_ice(c, i)) / frac_sno);
    } else {
      cv(c, i) = THIN_SFCLAYER;
    }
  }
}


// Thermal conductivity of h2osfc
template <typename ArrayD2>
ACCELERATE
double calc_h2osfc_tk(const int& c,
                      const double& h2osfc,
                      const ArrayD2 thk,
                      const ArrayD2 z)
{
  using ELMdims::nlevsno;
  using detail::TKWAT;

  double zh2osfc = 1.0e-3 * (0.5 * h2osfc);  // convert to [m] from [mm]
  return TKWAT * thk(c, nlevsno) * (z(c, nlevsno) + zh2osfc) / (TKWAT * z(c, nlevsno) + thk(c, nlevsno) * zh2osfc);
}


// heat capacity of h2osfc [J/(m2 K)]
ACCELERATE
double calc_h2osfc_heat_capacity(const int& snl, const double& h2osfc, const double& frac_h2osfc)
{
  using ELMconst::CPWAT;
  using detail::THIN_SFCLAYER;
  
  if ((h2osfc > THIN_SFCLAYER) && (frac_h2osfc > THIN_SFCLAYER)) {
    return std::max(THIN_SFCLAYER, CPWAT * h2osfc / frac_h2osfc);
  } else {
    return THIN_SFCLAYER;
  }
}

// height of standing surface water [m]
ACCELERATE
double calc_h2osfc_height(const int& snl, const double& h2osfc, const double& frac_h2osfc)
{
  using detail::THIN_SFCLAYER;

  if ((h2osfc > THIN_SFCLAYER) && (frac_h2osfc > THIN_SFCLAYER)) {
    return std::max(THIN_SFCLAYER, 1.0e-3 * h2osfc / frac_h2osfc);
  } else {
    return THIN_SFCLAYER;
  }
}



} // namespace ELM::soil_thermal
