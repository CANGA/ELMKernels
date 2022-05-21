
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

template <typename ArrayD1>
ACCELERATE
void calc_soil_tk(const int& ltype,
                  const ArrayD1 h2osoi_liq,
                  const ArrayD1 h2osoi_ice,
                  const ArrayD1 t_soisno,
                  const ArrayD1 dz,
                  const ArrayD1 watsat,
                  const ArrayD1 tkmg,
                  const ArrayD1 tkdry,
                  ArrayD1 thk)
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
      double satw = (h2osoi_liq(i) / DENH2O + h2osoi_ice(i) / DENICE) / (dz(i) * watsat(i));
      satw = std::min(1.0, satw);
      
      if (satw > 1.0e-6) {
        double dke;
        if (t_soisno(i) >= TFRZ) {                 // Unfrozen soil
          dke = std::max(0.0, std::log(satw) + 1.0);
        } else {                                   // Frozen soil
          dke = satw;
        }

        double fl = (h2osoi_liq(i) / (DENH2O * dz(i))) / (h2osoi_liq(i) / (DENH2O * dz(i)) +
             h2osoi_ice(i) / (DENICE * dz(i)));
        double dksat = tkmg(i) * pow(TKWAT, fl * watsat(i)) * pow(TKICE, (1.0 - fl) * watsat(i));
        thk(i) = dke * dksat + (1.0 - dke) * tkdry(i);  
      } else {
         thk(i) = tkdry(i);
      }

      if (i >= nlevsno + nlevbed) { thk(i) = TKBDRK; }

    } else if (ltype == istice || ltype == istice_mec) {

      thk(i) = TKWAT;
      if (t_soisno(i) < TFRZ) { thk(i) = TKICE; }

    } else if (ltype == istwet) {
      
      if (i >= nlevsno + nlevbed) { 
        thk(i) = TKBDRK;
      } else {
        thk(i) = TKWAT;
        if (t_soisno(i) < TFRZ) { thk(i) = TKICE; }
      }
    }
  }
}


template <typename ArrayD1>
ACCELERATE
void calc_snow_tk(const int& snl,
                  const double& frac_sno,
                  const ArrayD1 h2osoi_liq,
                  const ArrayD1 h2osoi_ice,
                  const ArrayD1 dz,
                  ArrayD1 thk)
{
  static constexpr double TKAIR{0.023};  // thermal conductivity of air   [W/m/K]
  using detail::TKICE;
  using ELMconst::TFRZ;
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;

  // zero out inactive layers
  const int top = nlevsno - snl;
  for (int i = 0; i < top; ++i) {
    thk(i) = 0.0;
  }

  // calculate for all active snow layers
  for (int i = top; i < nlevsno; ++i) {
    // partial density of water in the snow pack (ice + liquid) [kg/m3] 
    double bw = (h2osoi_ice(i) + h2osoi_liq(i)) / (frac_sno * dz(i));
    thk(i) = TKAIR + (7.75e-5 * bw + 1.105e-6 * bw * bw) * (TKICE - TKAIR);
  }
}


// Thermal conductivity at the layer interface
// tk[nlevgrnd+nlevsno]
// tk(i) is the interface between cells i and i+1
// this is different than zi, where zi(i) is between cells i-1 and i 
template <typename ArrayD1>
ACCELERATE
void calc_face_tk(const int& snl,
                  const ArrayD1 thk,
                  const ArrayD1 z,
                  const ArrayD1 zi,
                  ArrayD1 tk)
{
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;

  // zero out inactive interfaces
  const int top = nlevsno - snl;
  for (int i = 0; i < top; ++i) {
    tk(i) = 0.0;
  }

  // active interfaces above bottom interface
  for (int i = top; i < top - 1; ++i) {
    tk(i) = thk(i) * thk(i+1) * (z(i+1) - z(i)) / (thk(i) * (z(i+1) - zi(i+1)) + thk(i+1) * (zi(i+1) - z(i)));
  }

  // bottom interface
  tk(nlevgrnd + nlevsno - 1) = 0.0;
}


// Thermal conductivity of h2osfc
template <typename ArrayD1>
ACCELERATE
void calc_h2osfc_tk(const double& h2osfc,
                      const ArrayD1 thk,
                      const ArrayD1 z,
                      double& tk_h2osfc)
{
  using ELMdims::nlevsno;
  using detail::TKWAT;

  double zh2osfc = 1.0e-3 * (0.5 * h2osfc);  // convert to [m] from [mm]
  tk_h2osfc = TKWAT * thk(nlevsno) * (z(nlevsno) + zh2osfc) / (TKWAT * z(nlevsno) + thk(nlevsno) * zh2osfc);
}


// Soil heat capacity
template <typename ArrayD1>
ACCELERATE
void calc_soil_heat_capacity(const int& ltype,
                             const int& snl,
                             const double& h2osno,
                             const ArrayD1 watsat,
                             const ArrayD1 h2osoi_ice,
                             const ArrayD1 h2osoi_liq,
                             const ArrayD1 dz,
                             const ArrayD1 csol,
                             ArrayD1 cv)
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
      cv(i) = csol(i) * (1.0 - watsat(i-nlevsno)) * dz(i) + (h2osoi_ice(i) * CPICE + h2osoi_liq(i) * CPWAT);
    } else if (ltype == istwet) {
      cv(i) = (h2osoi_ice(i) * CPICE + h2osoi_liq(i) * CPWAT);
      if (i >= nlevsno + nlevbed) cv(i) = csol(i) * dz(i);
    } else if (ltype == istice || ltype == istice_mec) {
      cv(i) = (h2osoi_ice(i) * CPICE + h2osoi_liq(i) * CPWAT);
    }

    if (i == nlevsno && snl == 0 && h2osno > 0.0)
    {  cv(i) += CPICE * h2osno; }
  }
}


// Snow heat capacity
template <typename ArrayD1>
ACCELERATE
void calc_snow_heat_capacity(const int& snl,
                             const double& frac_sno,
                             const ArrayD1 h2osoi_ice,
                             const ArrayD1 h2osoi_liq,
                             ArrayD1 cv)
{
  static constexpr double thin_sfclayer{1.0e-6};   // Threshold for thin surface layer
  using ELMdims::nlevsno;

  using ELMconst::CPICE;
  using ELMconst::CPWAT;

   // zero out inactive layers
  const int top = nlevsno - snl;
  for (int i = 0; i < top; ++i) {
    cv(i) = 0.0;
  }

  // calculate for all active snow layers
  for (int i = top; i < nlevsno; ++i) {
    if (frac_sno < 0.0) {
      cv(i) = std::max(thin_sfclayer, (CPWAT * h2osoi_liq(i) + CPICE * h2osoi_ice(i)) / frac_sno);
    } else {
      cv(i) = thin_sfclayer;
    }
  }
}



} // namespace ELM::soil_thermal
