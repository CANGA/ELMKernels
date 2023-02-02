// functions from SurfaceResistanceMod.F90

#pragma once

namespace ELM::surface_resistance {

template <typename ArrayD1>
ACCELERATE
void calc_soilevap_stress(const LandType& Land, const double& frac_sno, const double& frac_h2osfc,
                          const ArrayD1 watsat, const ArrayD1 watfc, const ArrayD1 h2osoi_liq,
                          const ArrayD1 h2osoi_ice, const ArrayD1 dz, double& soilbeta)
{
  using ELMdims::nlevsno;
  
  if (!Land.lakpoi) {

    // local variables
    double fac, fac_fc, wx;

    if (Land.ltype != LND::istwet && Land.ltype != LND::istice && Land.ltype != LND::istice_mec) {
      if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
        wx = (h2osoi_liq(nlevsno()) / ELMconst::DENH2O() + h2osoi_ice(nlevsno()) / ELMconst::DENICE()) / dz(nlevsno());
        fac = std::min(1.0, wx / watsat(0));
        fac = std::max(fac, 0.01);
        // Lee and Pielke 1992 beta, added by K.Sakaguchi
        if (wx < watfc(0)) {                     //! when water content of ths top layer is less than that at F.C.
          fac_fc = std::min(1.0, wx / watfc(0)); // eqn5.66 but divided by theta at field capacity
          fac_fc = std::max(fac_fc, 0.01);
          // modify soil beta by snow cover. soilbeta for snow surface is one
          soilbeta =
              (1.0 - frac_sno - frac_h2osfc) * 0.25 * pow(1.0 - cos(ELMconst::ELM_PI() * fac_fc), 2.0) + frac_sno + frac_h2osfc;
        } else {
          soilbeta = 1.0;
        }
      } else if (Land.ltype == LND::icol_road_perv) {
        soilbeta = 0.0;
      } else if (Land.ltype == LND::icol_sunwall || Land.ltype == LND::icol_shadewall) {
        soilbeta = 0.0;
      } else if (Land.ltype == LND::icol_roof || Land.ltype == LND::icol_road_imperv) {
        soilbeta = 0.0;
      }
    } else {
      soilbeta = 1.0;
    }
  }
}

ACCELERATE
double getlblcef(const double& rho, const double& temp)
{
  static constexpr double C{120.0};      // K
  static constexpr double T0{291.25};    // K
  static constexpr double mu0{18.27e-6}; // Pa s
  static constexpr double prandtl{0.72};
  // compute the kinetic viscosity
  double mu = mu0 * (T0 + C) / (temp + C) * pow(temp / T0, 1.5) / rho; // m^2 s^-1
  double diffh2o = 0.229e-4 * pow(temp / 273.15, 1.75);                // m^2 s^-1
  double sc = mu / diffh2o;                                            // schmidt number
  double result = 2.0 / ELMconst::VKC() * pow(sc / prandtl, 2.0 / 3.0);
  return result;
}

} // namespace ELM::surface_resistance
