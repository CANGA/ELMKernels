// functions from SurfaceResistanceMod.F90

#pragma once

namespace ELM {
namespace surface_resistance {

template <class ArrayD1>
void calc_soilevap_stress(const LandType &Land, const double &frac_sno, const double &frac_h2osfc,
                          const ArrayD1 watsat, const ArrayD1 watfc, const ArrayD1 h2osoi_liq,
                          const ArrayD1 h2osoi_ice, const ArrayD1 dz, double &soilbeta) {
  if (!Land.lakpoi) {

    // local variables
    double fac, fac_fc, wx;

    if (Land.ltype != istwet && Land.ltype != istice && Land.ltype != istice_mec) {
      if (Land.ltype == istsoil || Land.ltype == istcrop) {
        wx = (h2osoi_liq[nlevsno] / denh2o + h2osoi_ice[nlevsno] / denice) / dz[nlevsno];
        fac = std::min(1.0, wx / watsat[0]);
        fac = std::max(fac, 0.01);
        // Lee and Pielke 1992 beta, added by K.Sakaguchi
        if (wx < watfc[0]) {                     //! when water content of ths top layer is less than that at F.C.
          fac_fc = std::min(1.0, wx / watfc[0]); // eqn5.66 but divided by theta at field capacity
          fac_fc = std::max(fac_fc, 0.01);
          // modify soil beta by snow cover. soilbeta for snow surface is one
          soilbeta =
              (1.0 - frac_sno - frac_h2osfc) * 0.25 * pow(1.0 - cos(ELM_PI * fac_fc), 2.0) + frac_sno + frac_h2osfc;
        } else {
          soilbeta = 1.0;
        }
      } else if (Land.ltype == icol_road_perv) {
        soilbeta = 0.0;
      } else if (Land.ltype == icol_sunwall || Land.ltype == icol_shadewall) {
        soilbeta = 0.0;
      } else if (Land.ltype == icol_roof || Land.ltype == icol_road_imperv) {
        soilbeta = 0.0;
      }
    } else {
      soilbeta = 1.0;
    }
  }
}

} // namespace surface_resistance
} // namespace ELM
