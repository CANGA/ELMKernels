// functions from SurfaceResistanceMod.F90

#include "clm_constants.h"
#include "landtype.h"
#include <algorithm>
#include <cmath>

namespace ELM {

/* calc_soilevap_stress()
DESCRIPTION:
Compute the Lee-Pielke beta factor to scale actual soil evaporation from potential evaporation.
This is the calc_beta_leepielke1992() function from CLM/ELM that gets called from calc_soilevap_stress().

code that gets called when use_vsfm == true is not currently included

INPUTS:
Land                         [LandType] struct containing information about landtypet
frac_sno                     [double] fraction of ground covered by snow (0 to 1)
frac_h2osfc                  [double] fraction of ground covered by surface water (0 to 1)
watsat[nlevgrnd]             [double] volumetric soil water at saturation (porosity)
watfc[nlevgrnd]              [double] volumetric soil water at field capacity
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
soilbeta [double] factor that reduces ground evaporation
*/
void calc_soilevap_stress(const LandType &Land, const double &frac_sno, const double &frac_h2osfc, const double *watsat,
                          const double *watfc, const double *h2osoi_liq, const double *h2osoi_ice, const double *dz,

                          double &soilbeta) {
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

/* getlblcef()
!compute the scaling paramter for laminar boundary resistance
     !the laminar boundary layer resistance is formulated as
     !Rb=2/(k*ustar)*(Sci/Pr)^(2/3)
     !cc = Rb*ustar
     !   = 2/k*(Sci/Pr)^(2/3)
     !   Pr=0.72, Prandtl number
     !   Sci = v/Di, Di is diffusivity of gas i
     !   v : kinetic viscosity
*/
double getlblcef(const double &rho, const double &temp) {
  const double C = 120.0;      // K
  const double T0 = 291.25;    // K
  const double mu0 = 18.27e-6; // Pa s
  const double prandtl = 0.72;
  // compute the kinetic viscosity
  double mu = mu0 * (T0 + C) / (temp + C) * pow(temp / T0, 1.5) / rho; // m^2 s^-1
  double diffh2o = 0.229e-4 * pow(temp / 273.15, 1.75);                // m^2 s^-1
  double sc = mu / diffh2o;                                            // schmidt number
  double result = 2.0 / vkc * pow(sc / prandtl, 2.0 / 3.0);
  return result;
}

} // namespace ELM
