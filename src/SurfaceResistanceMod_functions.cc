// functions from SurfaceResistanceMod.F90

#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

/* calc_soilevap_stress()
DESCRIPTION:
Compute the Lee-Pielke beta factor to scale actual soil evaporation from potential evaporation.
This is the calc_beta_leepielke1992() function from CLM/ELM that gets called from calc_soilevap_stress().

code that gets called when use_vsfm == true is not currently included 

INPUTS:
ltype                        [int]    landunit type                            
lakpoi                       [bool]   true => landunit is an urban point
frac_sno                     [double] fraction of ground covered by snow (0 to 1)
frac_h2osfc                  [double] fraction of ground covered by surface water (0 to 1)
watsat[nlevgrnd]             [double] volumetric soil water at saturation (porosity)
watfc[nlevgrnd]              [double] volumetric soil water at field capacity
h2osoi_liq[nlevgrnd+nlevsno] [double] ice lens (kg/m2)                       
h2osoi_ice[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)                  

OUTPUTS:
soilbeta [double] factor that reduces ground evaporation
*/
void calc_soilevap_stress(
  const int& ltype,
  const bool& lakpoi,
  const double& frac_sno,
  const double& frac_h2osfc,
  const double* watsat_,
  const double* watfc_,
  const double* h2osoi_liq_,
  const double* h2osoi_ice_,
  const double* dz_,
  

  double& soilbeta)
{
  // reference variables with pointer offset to allow use of CLM/ELM index without modification
  const double* watsat = &watsat_[-1];
  const double* watfc = &watfc_[-1];
  const double* h2osoi_liq = &h2osoi_liq_[nlevsno-1];
  const double* h2osoi_ice = &h2osoi_ice_[nlevsno-1];
  const double* dz = &dz_[nlevsno-1];

  // local variables
  double fac, fac_fc, wx;
  const double denice = SHR_CONST_RHOICE;
  const double denh2o = SHR_CONST_RHOFW;

  if (!lakpoi) {
    if (ltype != istwet && ltype != istice && ltype != istice_mec) {
      if (ltype == istsoil || ltype == istcrop) {
        wx   = (h2osoi_liq[1] / denh2o + h2osoi_ice[1] / denice) / dz[1];
        fac  = std::min(1.0, wx / watsat[1]);
        fac  = std::max(fac, 0.01);
        // Lee and Pielke 1992 beta, added by K.Sakaguchi
        if (wx < watfc[1] ) { //!when water content of ths top layer is less than that at F.C.
          fac_fc  = std::min(1.0, wx / watfc[1]);  // eqn5.66 but divided by theta at field capacity
          fac_fc  = std::max(fac_fc, 0.01);
          // modify soil beta by snow cover. soilbeta for snow surface is one
          soilbeta = (1.0 - frac_sno - frac_h2osfc) * 0.25 * pow(1.0 - cos(SHR_CONST_PI * fac_fc), 2.0)
          + frac_sno + frac_h2osfc;
        } else {
          soilbeta = 1.0;
        }
      } else if (ltype == icol_road_perv) {
        soilbeta = 0.0;
      } else if (ltype == icol_sunwall || ltype == icol_shadewall) {
        soilbeta = 0.0;
      } else if (ltype == icol_roof || ltype == icol_road_imperv) {
        soilbeta = 0.0;
      }
    } else {
      soilbeta = 1.0;
    }
  }
}
