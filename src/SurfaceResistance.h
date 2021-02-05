// functions from SurfaceResistanceMod.F90

#pragma once

#include "SurfaceResistance_impl.hh"
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
template <class dArray_type>
void calc_soilevap_stress(const LandType &Land, const double &frac_sno, const double &frac_h2osfc,
                          const dArray_type watsat, const dArray_type watfc, const dArray_type h2osoi_liq,
                          const dArray_type h2osoi_ice, const dArray_type dz, double &soilbeta);

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
double getlblcef(const double &rho, const double &temp);

} // namespace ELM
