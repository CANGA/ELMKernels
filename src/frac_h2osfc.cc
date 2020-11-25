/*
DESCRIPTION:
Determine fraction of land surfaces which are submerged  
based on surface microtopography and surface water storage.

!!h2osoi_liq - indexing in CLM, (-4,..,15) (top of snow,..,bottom of bedrock) -- for now, I decouple 
the snow layers CLM(-4,..,0), this code (5,..,0) from the subsurface CLM(1,..15) this code (0,..14)



INPUTS:
dtime              [double] land model time step (sec)
ltype              [int] landunit type
micro_sigma        [double] microtopography pdf sigma (m)
h2osno             [double] snow water (mm H2O)
no_update          [bool] flag to make calculation w/o updating variables

OUTPUTS:
h2osfc             [double] surface water (mm)
h2osoi_liq         [double] liquid water (kg/m2) --  just the top subsurface layer for this function
frac_sno           [double] fraction of ground covered by snow (0 to 1)
frac_sno_eff       [double] effective fraction of ground covered by snow (0 to 1)
qflx_h2osfc2topsoi [double] liquid water coming from surface standing water top soil (mm H2O/s)
frac_h2osfc        [double] fractional area with surface water greater than zero (0 to 1)


*/

#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

void FracH2OSfc(
  const double& dtime,
  const int& ltype,
  const double& micro_sigma,
  const double& h2osno,
  const bool& no_update,

  double& h2osfc,
  double& h2osoi_liq,
  double& frac_sno,
  double& frac_sno_eff,
  double& qflx_h2osfc2topsoi,
  double& frac_h2osfc)

{
  double d,fd,dfdd,sigma;
  // arbitrary lower limit on h2osfc for safer numerics...
  double min_h2osfc = 1.e-8;
  qflx_h2osfc2topsoi = 0.0;

  // h2osfc only calculated for soil vegetated land units
  if (ltype == istsoil || ltype == istcrop) {
    // Use newton-raphson method to iteratively determine frac_h20sfc
    // based on amount of surface water storage (h2osfc) and 
    // microtopography variability (micro_sigma)
    if (h2osfc > min_h2osfc) {
      // a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
      d = 0.0;
      sigma = 1.0e3 * micro_sigma; // convert to mm

      for (int l = 0; l < 10; l++) {
        fd = 0.5 * d * (1.0 + erf(d / (sigma * sqrt(2.0)))) + sigma / sqrt(2.0 * SHR_CONST_PI) * exp(-pow(d, 2) / (2.0 * pow(sigma, 2))) - h2osfc;
        dfdd = 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))));
        d = d - fd/dfdd;
      }
      // update the submerged areal fraction using the new d value
      frac_h2osfc = 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))));
    } else {
      frac_h2osfc = 0.0;
      h2osoi_liq = h2osoi_liq + h2osfc;
      qflx_h2osfc2topsoi = h2osfc / dtime;                
      h2osfc = 0.0;
    }
    if (!no_update) {
      // adjust fh2o, fsno when sum is greater than zero
      if (frac_sno > (1.0 - frac_h2osfc) && h2osno > 0.0) {
        if (frac_h2osfc > 0.01) {
          frac_h2osfc = std::max((1.0 - frac_sno), 0.01);
          frac_sno = 1.0 - frac_h2osfc;
        } else {
          frac_sno = 1.0 - frac_h2osfc;
        }
        frac_sno_eff = frac_sno;
      }
    }
  } else {
    frac_h2osfc = 0.0;
  }
}
