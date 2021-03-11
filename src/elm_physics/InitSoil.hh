// InitCold() from SoilStateType.F90
// sets watsat, sucsat, bsw, etc.

#include <RootBioPhys.hh>

namespace ELM {

void InitSoil() {
  double zsoifl[nlevsoi], zisoifl[nlevsoi + 1], dzsoifl[nlevsoi];

  smpmin = -1.0e8;

  if (urbpoi && ctype == icol_road_perv) {
    for (int i = 0; i < nlevgrnd; ++i) {
      rootfr_road_perv[i] = 0.0;
    }
    for (int i = 0; i < nlevsoi; ++i) {
      rootfr_road_perv[i] = 0.1;
    }
  }

  for (int i = nlevsoi; i < nlevgrnd; ++i) {
    rootfr[i] = 0.0;
  }

  init_vegrootfr(Land.vtype, rootfr);

  // --------------------------------------------------------------------
  // get original soil depths to be used in interpolation of sand and clay
  // --------------------------------------------------------------------
  for (int i = 0; i < nlevsoi; ++i) {
    zsoifl[i] = 0.025 * (exp(0.5 * (i - 0.5)) - 1.0);
  }                                           // node depths
  dzsoifl[0] = 0.5 * (zsoifl[0] + zsoifl[1]); // thickness b/n two interfaces
  for (int i = 1; i < nlevsoi - 1; ++i) {
    dzsoifl[i] = 0.5 * (zsoifl[i + 1] - zsoifl[i - 1]);
  }
  dzsoifl[nlevsoi - 1] = zsoifl[nlevsoi - 1] - zsoifl[nlevsoi - 2];
  zisoifl[0] = 0.0;
  for (int i = 1; i < nlevsoi; ++i) {
    zisoifl[i] = 0.5 * (zsoifl[i - 1] + zsoifl[i]);
  } // interface depths
  zisoifl[nlevsoi] = zsoifl[nlevsoi - 1] + 0.5 * dzsoifl[nlevsoi - 1];

  // --------------------------------------------------------------------
  // Set soil hydraulic and thermal properties: non-lake
  // --------------------------------------------------------------------
  // urban roof, sunwall and shadewall thermal properties used to
  // derive thermal conductivity and heat capacity are set to special
  // value because thermal conductivity and heat capacity for urban
  // roof, sunwall and shadewall are prescribed in SoilThermProp.F90
  // in SoilPhysicsMod.F90 -- I don't think those files exist
}

} // namespace ELM