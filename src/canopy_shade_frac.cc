/*
DESCRIPTION:
This subroutine calculates and returns:
1) absorbed PAR for sunlit leaves in canopy layer
2) absorbed PAR for shaded leaves in canopy layer
3) sunlit leaf area
4) shaded  leaf area
5) sunlit leaf area for canopy layer
6) shaded leaf area for canopy layer
7) sunlit fraction of canopy

INPUTS:
nrad       [int] number of canopy layers
elai       [double] one-sided leaf area index
tlai_z     [double] tlai increment for canopy layer
fsun_z     [double] sunlit fraction of canopy layer
forc_solad [double] direct beam radiation (W/m**2) 
forc_solai [double] diffuse radiation (W/m**2)
fabd_sun_z [double] absorbed sunlit leaf direct PAR 
fabd_sha_z [double] absorbed shaded leaf direct PAR
fabi_sun_z [double] absorbed sunlit leaf diffuse PAR
fabi_sha_z [double] absorbed shaded leaf diffuse PAR

OUTPUTS:
parsun_z   [double] absorbed PAR for sunlit leaves
parsha_z   [double] absorbed PAR for shaded leaves
laisun_z   [double] sunlit leaf area for canopy layer
laisha_z   [double] shaded leaf area for canopy layer
laisun,    [double] sunlit leaf area    
laisha,    [double] shaded  leaf area  
fsun)      [double] sunlit fraction of canopy
*/

void CanopySunShadeFrac(SurfAlb_type& surfalb_type)
 //const int& nrad,
 //const double& elai,
 //const double* tlai_z,
 //const double* fsun_z,
 //const double* forc_solad,
 //const double* forc_solai,
 //const double* fabd_sun_z,
 //const double* fabd_sha_z,
 //const double* fabi_sun_z,
 //const double* fabi_sha_z,
 //

 //double* parsun_z,
 //double* parsha_z,
 //double* laisun_z,
 //double* laisha_z,
 //double& laisun,  
 //double& laisha,
 //double& fsun
{
  int ipar = 1; // The band index for PAR
  for (int iv = 0; iv < nrad; iv++) {
    parsun_z[iv] = 0.0;
    parsha_z[iv] = 0.0;
    laisun_z[iv] = 0.0;
    laisha_z[iv] = 0.0;
  }
  // Loop over patches to calculate laisun_z and laisha_z for each layer.
  // Derive canopy laisun, laisha, and fsun from layer sums.
  // If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
  // SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
  laisun = 0.0;
  laisha = 0.0;

  for (int iv = 0; iv < nrad; iv++) {
    laisun_z[iv] = tlai_z[iv] * fsun_z[iv];
    laisha_z[iv] = tlai_z[iv] * (1.0 - fsun_z[iv]);
    laisun += laisun_z[iv];
    laisha += laisha_z[iv];     
  }
  if (elai > 0.0) {
    fsun = laisun / elai;
  } else {
    fsun = 0.0;
  }
  // Absorbed PAR profile through canopy
  // If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
  // are canopy integrated so that layer values equal big leaf values.
  for (int iv = 0; iv < nrad; iv++) {
    parsun_z[iv] = forc_solad[ipar] * fabd_sun_z[iv] + forc_solai[ipar] * fabi_sun_z[iv];
    parsha_z[iv] = forc_solad[ipar] * fabd_sha_z[iv] + forc_solai[ipar] * fabi_sha_z[iv];
  }
}
