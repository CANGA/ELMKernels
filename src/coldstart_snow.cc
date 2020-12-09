/*
//from initVerticalMod
DESCRIPTION:
Choose number of snow layers (max 5) based on depth of snow. Initialize snow layer thickness, depth, interface depth.

CLM uses negative indexing for the snow model - snl is the negative of the number of snow layers - i=0 is layer next to soil, i = snl+1 is top layer.
Here we use a positive snl - i=0 is layer next to soil, and i=snl-1 is top layer.

INPUTS:
snow_depth [double] snow height (m)
lakpoi     [bool] true => landunit is a lake point

OUTPUTS:
snl [integer] number of snow layers
dz  [double]  layer thickness (m)
z   [double]  layer depth (m)
zi  [double]  interface level below a "z" level (m) 
*/

void ColdStartSnow(
   const double& snow_depth,
   const bool& lakpoi,

   int& snl,
   double* dz,
   double* z,
   double* zi)
{

   if (!lakpoi) {
      if (snow_depth < 0.01) {
         snl = 0;
      } else {
         if ((snow_depth >= 0.01) && (snow_depth <= 0.03)) {
            snl   = 1;
            dz[0] = snow_depth;
         } else if ((snow_depth > 0.03) && (snow_depth <= 0.04)) {
            snl   = 2;
            dz[1] = snow_depth / 2.0;
            dz[0] = dz[1];
         } else if ((snow_depth > 0.04) && (snow_depth <= 0.07)) {
            snl   = 2;
            dz[1] = 0.02;
            dz[0] = snow_depth - dz[1];
         } else if ((snow_depth > 0.07) && (snow_depth <= 0.12)) {
            snl   = 3;
            dz[2] = 0.02;
            dz[1] = (snow_depth - 0.02) / 2.0;
            dz[0] = dz[1];
         } else if ((snow_depth > 0.12) && (snow_depth <= 0.18)) {
            snl   = 3;
            dz[2] = 0.02;
            dz[1] = 0.05;
            dz[0] = snow_depth - dz[2] - dz[1];
         } else if ((snow_depth > 0.18) && (snow_depth <= 0.29)) {
            snl   = 4;
            dz[3] = 0.02;
            dz[2] = 0.05;
            dz[1] = (snow_depth - dz[3] - dz[2]) / 2.0;
            dz[0] = dz[1];
         } else if ((snow_depth > 0.29) && (snow_depth <= 0.41)) {
            snl   = 4;
            dz[3] = 0.02;
            dz[2] = 0.05;
            dz[1] = 0.11;
            dz[0] = snow_depth - dz[3] - dz[2] - dz[1];
         } else if ((snow_depth > 0.41) && (snow_depth <= 0.64)) {
            snl   = 5;
            dz[4] = 0.02;
            dz[3] = 0.05;
            dz[2] = 0.11;
            dz[1] = (snow_depth - dz[4] - dz[3] - dz[2]) / 2.0;
            dz[0] = dz[1];
         } else if (snow_depth > 0.64) {
            snl   = 5;
            dz[4] = 0.02;
            dz[3] = 0.05;
            dz[2] = 0.11;
            dz[1] = 0.23;
            dz[0] = snow_depth - dz[4] - dz[3] - dz[2] - dz[1];
         }
      }
      for (int j = 0; j < snl; j++) {
         // calc z (cell centers - 0,..,snl-1) and zi (cell interface elevations - 0,..,snl)
         // zi is initialized as 0.0, which will always be the reference elevation of the ground/snow interface
         // zi[0] will always equal 0.0, and zi[1,..,snl] are calculated below 
         z[j]    = zi[j] - 0.5 * dz[j];
         zi[j+1] = zi[j] - dz[j];
      }
   } else {
      snl = 0;
   }
}
