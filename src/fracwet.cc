/*
DESCRIPTION:
Determine fraction of vegetated surfaces which are wet and
fraction of elai which is dry. The variable ``fwet'' is the
fraction of all vegetation surfaces which are wet including
stem area which contribute to evaporation. The variable ``fdry''
is the fraction of elai which is dry because only leaves
can transpire.  Adjusted for stem area which does not transpire.

INPUTS:
frac_veg_nosno [int]    fraction of veg not covered by snow (0/1 now) [-]
dewmx          [double] Maximum allowed dew [mm]                
elai           [double] one-sided leaf area index with burying by snow
esai           [double] one-sided stem area index with burying by snow
h2ocan         [double] total canopy water (mm H2O)

OUTPUTS:          
fwet           [double] fraction of canopy that is wet (0 to 1) 
fdry           [double] fraction of foliage that is green and dry [-] (new)

*/
#include <algorithm>
#include <cmath>

void FracWet(const int& frac_veg_nosno,
    const double& dewmx,
    const double& elai,
    const double& esai,
    const double& h2ocan,
    double& fwet,
    double& fdry)

{
    if (frac_veg_nosno == 1) {
        if (h2ocan > 0.0) {
            double vegt = frac_veg_nosno * (elai + esai);
            double dewmxi = 1.0 / dewmx;

            fwet = pow(((dewmxi/vegt)*h2ocan), 2.0/3.0);
            fwet = std::min(fwet, 1.0);  //Check for maximum limit of fwet
        } else {
            fwet = 0.0;
        }
        fdry = (1.0 - fwet) * elai / (elai + esai);
    } else {
        fwet = 0.0;
        fdry = 0.0;
    }
}
