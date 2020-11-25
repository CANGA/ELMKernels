#include <netcdf.h>
#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <iostream>
#include "physics_dummy.hh"

int main(int argc, char** argv)
{

//gather and format all input data
    //read_input();

//call physics
    physics_dummy();

    /* Calling tree currently looks like:
    read input, process data and initialize
    ReadInput() - need to write
    TopoSlopes() - set 0.2 (-) as minimum slope
    ColdStartMicroTopo() - sets n_melt and micro_sigma for SCA calculations
    ColdStartSnow() - sets number of snow layers (snl), snow dz, z, zi


    Time stepping:
    each of these functions operate on a single non-lake grid cell with !one! PFT
    CanopyInterception(); 
    CanopyIrrigation();
    CanopyGroundFlux();
    Canopy FracWet();
    SnowInit();
    FracH2OSfc();
    CanopySunShadeFrac();

    */



return 0;
}

