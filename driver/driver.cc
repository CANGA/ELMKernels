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
    CanopyInterception(); 
    CanopyIrrigation();
    CanopyGroundFlux();
    Canopy FracWet();
    */



return 0;
}

