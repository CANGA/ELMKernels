#include <algorithm>
#include <iostream>
#include <cmath>
#include <string>
#include "CanopyHydrology_cpp.hh"

using namespace std;
using std::min ;
using std::max ;

namespace ELM {

NATURE void CanopyHydrology_FracWet(const int& frac_veg_nosno,
        const double& h2ocan,
        const double& elai, 
        const double& esai,
        const double& dewmx,
        double& fwet,
        double& fdry)
{  

  double vegt, dewmxi ;    

  if (frac_veg_nosno == 1) {
    if (h2ocan > 0.0) {
        vegt    = frac_veg_nosno*(elai + esai);
        dewmxi  = 1.00/dewmx;
        fwet = pow(((dewmxi/vegt)*h2ocan), 2.0/3);
        fwet = min(fwet,1.00);   
        fdry = (1.0-fwet)*elai/(elai+esai);
      }
    else{
        fwet = 0.0;
        fdry = 0.0 ;
      } }
  else{
     fwet = 0.0;
     fdry = 0.0;
  }

}

}