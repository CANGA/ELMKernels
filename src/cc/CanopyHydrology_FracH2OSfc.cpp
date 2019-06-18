#include <algorithm>
#include <stdio.h>     
#include <cmath>   
#include <iostream>
#include <string>
#include "landunit_varcon.h"
#include "CanopyHydrology_cpp.hh"

using namespace std;
using std::min ;
using std::max ;

namespace ELM {

void CanopyHydrology_FracH2OSfc(const double& dtime,
        const double& min_h2osfc,
        const int& ltype,
        const double& micro_sigma,
        const double& h2osno,
        double& h2osfc,
        double& h2osoi_liq,
        double& frac_sno,
        double& frac_sno_eff,
        double& qflx_h2osfc2topsoi,
        double& frac_h2osfc)
  {
    bool no_update = false;
    double shr_const_pi=4.0e0*atan(1.0e0) ;
    bool no_update_l ;

    
    double d,fd,dfdd,sigma   ;

    if (!no_update) { 
      no_update_l = false; }
    else { no_update_l = no_update; }
    
    qflx_h2osfc2topsoi = 0.0 ;
    
    if ( ltype  == istsoil || ltype == istcrop) {

       

       if (h2osfc > min_h2osfc) {
          
          d=0.0 ;

          sigma=1.0e3 * micro_sigma ;
          for(int l = 0 ; l < 10; l++) {
             fd = 0.5*d*(1.00+erf(d/(sigma*sqrt(2.0)))) + sigma/sqrt(2.0*shr_const_pi)*exp(-pow(d,2)/(2.0*pow(sigma,2))) -h2osfc;
             dfdd = 0.5*(1.00+erf(d/(sigma*sqrt(2.0))));

             d = d - fd/dfdd;
          }
          
          frac_h2osfc = 0.5*(1.00+erf(d/(sigma*sqrt(2.0)))) ;  }

       else {
          frac_h2osfc = 0.0 ;
          h2osoi_liq = h2osoi_liq + h2osfc ;
          qflx_h2osfc2topsoi = h2osfc/dtime ;
          h2osfc=0.0 ;
        }

       if (!no_update_l) {

          
          if (frac_sno > (1.0 - frac_h2osfc) && h2osno > 0) {

             if (frac_h2osfc > 0.010) {
                frac_h2osfc = max(1.00 - frac_sno,0.010) ;
                frac_sno = 1.00 - frac_h2osfc; }
             else {
                frac_sno = 1.00 - frac_h2osfc;
              }
             
             frac_sno_eff=frac_sno;

          }

       } 
    }  
    else {

       frac_h2osfc = 0.0;

    }

  }

}