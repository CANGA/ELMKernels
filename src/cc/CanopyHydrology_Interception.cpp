#include <algorithm>
#include <stdio.h>     
#include <cmath>   
#include <iostream>
#include <string>
#include "landunit_varcon.h"
#include "column_varcon.h"     
#include "CanopyHydrology_cpp.hh"

using namespace std;
using std::min ;
using std::max ;

namespace ELM {
 NATURE void CanopyHydrology_Interception(double dtime,
        const double& forc_rain,
        const double& forc_snow,
        const double& irrig_rate,
        const int& ltype, const int& ctype,
        const bool& urbpoi, const bool& do_capsnow,
        const double& elai, const double& esai,
        const double& dewmx, const int& frac_veg_nosno,
        double& h2ocan,
        int n_irrig_steps_left, //fix it
        double& qflx_prec_intr,
        double& qflx_irrig,
        double& qflx_prec_grnd,
        double& qflx_snwcp_liq,
        double& qflx_snwcp_ice,
        double& qflx_snow_grnd_patch,
        double& qflx_rain_grnd)

 {  

  
      double  fpi, xrun, h2ocanmx   ;
      double  qflx_candrip, qflx_through_snow, qflx_through_rain ;
      double  qflx_prec_grnd_snow;
      double  qflx_prec_grnd_rain ;
      double  fracsnow ;
      double  fracrain , forc_irrig;


      if (ltype==istsoil || ltype==istwet || urbpoi || ltype==istcrop) {

         qflx_candrip = 0.0      ;
         qflx_through_snow = 0.0 ;
         qflx_through_rain = 0.0 ;
         qflx_prec_intr = 0.0    ;
         fracsnow = 0.0          ;
         fracrain = 0.0          ;
         forc_irrig = 0.0;


         if (ctype != icol_sunwall && ctype != icol_shadewall) {
            if (frac_veg_nosno == 1 && (forc_rain + forc_snow) > 0.0) {

              
               fracsnow = forc_snow/(forc_snow + forc_rain);
               fracrain = forc_rain/(forc_snow + forc_rain);

               
               h2ocanmx = dewmx * (elai + esai);

               
               fpi = 0.250*(1.0 - exp(-0.50*(elai + esai)));

              
               qflx_through_snow = forc_snow * (1.0-fpi);
               qflx_through_rain = forc_rain * (1.0-fpi);

               
               qflx_prec_intr = (forc_snow + forc_rain) * fpi;
               


               
               h2ocan = max(0.0, h2ocan + dtime*qflx_prec_intr);

               
               qflx_candrip = 0.0;

               
               xrun = (h2ocan - h2ocanmx)/dtime;

               
               if (xrun > 0.0) {
                  qflx_candrip = xrun;
                  h2ocan = h2ocanmx;
               }

            }
         }

      else if (ltype==istice || ltype==istice_mec) {
         
         h2ocan            = 0.0;
         qflx_candrip      = 0.0;
         qflx_through_snow = 0.0;
         qflx_through_rain = 0.0;
         qflx_prec_intr    = 0.0;
         fracsnow          = 0.0;
         fracrain          = 0.0;

      }

      

      if (ctype != icol_sunwall && ctype != icol_shadewall) {
         if (frac_veg_nosno == 0) {
            qflx_prec_grnd_snow = forc_snow;
            qflx_prec_grnd_rain = forc_rain;  }
         else{
            qflx_prec_grnd_snow = qflx_through_snow + (qflx_candrip * fracsnow);
            qflx_prec_grnd_rain = qflx_through_rain + (qflx_candrip * fracrain);
          }
      }   
      else{
         qflx_prec_grnd_snow = 0.;
         qflx_prec_grnd_rain = 0.;
        }

      
      if (n_irrig_steps_left > 0) {
         qflx_irrig         = forc_irrig;
         n_irrig_steps_left = n_irrig_steps_left - 1; }
      else{
         qflx_irrig = 0.0;
        }

      
      qflx_prec_grnd_rain = qflx_prec_grnd_rain + qflx_irrig;

      

      qflx_prec_grnd = qflx_prec_grnd_snow + qflx_prec_grnd_rain;

      if (do_capsnow) {
         qflx_snwcp_liq = qflx_prec_grnd_rain;
         qflx_snwcp_ice = qflx_prec_grnd_snow;

         qflx_snow_grnd_patch = 0.0;
         qflx_rain_grnd = 0.0;  }
      else{

         qflx_snwcp_liq = 0.0;
         qflx_snwcp_ice = 0.0;
         qflx_snow_grnd_patch = qflx_prec_grnd_snow   ;      //ice onto ground (mm/s)
         qflx_rain_grnd     = qflx_prec_grnd_rain      ;   //liquid water onto ground (mm/s)
        }

    }
  }

}