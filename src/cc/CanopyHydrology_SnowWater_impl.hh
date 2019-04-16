#include <stdio.h>  
#include <math.h>     
#include "landunit_varcon.h"    
#include "column_varcon.h" 
#include "clm_varpar.h"         
#include "clm_varctl.h"     
#include "CanopyHydrology_cc.hh"    

using namespace std;

namespace ELM {

template<typename Array_d>
void CanopyHydrology_SnowWater(const double& dtime,
        const double& qflx_floodg,
        const int& ltype,
        const int& ctype,
        const bool& urbpoi,
        const bool& do_capsnow,                            
        const int& oldfflag,
        const double& forc_air_temp,
        const double& t_grnd,
        const double& qflx_snow_grnd_col,
        const double& qflx_snow_melt,
        const double& n_melt,
        const double& frac_h2osfc,
        double& snow_depth,
        double& h2osno,
        double& integrated_snow,
        Array_d swe_old,
        Array_d h2osoi_liq,
        Array_d h2osoi_ice,
        Array_d t_soisno,
        Array_d frac_iceold,
        int& snow_level,
        Array_d dz,
        Array_d z,
        Array_d zi,
        int& newnode,
        double& qflx_floodc,
        double& qflx_snow_h2osfc,
        double& frac_sno_eff,
        double& frac_sno)
{       
  
  
//parameters
  double rpi=4.0e0*atan(1.0e0)  ;
  double tfrz=273.15;
  double zlnd = 0.010;

  
  // real(r8), intent(inout), dimension(-nlevsno+1:0)  :: swe_old 
  // real(r8), intent(inout), dimension(-nlevsno+1:0) :: h2osoi_liq, h2osoi_ice
  // real(r8), intent(inout), dimension(-nlevsno+1:0)  :: t_soisno, frac_iceold
  // real(r8), intent(inout), dimension(-nlevsno+1:0)  :: dz, z, zi
  
//local variables 
  double  temp_intsnow, temp_snow_depth, z_avg, fmelt, dz_snowf, snowmelt ;
  double  newsnow, bifall, accum_factor, fsno_new, smr ;
  int j ;

//apply gridcell flood water flux to non-lake columns
  if (ctype != icol_sunwall && ctype != icol_shadewall) {      
     qflx_floodc = qflx_floodg; }
  else{
     qflx_floodc = 0.0;
  }

//Determine snow height and snow water

//Use Alta relationship, Anderson(1976); LaChapelle(1961),
//U.S.Department of Agriculture Forest Service, Project F,
//Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

  qflx_snow_h2osfc = 0.0;
//set temporary variables prior to updating
  temp_snow_depth=snow_depth;
//save initial snow content
  for(j = -nlevsno+1; j < snow_level; j++) {
     swe_old[j] = 0.00;
  }
  for(j = snow_level+1; j < 0; j++) {
     swe_old[j]=h2osoi_liq[j]+h2osoi_ice[j];
  }

  if (do_capsnow) {
     dz_snowf = 0.;
     newsnow = (1. - frac_h2osfc) * qflx_snow_grnd_col * dtime;
     frac_sno=1.;
     integrated_snow = 5.e2; 
  } else {

     if (forc_air_temp > tfrz + 2.) {
        bifall=50. + 1.7*std::pow((17.0),1.5);
     } else if (forc_air_temp > tfrz - 15.) {
        bifall=50. + 1.7*std::pow((forc_air_temp - tfrz + 15.),1.5);
     } else {
        bifall=50.;
     }

     // newsnow is all snow that doesn't fall on h2osfc
     newsnow = (1. - frac_h2osfc) * qflx_snow_grnd_col * dtime;

     // update integrated_snow
     integrated_snow = std::max(integrated_snow,h2osno) ; //h2osno could be larger due to frost

     // snowmelt from previous time step * dtime
     snowmelt = qflx_snow_melt * dtime;

     // set shape factor for accumulation of snow
     accum_factor=0.1;

     if (h2osno > 0.0) {

        //======================  FSCA PARAMETERIZATIONS  ======================
        // fsca parameterization based on *changes* in swe
        // first compute change from melt during previous time step
        if(snowmelt > 0.) {

           smr=std::min(1.,(h2osno)/(integrated_snow));

           frac_sno = 1. - std::pow((acos(min(1.,(2.*smr - 1.)))/rpi),(n_melt)) ;

        }

        // update fsca by new snow event, add to previous fsca
        if (newsnow > 0.) {
           fsno_new = 1. - (1. - tanh(accum_factor*newsnow))*(1. - frac_sno);
           frac_sno = fsno_new;

           // reset integrated_snow after accumulation events
           temp_intsnow= (h2osno + newsnow) / (0.5*(cos(rpi*std::pow((1.0-std::max(frac_sno,1e-6)),(1.0/n_melt))+1.0))) ;
           integrated_snow = std::min(1.e8,temp_intsnow) ;
        }

        //====================================================================

        // for subgrid fluxes
        if (subgridflag ==1 && ! urbpoi) {
           if (frac_sno > 0.){
              snow_depth=snow_depth + newsnow/(bifall * frac_sno);
           } else {
              snow_depth=0.;
           }
        } else {
           // for uniform snow cover
           snow_depth=snow_depth+newsnow/bifall;
        }

        // use original fsca formulation (n&y 07)
        if (oldfflag == 1) { 
           // snow cover fraction in Niu et al. 2007
           if(snow_depth > 0.0)  {
              frac_sno = tanh(snow_depth/(2.5*zlnd*std::pow((std::min(800.0,(h2osno+ newsnow)/snow_depth)/100.0),1.0)) ) ;
           }
           if(h2osno < 1.0)  {
              frac_sno=std::min(frac_sno,h2osno);
           }
        }

     } else { //h2osno == 0
        // initialize frac_sno and snow_depth when no snow present initially
        if (newsnow > 0.) { 
           z_avg = newsnow/bifall;
           fmelt=newsnow;
           frac_sno = tanh(accum_factor*newsnow);

           // make integrated_snow consistent w/ new fsno, h2osno
           integrated_snow = 0. ;//reset prior to adding newsnow below
           temp_intsnow= (h2osno + newsnow) / (0.5*(cos(rpi*std::pow((1.0-std::max(frac_sno,1e-6)),(1.0/n_melt)))+1.0));
           integrated_snow = std::min(1.e8,temp_intsnow);

           // update snow_depth and h2osno to be consistent with frac_sno, z_avg
           if (subgridflag ==1 && !urbpoi) {
              snow_depth=z_avg/frac_sno;
           } else {
              snow_depth=newsnow/bifall;
           }
           // use n&y07 formulation
           if (oldfflag == 1) { 
              // snow cover fraction in Niu et al. 2007
              if(snow_depth > 0.0)  {
                 frac_sno = tanh(snow_depth/(2.5*zlnd*std::pow((std::min(800.0,newsnow/snow_depth)/100.0),1.0)) );
              }
           }
        } else {
           z_avg = 0.;
           snow_depth = 0.;
           frac_sno = 0.;
        }
     } // end of h2osno > 0

     // snow directly falling on surface water melts, increases h2osfc
     qflx_snow_h2osfc = frac_h2osfc*qflx_snow_grnd_col;

     // update h2osno for new snow
     h2osno = h2osno + newsnow ;
     integrated_snow = integrated_snow + newsnow;

     // update change in snow depth
     dz_snowf = (snow_depth - temp_snow_depth) / dtime;

  } //end of do_capsnow construct

  // set frac_sno_eff variable
  if (ltype == istsoil || ltype == istcrop) {
     if (subgridflag ==1) { 
        frac_sno_eff = frac_sno;
     } else {
        frac_sno_eff = 1.;
     }
  } else {
     frac_sno_eff = 1.;
  }

  if (ltype==istwet && t_grnd>tfrz) {
     h2osno=0.;
     snow_depth=0.;
  }

//When the snow accumulation exceeds 10 mm, initialize snow layer
//Currently, the water temperature for the precipitation is simply set
//as the surface air temperature
  newnode = 0 ; //flag for when snow node will be initialized
        if (snow_level == 0 && qflx_snow_grnd_col > 0.00 && frac_sno*snow_depth >= 0.010) {
           newnode = 1;
           snow_level = -1;
           dz[0] = snow_depth ;                    //meter
           z[0] = -0.50*dz[0];
           zi[-1] = -dz[0];
           t_soisno[0] = min(tfrz, forc_air_temp) ;   //K
           h2osoi_ice[0] = h2osno ;            //kg/m2
           h2osoi_liq[0] = 0.0  ;               //kg/m2
           frac_iceold[0] = 1.0;
        }

//The change of ice partial density of surface node due to precipitation.
//Only ice part of snowfall is added here, the liquid part will be added
//later.
        if (snow_level < 0 && newnode == 0) {
        h2osoi_ice[snow_level+1] = h2osoi_ice[snow_level+1]+newsnow;
        dz[snow_level+1] = dz[snow_level+1]+dz_snowf*dtime;
        }
  }
 }

