#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <Kokkos_Core.hpp>
#include "landunit_varcon.h"    
#include "column_varcon.h" 
#include "clm_varpar.h"         
#include "clm_varctl.h"   
#include "utils.hh"
#include "readers.hh"

namespace ELM {
namespace Utils {

static const int n_months = 12;
static const int n_pfts = 17;
static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                            // day * half hour timestep
static const int n_grid_cells = 24;
static const int n_levels_snow = 5;

using MatrixStatePFT = MatrixStatic<n_grid_cells, n_pfts>;
using MatrixStateSoilColumn = MatrixStatic<n_grid_cells, n_levels_snow>;
using MatrixForc = MatrixStatic<n_max_times,n_grid_cells>;
using VectorColumn = VectorStatic<n_grid_cells>;
using VectorColumnInt = VectorStatic<n_grid_cells,int>;

} // namespace
} // namespace

namespace ELM {
KOKKOS_INLINE_FUNCTION void CanopyHydrology_Interception(double dtime,
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

namespace ELM {

KOKKOS_INLINE_FUNCTION void CanopyHydrology_FracWet(const int& frac_veg_nosno,
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


namespace ELM {

template<typename Array_d>
KOKKOS_INLINE_FUNCTION void CanopyHydrology_SnowWater(const double& dtime,
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
        bifall=50. + 1.7*pow((17.0),1.5);
     } else if (forc_air_temp > tfrz - 15.) {
        bifall=50. + 1.7*pow((forc_air_temp - tfrz + 15.),1.5);
     } else {
        bifall=50.;
     }

     // newsnow is all snow that doesn't fall on h2osfc
     newsnow = (1. - frac_h2osfc) * qflx_snow_grnd_col * dtime;

     // update integrated_snow
     integrated_snow = max(integrated_snow,h2osno) ; //h2osno could be larger due to frost

     // snowmelt from previous time step * dtime
     snowmelt = qflx_snow_melt * dtime;

     // set shape factor for accumulation of snow
     accum_factor=0.1;

     if (h2osno > 0.0) {

        //======================  FSCA PARAMETERIZATIONS  ======================
        // fsca parameterization based on *changes* in swe
        // first compute change from melt during previous time step
        if(snowmelt > 0.) {

           smr=min(1.,(h2osno)/(integrated_snow));

           frac_sno = 1. - pow((acos(fmin(1.,(2.*smr - 1.)))/rpi),(n_melt)) ;

        }

        // update fsca by new snow event, add to previous fsca
        if (newsnow > 0.) {
           fsno_new = 1. - (1. - tanh(accum_factor*newsnow))*(1. - frac_sno);
           frac_sno = fsno_new;

           // reset integrated_snow after accumulation events
           temp_intsnow= (h2osno + newsnow) / (0.5*(cos(rpi*pow((1.0-max(frac_sno,1e-6)),(1.0/n_melt))+1.0))) ;
           integrated_snow = min(1.e8,temp_intsnow) ;
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
              frac_sno = tanh(snow_depth/(2.5*zlnd*pow((min(800.0,(h2osno+ newsnow)/snow_depth)/100.0),1.0)) ) ;
           }
           if(h2osno < 1.0)  {
              frac_sno=min(frac_sno,h2osno);
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
           temp_intsnow= (h2osno + newsnow) / (0.5*(cos(rpi*pow((1.0-max(frac_sno,1e-6)),(1.0/n_melt)))+1.0));
           integrated_snow = min(1.e8,temp_intsnow);

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
                 frac_sno = tanh(snow_depth/(2.5*zlnd*pow((min(800.0,newsnow/snow_depth)/100.0),1.0)) );
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
           t_soisno[0] = fmin(tfrz, forc_air_temp) ;   //K
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

namespace ELM {

KOKKOS_INLINE_FUNCTION void CanopyHydrology_FracH2OSfc(const double& dtime,
        const double& min_h2osfc,
        const int& ltype,
        const double& micro_sigma,
        double& h2osno,
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




int main(int argc, char ** argv)
{
  using ELM::Utils::n_months;
  using ELM::Utils::n_pfts;
  using ELM::Utils::n_grid_cells;
  using ELM::Utils::n_max_times;
  using ELM::Utils::n_levels_snow;
  
  // fixed magic parameters for now
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  int n_irrig_steps_left = 0;

  const double dewmx = 0.1;
  const double dtime = 1800.0;

  // fixed magic parameters for SnowWater
  const double qflx_snow_melt = 0.;

  // fixed magic parameters for fracH2Osfc  
  const int oldfflag = 0;
  const double micro_sigma = 0.1;
  const double min_h2osfc = 1.0e-8;
  const double n_melt = 0.7;
  double qflx_floodg = 0.0;
  
  Kokkos::initialize( argc, argv );
  {                             
  // phenology input
  typedef Kokkos::View<double*>   ViewVectorType;
  typedef Kokkos::View<double**>  ViewMatrixType;
  typedef Kokkos::View<int**>  ViewMatrixType1;
  typedef Kokkos::View<int*>   ViewVectorType1;
  // ELM::Utils::MatrixState elai;
  // ELM::Utils::MatrixState esai;
  ViewMatrixType elai( "elai", n_grid_cells, n_pfts );
  ViewMatrixType esai( "esai", n_grid_cells, n_pfts );
  ViewMatrixType::HostMirror h_elai = Kokkos::create_mirror_view( elai );
  ViewMatrixType::HostMirror h_esai = Kokkos::create_mirror_view( esai );
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, h_elai, h_esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, h_elai, h_esai);

  // forcing input
  ViewMatrixType forc_rain( "forc_rain", n_max_times,n_grid_cells );
  ViewMatrixType forc_snow( "forc_snow", n_max_times,n_grid_cells );
  ViewMatrixType forc_air_temp( "forc_air_temp", n_max_times,n_grid_cells );
  ViewMatrixType::HostMirror h_forc_rain = Kokkos::create_mirror_view( forc_rain );
  ViewMatrixType::HostMirror h_forc_snow = Kokkos::create_mirror_view( forc_snow );
  ViewMatrixType::HostMirror h_forc_air_temp = Kokkos::create_mirror_view( forc_air_temp );
  const int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells, h_forc_rain, h_forc_snow, h_forc_air_temp);
  // ELM::Utils::MatrixForc forc_irrig; forc_irrig = 0.;
  ViewMatrixType forc_irrig( "forc_irrig", n_max_times,n_grid_cells );
  ViewMatrixType::HostMirror h_forc_irrig = Kokkos::create_mirror_view( forc_irrig );
  

  
  // mesh input (though can also change as snow layers evolve)
  //
  // NOTE: in a real case, these would be populated, but we don't actually
  // // need them to be for these kernels. --etc
  // auto z = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto zi = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto dz = ELM::Utils::MatrixStateSoilColumn(0.);
  ViewMatrixType z( "z", n_grid_cells, n_levels_snow );
  ViewMatrixType zi( "zi", n_grid_cells, n_levels_snow );
  ViewMatrixType dz( "dz", n_grid_cells, n_levels_snow );
  ViewMatrixType::HostMirror h_z = Kokkos::create_mirror_view( z );
  ViewMatrixType::HostMirror h_zi = Kokkos::create_mirror_view( zi );
  ViewMatrixType::HostMirror h_dz = Kokkos::create_mirror_view( dz );

  // state variables that require ICs and evolve (in/out)
  // auto h2ocan = ELM::Utils::MatrixStatePFT(0.);
  // auto swe_old = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto h2osoi_liq = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto h2osoi_ice = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto t_soisno = ELM::Utils::MatrixStateSoilColumn(0.);
  // auto frac_iceold = ELM::Utils::MatrixStateSoilColumn(0.);
  ViewMatrixType h2ocan( "h2ocan", n_grid_cells, n_pfts );
  ViewMatrixType swe_old( "swe_old", n_grid_cells, n_levels_snow );
  ViewMatrixType h2osoi_liq( "h2osoi_liq", n_grid_cells, n_levels_snow );
  ViewMatrixType h2osoi_ice( "h2osoi_ice", n_grid_cells, n_levels_snow );
  ViewMatrixType t_soisno( "t_soisno", n_grid_cells, n_levels_snow );
  ViewMatrixType frac_iceold( "frac_iceold", n_grid_cells, n_levels_snow );
  ViewMatrixType::HostMirror h_h2ocan = Kokkos::create_mirror_view( h2ocan );
  ViewMatrixType::HostMirror h_swe_old = Kokkos::create_mirror_view( swe_old );
  ViewMatrixType::HostMirror h_h2osoi_liq = Kokkos::create_mirror_view( h2osoi_liq );
  ViewMatrixType::HostMirror h_h2osoi_ice = Kokkos::create_mirror_view( h2osoi_ice );
  ViewMatrixType::HostMirror h_t_soisno = Kokkos::create_mirror_view( t_soisno );
  ViewMatrixType::HostMirror h_frac_iceold = Kokkos::create_mirror_view( frac_iceold );

  // auto t_grnd = ELM::Utils::VectorColumn(0.);
  // auto h2osno = ELM::Utils::VectorColumn(0.);
  // auto snow_depth = ELM::Utils::VectorColumn(0.);
  // auto snl = ELM::Utils::VectorColumnInt(0.); // note this tracks the snow_depth
  ViewVectorType t_grnd( "t_grnd", n_grid_cells );
  ViewVectorType h2osno( "h2osno", n_grid_cells );
  ViewVectorType snow_depth( "snow_depth", n_grid_cells );
  ViewVectorType1 snow_level( "snow_level", n_grid_cells );
  ViewVectorType::HostMirror h_t_grnd = Kokkos::create_mirror_view(  t_grnd);
  ViewVectorType::HostMirror h_h2osno = Kokkos::create_mirror_view( h2osno);
  ViewVectorType::HostMirror h_snow_depth = Kokkos::create_mirror_view(  snow_depth);
  ViewVectorType1::HostMirror h_snow_level = Kokkos::create_mirror_view( snow_level);

  // auto h2osfc = ELM::Utils::VectorColumn(0.);
  // auto frac_h2osfc = ELM::Utils::VectorColumn(0.);
  ViewVectorType h2osfc( "h2osfc", n_grid_cells );
  ViewVectorType frac_h2osfc( "frac_h2osfc", n_grid_cells );
  ViewVectorType::HostMirror h_h2osfc = Kokkos::create_mirror_view(  h2osfc);
  ViewVectorType::HostMirror h_frac_h2osfc = Kokkos::create_mirror_view( frac_h2osfc);

  
  // output fluxes by pft
  // auto qflx_prec_intr = ELM::Utils::MatrixStatePFT();
  // auto qflx_irrig = ELM::Utils::MatrixStatePFT();
  // auto qflx_prec_grnd = ELM::Utils::MatrixStatePFT();
  // auto qflx_snwcp_liq = ELM::Utils::MatrixStatePFT();
  // auto qflx_snwcp_ice = ELM::Utils::MatrixStatePFT();
  // auto qflx_snow_grnd_patch = ELM::Utils::MatrixStatePFT();
  // auto qflx_rain_grnd = ELM::Utils::MatrixStatePFT();
  ViewMatrixType qflx_prec_intr( "qflx_prec_intr", n_grid_cells, n_pfts );
  ViewMatrixType qflx_irrig( "qflx_irrig", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_prec_grnd( "qflx_prec_grnd", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_snwcp_liq( "qflx_snwcp_liq", n_grid_cells, n_pfts );
  ViewMatrixType qflx_snwcp_ice ( "qflx_snwcp_ice ", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_snow_grnd_patch( "qflx_snow_grnd_patch", n_grid_cells, n_pfts  );
  ViewMatrixType qflx_rain_grnd( "qflx_rain_grnd", n_grid_cells, n_pfts  );
  ViewMatrixType::HostMirror h_qflx_prec_intr = Kokkos::create_mirror_view( qflx_prec_intr );
  ViewMatrixType::HostMirror h_qflx_irrig = Kokkos::create_mirror_view( qflx_irrig);
  ViewMatrixType::HostMirror h_qflx_prec_grnd = Kokkos::create_mirror_view( qflx_prec_grnd);
  ViewMatrixType::HostMirror h_qflx_snwcp_liq = Kokkos::create_mirror_view(  qflx_snwcp_liq);
  ViewMatrixType::HostMirror h_qflx_snwcp_ice = Kokkos::create_mirror_view( qflx_snwcp_ice   );
  ViewMatrixType::HostMirror h_qflx_snow_grnd_patch = Kokkos::create_mirror_view( qflx_snow_grnd_patch );
  ViewMatrixType::HostMirror h_qflx_rain_grnd = Kokkos::create_mirror_view(  qflx_rain_grnd  );

  // FIXME: I have no clue what this is... it is inout on WaterSnow.  For now I
  // am guessing the data structure. Ask Scott.  --etc
  //auto integrated_snow = ELM::Utils::VectorColumn(0.);
  ViewVectorType integrated_snow( "integrated_snow", n_grid_cells );
  ViewVectorType::HostMirror h_integrated_snow = Kokkos::create_mirror_view(  integrated_snow);
  
  // output fluxes, state by the column
  // auto qflx_snow_grnd_col = ELM::Utils::VectorColumn();
  // auto qflx_snow_h2osfc = ELM::Utils::VectorColumn();
  // auto qflx_h2osfc2topsoi = ELM::Utils::VectorColumn();
  // auto qflx_floodc = ELM::Utils::VectorColumn();
  ViewVectorType qflx_snow_grnd_col( "qflx_snow_grnd_col", n_grid_cells );
  ViewVectorType qflx_snow_h2osfc( "qflx_snow_h2osfc", n_grid_cells );
  ViewVectorType qflx_h2osfc2topsoi( "qflx_h2osfc2topsoi", n_grid_cells );
  ViewVectorType qflx_floodc( "qflx_floodc", n_grid_cells );
  ViewVectorType::HostMirror h_qflx_snow_grnd_col = Kokkos::create_mirror_view(  qflx_snow_grnd_col);
  ViewVectorType::HostMirror h_qflx_snow_h2osfc = Kokkos::create_mirror_view( qflx_snow_h2osfc);
  ViewVectorType::HostMirror h_qflx_h2osfc2topsoi = Kokkos::create_mirror_view(  qflx_h2osfc2topsoi);
  ViewVectorType::HostMirror h_qflx_floodc = Kokkos::create_mirror_view( qflx_floodc);

  // auto frac_sno_eff = ELM::Utils::VectorColumn();
  // auto frac_sno = ELM::Utils::VectorColumn();
  ViewVectorType frac_sno_eff( "frac_sno_eff", n_grid_cells );
  ViewVectorType frac_sno( "frac_sno", n_grid_cells );
  ViewVectorType::HostMirror h_frac_sno_eff = Kokkos::create_mirror_view(  frac_sno_eff);
  ViewVectorType::HostMirror h_frac_sno = Kokkos::create_mirror_view( frac_sno);
  

  // std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  // auto min_max = std::minmax_element(&h_h2ocan(0,0), end1+1);
  // std::cout << std::minmax_element(16)
  //           << 0 << "\t" << std::accumulate(&h_h2ocan(0,0), end1+1, 0.)
  //           << "\t" << *min_max.first
  //           << "\t" << *min_max.second << std::endl;

    
  
  Kokkos::deep_copy( elai, h_elai);
  Kokkos::deep_copy( esai, h_esai);
  Kokkos::deep_copy( forc_rain, h_forc_rain);
  Kokkos::deep_copy( forc_snow, h_forc_snow);
  Kokkos::deep_copy( forc_irrig, h_forc_irrig);
  Kokkos::deep_copy( forc_air_temp, h_forc_air_temp);
  Kokkos::deep_copy( z, h_z);
  Kokkos::deep_copy( zi, h_zi);
  Kokkos::deep_copy( dz, h_dz);
  Kokkos::deep_copy( h2ocan, h_h2ocan);
  Kokkos::deep_copy( swe_old, h_swe_old);
  Kokkos::deep_copy( h2osoi_liq, h_h2osoi_liq);
  Kokkos::deep_copy( h2osoi_ice, h_h2osoi_ice);
  Kokkos::deep_copy( t_soisno, h_t_soisno);
  Kokkos::deep_copy( frac_iceold, h_frac_iceold);
  Kokkos::deep_copy( t_grnd, h_t_grnd);
  Kokkos::deep_copy( h2osno, h_h2osno);
  Kokkos::deep_copy( snow_depth, h_snow_depth);
  Kokkos::deep_copy( snow_level, h_snow_level);
  Kokkos::deep_copy( h2osfc, h_h2osfc);
  Kokkos::deep_copy( frac_h2osfc, h_frac_h2osfc);
  Kokkos::deep_copy( qflx_prec_intr,h_qflx_prec_intr);
  Kokkos::deep_copy( qflx_irrig,h_qflx_irrig);
  Kokkos::deep_copy( qflx_prec_grnd,h_qflx_prec_grnd);
  Kokkos::deep_copy( qflx_snwcp_liq,h_qflx_snwcp_liq);
  Kokkos::deep_copy( qflx_snwcp_ice,h_qflx_snwcp_ice);
  Kokkos::deep_copy( qflx_snow_grnd_patch,h_qflx_snow_grnd_patch);
  Kokkos::deep_copy( qflx_rain_grnd,h_qflx_rain_grnd);
  Kokkos::deep_copy( integrated_snow,h_integrated_snow);
  Kokkos::deep_copy( qflx_snow_grnd_col, h_qflx_snow_grnd_col);
  Kokkos::deep_copy( qflx_snow_h2osfc, h_qflx_snow_h2osfc);
  Kokkos::deep_copy( qflx_h2osfc2topsoi, h_qflx_h2osfc2topsoi);
  Kokkos::deep_copy( qflx_floodc, h_qflx_floodc);
  Kokkos::deep_copy( frac_sno_eff, h_frac_sno_eff);
  Kokkos::deep_copy( frac_sno, h_frac_sno);

  double* end1 = &h_h2ocan(n_grid_cells-1, n_pfts-1) ;
  double* end2 = &h_h2osno(n_grid_cells-1) ;
  double* end3 = &h_frac_h2osfc(n_grid_cells-1) ;
  std::cout << "Time\t Total Canopy Water\t Min Water\t Max Water\t Total Snow\t Min Snow\t Max Snow\t Avg Frac Sfc\t Min Frac Sfc\t Max Frac Sfc" << std::endl;
  auto min_max_water = std::minmax_element(&h_h2ocan(0,0), end1+1);
  auto sum_water = std::accumulate(&h_h2ocan(0,0), end1+1, 0.);

  auto min_max_snow = std::minmax_element(&h_h2osno(0), end2+1);
  auto sum_snow = std::accumulate(&h_h2osno(0), end2+1, 0.);

  auto min_max_frac_sfc = std::minmax_element(&h_frac_h2osfc(0), end3+1);
  auto avg_frac_sfc = std::accumulate(&h_frac_h2osfc(0), end3+1, 0.) / (end3+1 - &h_frac_h2osfc(0));

  std::cout << std::setprecision(16)
            << 0 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
            << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
            << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;


  // main loop
  // -- the timestep loop cannot/should not be parallelized
  for (size_t t = 0; t != n_times; ++t) {

    // grid cell and/or pft loop can be parallelized
    //for (size_t g = 0; g != n_grid_cells; ++g) {

      // PFT level operations
      //for (size_t p = 0; p != n_pfts; ++p) {
    Kokkos::parallel_for("n_grid_cells", n_grid_cells, KOKKOS_LAMBDA (const size_t& g) {
      for (size_t p = 0; p != n_pfts; ++p) {
        //
        // Calculate interception
        //
        // NOTE: this currently punts on what to do with the qflx variables!
        // Surely they should be either std::accumulated or stored on PFTs as well.
        // --etc
        ELM::CanopyHydrology_Interception(dtime,
                forc_rain(t,g), forc_snow(t,g), forc_irrig(t,g),
                ltype, ctype, urbpoi, do_capsnow,
                elai(g,p), esai(g,p), dewmx, frac_veg_nosno,
                h2ocan(g,p), n_irrig_steps_left,
                qflx_prec_intr(g,p), qflx_irrig(g,p), qflx_prec_grnd(g,p),
                qflx_snwcp_liq(g,p), qflx_snwcp_ice(g,p),
                qflx_snow_grnd_patch(g,p), qflx_rain_grnd(g,p));
        //printf("%i %i %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", g, p, forc_rain(t,g), forc_snow(t,g), elai(g,p), esai(g,p), h2ocan(g,p), qflx_prec_intr(g));

        //
        // Calculate fraction of LAI that is wet vs dry.
        //
        // FIXME: this currently punts on what to do with the fwet/fdry variables.
        // Surely they should be something, as such this is dead code.
        // By the PFT?
        // --etc
        double fwet = 0., fdry = 0.;
        ELM::CanopyHydrology_FracWet(frac_veg_nosno, h2ocan(g,p), elai(g,p), esai(g,p), dewmx, fwet, fdry);
      } // end PFT loop

      // Column level operations
      
      double* qpatch = &qflx_snow_grnd_patch(n_grid_cells-1, n_pfts-1);
      // NOTE: this is effectively an accumulation kernel/task! --etc
      qflx_snow_grnd_col(g) = std::accumulate(&qflx_snow_grnd_patch(0,0), qpatch+1, 0.);

      // Calculate ?water balance? on the snow column, adding throughfall,
      // removing melt, etc.
      //
      // local outputs
      int newnode;

      ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
              ltype, ctype, urbpoi, do_capsnow, oldfflag,
              forc_air_temp(t,g), t_grnd(g),
              qflx_snow_grnd_col(g), qflx_snow_melt, n_melt, frac_h2osfc(g),
              snow_depth(g), h2osno(g), integrated_snow(g), Kokkos::subview(swe_old, g , Kokkos::ALL),
              Kokkos::subview(h2osoi_liq, g , Kokkos::ALL), Kokkos::subview(h2osoi_ice, g , Kokkos::ALL), Kokkos::subview(t_soisno, g , Kokkos::ALL), Kokkos::subview(frac_iceold, g , Kokkos::ALL),
              snow_level(g), Kokkos::subview(dz, g , Kokkos::ALL), Kokkos::subview(z, g , Kokkos::ALL), Kokkos::subview(zi, g , Kokkos::ALL), newnode,
              qflx_floodc(g), qflx_snow_h2osfc(g), frac_sno_eff(g), frac_sno(g));

      // Calculate Fraction of Water to the Surface?
      //
      // FIXME: Fortran black magic... h2osoi_liq is a vector, but the
      // interface specifies a single double.  For now passing the 0th
      // entry. --etc
      ELM::CanopyHydrology_FracH2OSfc(dtime, min_h2osfc, ltype, micro_sigma,
              h2osno(g), h2osfc(g), h2osoi_liq(g,0), frac_sno(g), frac_sno_eff(g),
              qflx_h2osfc2topsoi(g), frac_h2osfc(g));
      
    }); // end grid cell loop

    
    // auto min_max = std::minmax_element(&h_h2ocan(0,0), end1+1);
    // std::cout << std::minmax_element(16)
    //           << t+1 << "\t" << std::accumulate(&h_h2ocan(0,0), end1+1, 0.)
    //           << "\t" << *min_max.first
    //           << "\t" << *min_max.second << std::endl;

    auto min_max_water = std::minmax_element(&h_h2ocan(0,0), end1+1);
    auto sum_water = std::accumulate(&h_h2ocan(0,0), end1+1, 0.);

    auto min_max_snow = std::minmax_element(&h_h2osno(0), end2+1);
    auto sum_snow = std::accumulate(&h_h2osno(0), end2+1, 0.);

    auto min_max_frac_sfc = std::minmax_element(&h_frac_h2osfc(0), end3+1);
    auto avg_frac_sfc = std::accumulate(&h_frac_h2osfc(0), end3+1, 0.) / (end3+1 - &h_frac_h2osfc(0));
                  
    std::cout << std::setprecision(16)
              << 0 << "\t" << sum_water << "\t" << *min_max_water.first << "\t" << *min_max_water.second
              << "\t" << sum_snow << "\t" << *min_max_snow.first << "\t" << *min_max_snow.second
              << "\t" << avg_frac_sfc << "\t" << *min_max_frac_sfc.first << "\t" << *min_max_frac_sfc.second << std::endl;

  } // end timestep loop
  }
  Kokkos::finalize();
  return 0;
}
