
#pragma once
// InitCold() from SoilStateType.F90
// sets watsat, sucsat, bsw, etc.

// also contains code from ColumnDataType.F90 col_ws_init, col_es_init
// should maybe also include some reads from col_*_restart -- later

//#include "RootBioPhys.hh"
#include "array.hh"
#include "elm_constants.h"
#include "land_data.h"

namespace ELM {

//void InitSoil(const int vtype, ArrayD1 rootfr) {
//
//  //smpmin = -1.0e8; global for now - maybe change
//
//  for (int i = 0; i < ELM::nlevgrnd; ++i) {
//    rootfr(i) = 0.0;
//  }
//
//  init_vegrootfr(Land.vtype, rootfr);
//}



void pedotransfer(const double& pct_sand, const double& pct_clay, double& watsat, double& bsw, double& sucsat, double& xksat) {
  //compute hydraulic properties based on functions derived from Table 5 in cosby et al, 1984
   
  //Cosby et al. Table 5     
  watsat = 0.489 - 0.00126 * pct_sand;
  bsw = 2.91 + 0.159 * pct_clay;
  sucsat = 10.0 * pow(10.0, (1.88 - 0.0131 * pct_sand));
  xksat = 0.0070556 * pow (10.0, (-0.884 + 0.0153 * pct_sand)); // mm/s, from table 5 
}


void soil_hydraulic_params(const double& pct_sand, const double& pct_clay, const double& zsoi, const double& om_frac, 
  double& watsat, double& bsw, double& sucsat, double& watdry, double& watopt, double& watfc) {

static const double zsapric = 0.5; // depth (m) that organic matter takes on characteristics of sapric peat
static const double pcalpha = 0.5; // percolation threshold
static const double pcbeta  = 0.139; // percolation exponent

double xksat;
pedotransfer(pct_sand, pct_clay, watsat, bsw, sucsat, xksat);
const double om_watsat = std::max(0.93 - 0.1 * (zsoi / zsapric), 0.83);
const double om_b = std::min(2.7 + 9.3 * (zsoi / zsapric), 12.0);
const double om_sucsat = std::min(10.3 - 0.2 * (zsoi / zsapric), 10.1);
const double om_hksat = std::max(0.28 - 0.2799 * (zsoi /zsapric), 0.0001);

//const double bulk_den = (1.0 - watsat) * 2.7e3;
//const double tkm = (1.0 - om_frac) * (8.8 * sand + 2.92 * clay) / (sand + clay) + om_tkm * om_frac; // W/(m K)
watsat = (1.0 - om_frac) * watsat + om_watsat * om_frac;
bsw = (1.0 - om_frac) * (2.91 + 0.159 * pct_clay) + om_frac * om_b;
sucsat = (1.0 - om_frac) * sucsat + om_sucsat * om_frac;
//hksat_min(i) = xksat;

// perc_frac is zero unless perf_frac greater than percolation threshold
double perc_frac;
if (om_frac > pcalpha) {
  double perc_norm = pow((1.0 - pcalpha), -pcbeta);
  perc_frac = perc_norm * pow((om_frac - pcalpha), pcbeta);
} else {
  perc_frac = 0.0;
}

// uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
double uncon_frac = (1.0 - om_frac) + (1.0 - perc_frac) * om_frac;

// uncon_hksat is series addition of mineral/organic conductivites
double uncon_hksat;
if (om_frac < 1.0) {
   uncon_hksat = uncon_frac / ((1.0 - om_frac) / xksat
        + ((1.0 - perc_frac) * om_frac) / om_hksat);
} else {
   uncon_hksat = 0.0;
}

double hksat  = uncon_frac * uncon_hksat + (perc_frac * om_frac) * om_hksat;

//this%tkmg_col(c,lev)   = tkm ** (1._r8- this%watsat_col(c,lev))           

//this%tksatu_col(c,lev) = this%tkmg_col(c,lev)*0.57_r8**this%watsat_col(c,lev)

//this%tkdry_col(c,lev)  = ((0.135_r8*this%bd_col(c,lev) + 64.7_r8) / &
//     (2.7e3_r8 - 0.947_r8*this%bd_col(c,lev)))*(1._r8-om_frac) + om_tkd*om_frac  

//this%csol_col(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) + &
//     om_csol*om_frac)*1.e6_r8  ! J/(m3 K)

//if (lev > nlevbed) then
//   this%csol_col(c,lev) = csol_bedrock
//endif

watdry = watsat *
     pow((316230.0 / sucsat), (-1.0 / bsw));
watopt = watsat *
     pow((158490.0 / sucsat), (-1.0 / bsw));



// added by K.Sakaguchi for beta from Lee and Pielke, 1992
// water content at field capacity, defined as hk = 0.1 mm/day
// used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / secspday (day/sec)
watfc = watsat *
     pow((0.1 / (hksat * ELM::secspday)), (1.0 / (2.0 * bsw + 3.0)));

//this%sucmin_col(c,lev) = min_liquid_pressure

//this%watmin_col(c,lev) = &
//     this%watsat_col(c,lev)*(-min_liquid_pressure/this%sucsat_col(c,lev))**(-1._r8/this%bsw_col(c,lev))
}


template <typename ArrayD1>
void init_soil_hydraulics(const ArrayD1& pct_sand, const ArrayD1& pct_clay, const ArrayD1& organic, const ArrayD1& zsoi, 
  ArrayD1 watsat, ArrayD1 bsw, ArrayD1 sucsat, ArrayD1 watdry, ArrayD1 watopt, ArrayD1 watfc) {

double om_frac;
for (int i = 0; i < ELM::nlevsoi; ++i) {
  om_frac = pow((organic(i) / ELM::organic_max), 2.0);
  soil_hydraulic_params(pct_sand(i), pct_clay(i), zsoi(i+ELM::nlevsno), om_frac,
  watsat(i), bsw(i), sucsat(i), watdry(i), watopt(i), watfc(i));
}


for (int i = ELM::nlevsoi; i < ELM::nlevgrnd; ++i) {
  om_frac = 0.0;
  soil_hydraulic_params(pct_sand(ELM::nlevsoi-1), pct_clay(ELM::nlevsoi-1), zsoi(i+ELM::nlevsno), om_frac,
  watsat(i), bsw(i), sucsat(i), watdry(i), watopt(i), watfc(i));
}

}



// from ColumnDataType.F90
//-----------------------------------------------------------------------
// set cold-start initial values for select members of col_es
//-----------------------------------------------------------------------
template <typename ArrayD1>
void init_soil_temp(const LandType& Land, const int& snl, ArrayD1 t_soisno, double& t_grnd) {

  // Snow level temperatures - all land points
  if (snl > 0) {
    for (int i = ELM::nlevsno-snl; i < ELM::nlevsno; ++i) { t_soisno(i) = 250.0; }
  }

  // Below snow temperatures - nonlake points (lake points are set below)
  if (!Land.lakpoi) {
    if (Land.ltype == istice || Land.ltype == istice_mec) {
      for (int i = ELM::nlevsno; i < ELM::nlevgrnd+ELM::nlevsno; ++i) { t_soisno(i) = 250.0; }
    } else if (Land.ltype == istwet) {
        for (int i = ELM::nlevsno; i < ELM::nlevgrnd+ELM::nlevsno; ++i) { t_soisno(i) = 277.0; }
    } else if (Land.urbpoi) {
      if (Land.ctype == icol_road_perv || Land.ctype == icol_road_imperv) {
        for (int i = ELM::nlevsno; i < ELM::nlevgrnd+ELM::nlevsno; ++i) { t_soisno(i) = 274.0; }
      } else if (Land.ctype == icol_sunwall || Land.ctype == icol_shadewall || Land.ctype == icol_roof) {
        // Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
        // shock from large heating/air conditioning flux
        for (int i = ELM::nlevsno; i < ELM::nlevurb+ELM::nlevsno; ++i) { t_soisno(i) = 292.0; }
      }
    } else {
      for (int i = ELM::nlevsno; i < ELM::nlevgrnd+ELM::nlevsno; ++i) { t_soisno(i) = 274.0; }
    }
    t_grnd = t_soisno(ELM::nlevsno - snl);
  }
}


// from ColumnDataType.F90 and WaterStateType.F90
//-----------------------------------------------------------------------
// set cold-start initial values for select members of col_ws
//-----------------------------------------------------------------------
template <typename ArrayD1>
void init_snow_state(const bool& urbpoi, const int& snl, double& h2osno, double& int_snow, double& snow_depth, double& h2osfc, double& h2ocan, double& frac_h2osfc, double& fwet, double& fdry, double& frac_sno, ArrayD1 snw_rds) {

  // this could be intitialized from input data - 0.0 for now
  h2osno = 0.0;
  int_snow = 0.0;
  snow_depth  = 0.0;
  h2osfc = 0.0;
  h2ocan = 0.0;
  frac_h2osfc = 0.0;
  fwet = 0.0;
  fdry = 0.0;

  // initial snow fraction
  if (urbpoi) {
    // From Bonan 1996 (LSM technical note)
    frac_sno = std::min(snow_depth / 0.05, 1.0);
  } else {
    frac_sno = 0.0;
    // snow cover fraction as in Niu and Yang 2007
    if(snow_depth > 0.0) {
      const double snowbd = std::min(400.0, h2osno / snow_depth); // bulk density of snow (kg/m3)
      const double fmelt = pow(snowbd/100.0, 1.0);
      // 100 is the assumed fresh snow density; 1 is a melting factor that could be
      // reconsidered, optimal value of 1.5 in Niu et al., 2007
      frac_sno = tanh(snow_depth / (2.5 * zlnd * fmelt));
    }
  }
  
  // initial snow radius
  if (snl > 0) {
    for (int i = 0; i < ELM::nlevsno-snl; ++i) { snw_rds(i) = 0.0; }
    for (int i = ELM::nlevsno-snl; i < ELM::nlevsno; ++i) { snw_rds(i) = snw_rds_min; }
  } else if (h2osno > 0.0) {
    snw_rds(ELM::nlevsno-1) = snw_rds_min;
    for (int i = 0; i < ELM::nlevsno-1; ++i) { snw_rds(i) = 0.0; }
  } else {
    for (int i = 0; i < ELM::nlevsno; ++i) { snw_rds(i) = 0.0; }
  }
}




      

// from ColumnDataType.F90 and WaterStateType.F90
//--------------------------------------------
// Set soil water
//--------------------------------------------
// volumetric water is set first and liquid content and ice lens are obtained
// NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
// and urban pervious road (other urban columns have zero soil water)
template <typename ArrayD1>
void init_soilh2o_state(const LandType& Land, const int& snl, const ArrayD1& watsat, const ArrayD1& t_soisno, const ArrayD1& dz, ArrayD1 h2osoi_vol, ArrayD1 h2osoi_liq, ArrayD1 h2osoi_ice) 

{
  for (int i = 0; i < ELM::nlevgrnd; ++i) { h2osoi_vol(i) = spval; }
  for (int i = 0; i < ELM::nlevgrnd+ELM::nlevsno; ++i) { h2osoi_liq(i) = spval; }
  for (int i = 0; i < ELM::nlevgrnd+ELM::nlevsno; ++i) { h2osoi_ice(i) = spval; }



  int nlevs = ELM::nlevgrnd;
  if (!Land.lakpoi) {
  
    if (Land.ltype == istsoil || Land.ltype == istcrop) {
      for (int i = 0; i < ELM::nlevgrnd; ++i) {
        if (i >= ELM::nlevbed) {
          h2osoi_vol(i) = 0.0;
        } else {
          //if (use_fates_planthydro .or. use_hydrstress) {
          //  h2osoi_vol(i) = 0.70_r8*watsat(i) !0.15_r8 to avoid very dry conditions that cause errors in FATES HYDRO
          //} else {
          h2osoi_vol(i) = 0.15;
        }
      } 
    } else if (Land.urbpoi) {
      if (Land.ctype == icol_road_perv) {
        for (int i = 0; i < ELM::nlevgrnd; ++i) {
          if (i < ELM::nlevbed) {
            h2osoi_vol(i) = 0.3;
          } else {
            h2osoi_vol(i) = 0.0;
          }
        }
      } else if (Land.ctype == icol_road_imperv) {
        for (int i = 0; i < ELM::nlevgrnd; ++i) {
          h2osoi_vol(i) = 0.0;
        }
      } else {
        nlevs = ELM::nlevurb;
        for (int i = 0; i < ELM::nlevurb; ++i) {
          h2osoi_vol(i) = 0.0;
        }
      }
    } else if (Land.ltype == istwet) {
      for (int i = 0; i < ELM::nlevgrnd; ++i) {
        if (i >= ELM::nlevbed) {
          h2osoi_vol(i) = 0.0;
        } else {
          h2osoi_vol(i) = 1.0;
        }
      }
    } else if (Land.ltype == istice || Land.ltype == istice_mec) {
      for (int i = 0; i < ELM::nlevgrnd; ++i) {
        h2osoi_vol(i) = 1.0;
      }
    }
    for (int i = 0; i < nlevs; ++i) {
      const int snw_offset = i+ELM::nlevsno;
      h2osoi_vol(i) = std::min(h2osoi_vol(i), watsat(i));
      if (t_soisno(snw_offset) <= ELM::tfrz) {
        h2osoi_ice(snw_offset) = dz(snw_offset) * ELM::denice * h2osoi_vol(i);
        h2osoi_liq(snw_offset) = 0.0;
      } else {
        h2osoi_ice(snw_offset) = 0.0;
        h2osoi_liq(snw_offset) = dz(snw_offset) * ELM::denh2o * h2osoi_vol(i);
      }
    }
  
    for (int i = 0; i < ELM::nlevsno; ++i) {
      if (i >= ELM::nlevsno-snl) {
        h2osoi_ice(i) = dz(i) * 250.0;
        h2osoi_liq(i) = 0.0;
      }
    }
  } else {
    //--------------------------------------------
    // Set Lake water
    //--------------------------------------------
    for (int i = 0; i < ELM::nlevsno; ++i) {
      if (i >= ELM::nlevsno-snl) {
        h2osoi_ice(i) = dz(i) * bdsno;
        h2osoi_liq(i) = 0.0;
      }
    }
    for (int i = 0; i < ELM::nlevgrnd; ++i) {
      const int snw_offset = i+ELM::nlevsno;
      if (i < ELM::nlevsoi) { // soil
         h2osoi_vol(i) = watsat(i);
         h2osoi_liq(snw_offset) = spval;
         h2osoi_ice(snw_offset) = spval;
      } else { // bedrock
         h2osoi_vol(i) = 0.0;
      }
    }
  }

  //--------------------------------------------
  // For frozen layers !TODO - does the following make sense ???? it seems to overwrite everything
  //--------------------------------------------
  for (int i = 0; i < ELM::nlevgrnd; ++i) {
    const int snw_offset = i+ELM::nlevsno;
    if (t_soisno(snw_offset) <= ELM::tfrz) {
      h2osoi_ice(snw_offset) = dz(snw_offset) * denice * h2osoi_vol(i);
      h2osoi_liq(snw_offset) = 0.0;
    } else {
      h2osoi_ice(snw_offset) = 0.0;
      h2osoi_liq(snw_offset) = dz(snw_offset) * denh2o * h2osoi_vol(i);
    }
  }
  //h2osoi_liq_old(c,:) = h2osoi_liq(c,:)
  //h2osoi_ice_old(c,:) = h2osoi_ice(c,:)
}
  



} // namespace ELM
