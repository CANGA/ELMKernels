
// h2ocan update - eq 7.9 in CLM 4.5 tech note - 
//h2ocan(n+1) = h2ocan(n) + intercepted_precip * dt - qflx_candrip * dt - evap * dt
// - evap * dt portion gets calculated in CanopyFluxes


#include "ELMConstants.h"
#include "LandType.h"
#include <algorithm>
#include <cmath>

namespace ELM {

// determine fraction of input precipitation that is snow and rain
double phase_precip_fraction(const double &total_precip, const double &phase_precip) {
  return total_precip > 0.0 ? phase_precip / total_precip : 0.0;
}

// set fraction of potential interception to max 0.25
double fraction_of_potential_interception(const double &elai, const double &esai) {
  return 0.25 * (1.0 - exp(-0.5 * (elai + esai)));
}

// Direct throughfall
double phase_throughfall(const double &phase_precip, const double &fpi) {
  return phase_precip > 0.0 ? phase_precip * (1.0 - fpi) : 0.0;
}


// maximum allowed water on canopy [mm]
double max_canopy_water(const double &dewmx, const double &elai, const double &esai) {
  return dewmx * (elai + esai);
}

// Intercepted precipitation [mm/s]
double intercepted_precip(const double total_precip, const double &fpi) {
  return total_precip * fpi;
}

// Excess water that exceeds the leaf capacity
double excess_canopy_water(const double &h2ocan, const double &h2ocanmx, const double &dtime) {
  return (h2ocan - h2ocanmx) / dtime;
}


// Canopy storage of interception and canopy drainage
void canopy_storage_and_drainage(const LandType &Land, const int &frac_veg_nosno, const double &total_precip, 
  const double &fpi, const double &dtime, const double &dewmx, const double &elai, const double &esai, double &h2ocan, 
  double &qflx_candrip) {
  if (!Land.lakpoi) {
    if (Land.ltype != istice && Land.ltype != istice_mec) {
      if (Land.ctype != icol_sunwall && Land.ctype != icol_shadewall) {
        if (frac_veg_nosno == 1 && total_precip > 0.0) {
          h2ocan += std::max(0.0, (dtime * intercepted_precip(total_precip, fpi) ));
          qflx_candrip = 0.0;
          double h2ocanmx = max_canopy_water(dewmx, elai, esai);
          double xrun = excess_canopy_water(h2ocan, h2ocanmx, dtime);
          if (xrun > 0.0) {
            qflx_candrip = xrun;
            h2ocan = h2ocanmx;
          }
        }
      }
    } else {
      h2ocan = 0.0;
      qflx_candrip = 0.0;
    }
  }
}

double phase_ground_precip(const LandType &Land, const int &frac_veg_nosno, const double &total_precip, 
  const double &phase_precip, const double &fpi, const double &qflx_candrip) {
  double phase_qflx_prec_grnd = 0.0;
  if ((Land.ctype != icol_sunwall) && (Land.ctype != icol_shadewall)) {
    if (frac_veg_nosno == 0) {
      phase_qflx_prec_grnd = phase_precip;
    } else {
      phase_qflx_prec_grnd = phase_throughfall(phase_precip, fpi) + (qflx_candrip * 
        phase_precip_fraction(total_precip, phase_precip));
    }
  } 
  return phase_qflx_prec_grnd;
}


void ground_precip_fluxes(const LandType &Land, const bool &do_capsnow, const int &frac_veg_nosno, 
  const double &total_precip, const double &forc_rain, const double &forc_snow, const double &fpi, 
  const double &qflx_candrip, const double &qflx_irrig,double &qflx_snwcp_liq, double &qflx_snwcp_ice, 
  double &qflx_snow_grnd, double &qflx_rain_grnd, double &qflx_prec_grnd) {
  if (!Land.lakpoi) {
    if (do_capsnow) {
      // Add irrigation water directly onto ground (bypassing canopy interception)
      qflx_snwcp_liq = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_rain, fpi, qflx_candrip) + qflx_irrig;
      qflx_snwcp_ice = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_snow, fpi, qflx_candrip);
      qflx_snow_grnd = 0.0;
      qflx_rain_grnd = 0.0;
      // Total water onto ground
      qflx_prec_grnd = qflx_snwcp_ice + qflx_snwcp_liq;
    } else {
      qflx_snwcp_liq = 0.0;
      qflx_snwcp_ice = 0.0;
      // Add irrigation water directly onto ground (bypassing canopy interception)
      qflx_rain_grnd = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_rain, fpi, qflx_candrip) + qflx_irrig;;
      qflx_snow_grnd = phase_ground_precip(Land, frac_veg_nosno, total_precip, forc_snow, fpi, qflx_candrip);
      // Total water onto ground
      qflx_prec_grnd = qflx_snow_grnd + qflx_rain_grnd;
    }
  }
}




void interception_physics(const LandType &Land, const bool &do_capsnow, const int &frac_veg_nosno, 
  const double &forc_rain, const double &forc_snow, const double &dtime, const double &dewmx, const double &elai, 
  const double &esai, double &h2ocan, const double &qflx_irrig, double &qflx_snwcp_liq, 
  double &qflx_snwcp_ice, double &qflx_snow_grnd, double &qflx_rain_grnd, double &qflx_prec_grnd) {

  double qflx_candrip;
  double total_precip = forc_snow + forc_rain;
  
  double fpi = fraction_of_potential_interception(elai, esai);
  
  canopy_storage_and_drainage(Land, frac_veg_nosno, total_precip, 
    fpi, dtime, dewmx, elai, esai, h2ocan, qflx_candrip);
  
  ground_precip_fluxes(Land, do_capsnow, frac_veg_nosno, total_precip, forc_rain, forc_snow, fpi, qflx_candrip, 
    qflx_irrig, qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd, qflx_rain_grnd, qflx_prec_grnd);

}


//in:
//Land
//qflx_snow_grnd
//dtime
//do_capsnow
//forc_t
//n_melt
//
//out:
//dz_snowf  - local
//newsnow  - local
//frac_sno - state
//frac_sno_eff - state
//int_snow - state
//h2osno - state
//snow_depth - state

void xxx(const LandType &Land, const bool &do_capsnow, const double &qflx_snow_grnd, const double &dtime, 
const double &forc_t, const double &n_melt, double &dz_snowf , double &newsnow , double &frac_sno, 
double &frac_sno_eff, double &int_snow, double &h2osno, double &snow_depth) {
  if (!Land.lakpoi) {
    if (do_capsnow) {
      dz_snowf = 0.0;
      newsnow = qflx_snow_grnd * dtime;
      frac_sno = 1.0;
      int_snow = 5.e2;
    } else {
      double bifall;
      if (forc_t > tfrz + 2.0) {
        bifall = 50.0 + 1.7 * pow(17.0, 1.5);
      } else if (forc_t > tfrz - 15.0) {
        bifall = 50.0 + 1.7 * pow((forc_t - tfrz + 15.0), 1.5);
      } else {
        bifall = 50.0;
      }
      // all snow falls on ground, no snow on h2osfc
      newsnow = qflx_snow_grnd * dtime;
      // update int_snow
      int_snow = std::max(int_snow, h2osno); // h2osno could be larger due to frost
      // snowmelt from previous time step * dtime
      double snowmelt = qflx_snow_melt * dtime;
      // set shape factor for accumulation of snow
      double accum_factor = 0.1;
      // save initial snow
      double temp_snow_depth = snow_depth;

      // calculate frac_sno and snow_depth
      update_sca_and_snow_depth(Land.urbpoi, oldfflag, h2osno, snowmelt,
      n_melt, accum_factor, newsnow, bifall, 
      frac_sno, snow_depth);

      h2osno = h2osno + newsnow; // update h2osno for new snow
      int_snow = int_snow + newsnow;
      dz_snowf = (snow_depth - temp_snow_depth) / dtime; // update change in snow depth
    } // end else do_capsnow
    if (Land.ltype == istwet && t_grnd > tfrz) {
      h2osno = 0.0;
      snow_depth = 0.0;
    }

    frac_sno_eff = effective_snow_fraction(Land.ltype, frac_sno);

  } // if !lakpoi
}


void save_initial_swe(const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, ArrayD1 swe_old) {
  // save initial snow content
  for (int j = 0; j < nlevsno - snl; j++) {
    swe_old[j] = 0.0;
  }
  for (int j = nlevsno - snl; j < nlevsno; j++) {
    swe_old[j] = h2osoi_liq[j] + h2osoi_ice[j];
  }
} 

//in:
//urbpoi
//h2osno
//snowmelt
//n_melt
//accum_factor
//newsnow
//bifall
//oldfflag
//
//out:
//frac_sno - state
//snow_depth - state


// SCA is based on CLM 4.5 tech note eqs 7.14 (accumulation) and 7.15 (depletion) 
void update_sca_and_snow_depth(const bool &urbpoi, const int &oldfflag, const double &h2osno, const double &snowmelt,
 const double &n_melt, const double &accum_factor, const double &newsnow, const double &bifall, 
 double &frac_sno, double &snow_depth) {
  /*======================  FSCA PARAMETERIZATIONS  ======================
  fsca parameterization based on *changes* in swe
  first compute change from melt during previous time step */
  if (h2osno > 0.0) {
    if (snowmelt > 0.0) {
      frac_sno = snow_area_depletion(h2osno, int_snow, n_melt);
    }
    // update fsca by new snow event, add to previous fsca
    if (newsnow > 0.0) {
      double fsno_new = 1.0 - (1.0 - tanh(accum_factor * newsnow)) * (1.0 - frac_sno);
      frac_sno = fsno_new;
      // reset int_snow after accumulation events
      int_snow = integrated_snow(h2osno, newsnow, frac_sno, n_melt);
    }
    /*====================================================================*/
    // for subgrid fluxes
    if (subgridflag == 1 && !urbpoi) {
      if (frac_sno > 0.0) {
        snow_depth = snow_depth + newsnow / (bifall * frac_sno);
      } else {
        snow_depth = 0.0;
      }
    } else {
      // for uniform snow cover
      snow_depth = snow_depth + newsnow / bifall;
    }
    if (oldfflag == 1) {
      // snow cover fraction in Niu et al. 2007
      if (snow_depth > 0.0) {
        frac_sno = snow_covered_area_fraction(h2osno, newsnow, snow_depth);
      }
      if (h2osno < 1.0) {
        frac_sno = std::min(frac_sno, h2osno);
      }
    }
  } else { // h2osno == 0
    // initialize frac_sno and snow_depth when no snow present initially
    if (newsnow > 0.0) {
      double z_avg = newsnow / bifall;
      frac_sno = tanh(accum_factor * newsnow);
      // make int_snow consistent w/ new fsno, h2osno
      // reset prior to adding newsnow below
      int_snow = integrated_snow(h2osno, newsnow, frac_sno, n_melt);
      // update snow_depth and h2osno to be consistent with frac_sno, z_avg
      if (subgridflag == 1 && !urbpoi) {
        snow_depth = z_avg / frac_sno;
      } else {
        snow_depth = newsnow / bifall;
      }
      // use n&y07 formulation
      if (oldfflag == 1) {
        // snow cover fraction in Niu et al. 2007
        if (snow_depth > 0.0) {
          frac_sno = snow_covered_area_fraction(h2osno, newsnow, snow_depth);
        }
      }
    } else {
      //z_avg = 0.0;
      snow_depth = 0.0;
      frac_sno = 0.0;
    }
  }
}



int initialize_new_snow_layer(const int &snl, const double &qflx_snow_grnd, const double &frac_sno, 
  const double &snow_depth, const double &snow_depth, const double &forc_t, const double &h2osno, 
  ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
  ArrayD1 frac_iceold, ArrayD1 snw_rds) {
  // When the snow accumulation exceeds 10 mm, initialize snow layer
  // Currently, the water temperature for the precipitation is simply set
  // as the surface air temperature
  int newnode = 0; // flag for when snow node will be initialized
  if (snl == 0 && qflx_snow_grnd > 0.0 && (frac_sno * snow_depth) >= 0.01) {
    newnode = 1;
    snl = 1;
    dz[nlevsno - snl] = snow_depth; // meter
    z[nlevsno - snl] = -0.5 * dz[nlevsno - snl];
    zi[nlevsno - snl] = -dz[nlevsno - snl];
    t_soisno[nlevsno - snl] = std::min(tfrz, forc_t); // K
    h2osoi_ice[nlevsno - snl] = h2osno;               // kg/m2
    h2osoi_liq[nlevsno - snl] = 0.0;                  // kg/m2
    frac_iceold[nlevsno - snl] = 1.0;
    snw_rds[nlevsno - snl] = snw_rds_min;
  }
  return newnode;
}


void apply_new_snow_to_top_layer(const int &snl, const int &newnode, const double &newsnow, const double &dz_snowf, 
  const double &dtime, ArrayD1 h2osoi_ice, ArrayD1 dz) {
  // The change of ice partial density of surface node due to precipitation.
  // Only ice part of snowfall is added here, the liquid part will be added later.
  if (snl > 0 && newnode == 0) {
    h2osoi_ice[nlevsno - snl] = h2osoi_ice[nlevsno - snl] + newsnow;
    dz[nlevsno - snl] = dz[nlevsno - snl] + dz_snowf * dtime;
  }
}

double snow_covered_area_fraction(const double &h2osno, const double &newsnow, const double &snow_depth) {
return tanh(snow_depth / (2.5 * zlnd *
              pow((std::min(800.0, ((h2osno + newsnow) / snow_depth / 100.0))), 1.0))); // why to the power of 1.0??
}

// I'm not sure what this means - it matches ELM output, but the values are really large - gets used in SnowWater()
double integrated_snow(const double &h2osno, const double &newsnow, const double &frac_sno, 
  const double &n_melt) {
  double temp_intsnow =
  (h2osno + newsnow) / (0.5 * (cos(ELM_PI * pow((1.0 - std::max(frac_sno, 1.e-6)), (1.0 / n_melt))) + 1.0));
  return std::min(1.e8, temp_intsnow);
}


double snow_area_depletion(const double &h2osno, const double &int_snow, const double &n_melt) {
  double smr = std::min(1.0, (h2osno / int_snow));
  return 1.0 - pow((acos(std::min(1.0, (2.0 * smr - 1.0))) / ELM_PI), n_melt);
}


double effective_snow_fraction(const int &ltype, const double &frac_sno) { // set frac_sno_eff variable
  if (Land.ltype == istsoil || Land.ltype == istcrop) {
    if (subgridflag == 1) {
      return frac_sno;
    } else {
      return 1.0;
    }
  } else {
    return 1.0;
  }
}


} // namespace ELM

