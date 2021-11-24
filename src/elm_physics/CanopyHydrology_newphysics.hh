
// h2ocan update - eq 7.9 in CLM 4.5 tech note - 
//h2ocan(n+1) = h2ocan(n) + intercepted_precip * dt - qflx_candrip * dt - evap * dt
// - evap * dt portion gets calculated in CanopyFluxes

#include "array.hh"
#include "elm_constants.h"
#include "landtype.h"
#include <algorithm>
#include <cmath>
//#include "CanopyHydrology_newphysics.hh"

using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

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




void interception_physics(const LandType &Land, const bool &do_capsnow, const int &frac_veg_nosno, 
  const double &forc_rain, const double &forc_snow, const double &dtime, const double &dewmx, const double &elai, 
  const double &esai, double &h2ocan, const double &qflx_irrig, double &qflx_snwcp_liq, 
  double &qflx_snwcp_ice, double &qflx_snow_grnd, double &qflx_rain_grnd, double &qflx_prec_grnd) {

  if (!Land.lakpoi) {
    double qflx_candrip;
    double total_precip = forc_snow + forc_rain;
    
    double fpi = fraction_of_potential_interception(elai, esai);
    
    canopy_storage_and_drainage(Land, frac_veg_nosno, total_precip, 
      fpi, dtime, dewmx, elai, esai, h2ocan, qflx_candrip);
    
    ground_precip_fluxes(Land, do_capsnow, frac_veg_nosno, total_precip, forc_rain, forc_snow, fpi, qflx_candrip, 
      qflx_irrig, qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd, qflx_rain_grnd, qflx_prec_grnd);
  }
}












//double sca_fraction_oldfflag(const double &h2osno, const double &newsnow, const double &snow_depth) {
//  return tanh(snow_depth / (2.5 * zlnd *
//              pow((std::min(800.0, ((h2osno + newsnow) / snow_depth / 100.0))), 1.0))); // why to the power of 1.0??
//}
//
//double sca_fraction_accumulation_bare(const double &newsnow) {
//  return tanh(accum_factor * newsnow);
//}
//
//double sca_fraction_accumulation_snow(const double &newsnow, const double &frac_sno) {
//  return 1.0 - (1.0 - sca_fraction_accumulation_bare(newsnow)) * (1.0 - frac_sno);
//}
//
//double sca_fraction_depletion(const double &h2osno, const double &int_snow, const double &n_melt) {
//  double smr = std::min(1.0, (h2osno / int_snow));
//  return 1.0 - pow((acos(std::min(1.0, (2.0 * smr - 1.0))) / ELM_PI), n_melt);
//}
//
//// I'm not sure what this means - it matches ELM output, but the values are really large - gets used in SnowWater()
//double integrated_snow1(const double &h2osno, const double &newsnow, const double &frac_sno, 
//  const double &n_melt) {
//  double temp_intsnow =
//  (h2osno + newsnow) / (0.5 * (cos(ELM_PI * pow((1.0 - std::max(frac_sno, 1.e-6)), (1.0 / n_melt))) + 1.0));
//  return std::min(1.e8, temp_intsnow);
//}
//
//
//double effective_snow_fraction(const int &ltype, const double &frac_sno) { // set frac_sno_eff variable
//  if (ltype == istsoil || ltype == istcrop) {
//    if (subgridflag == 1) {
//      return frac_sno;
//    } else {
//      return 1.0;
//    }
//  } else {
//    return 1.0;
//  }
//}
//
//
//void save_initial_swe1(const int& snl, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, ArrayD1 swe_old) {
//  // save initial snow content
//  for (int j = 0; j < nlevsno - snl; j++) {
//    swe_old[j] = 0.0;
//  }
//  for (int j = nlevsno - snl; j < nlevsno; j++) {
//    swe_old[j] = h2osoi_liq[j] + h2osoi_ice[j];
//  }
//} 
//
//double save_initial_swe(const int& snl, const double &snow_depth, const ArrayD1 h2osoi_liq,
//  const ArrayD1 h2osoi_ice, double &snow_depth_old, ArrayD1 swe_old) {
//  // save initial snow content
//  for (int j = 0; j < nlevsno - snl; j++) {
//    swe_old[j] = 0.0;
//  }
//  for (int j = nlevsno - snl; j < nlevsno; j++) {
//    swe_old[j] = h2osoi_liq[j] + h2osoi_ice[j];
//  }
//  return snow_depth;
//} 
//
//
//
//double snow_bulk_density(const double &forc_t) {
//  if (forc_t > tfrz + 2.0) {
//    return 50.0 + 1.7 * pow(17.0, 1.5);
//  } else if (forc_t > tfrz - 15.0) {
//    return 50.0 + 1.7 * pow((forc_t - tfrz + 15.0), 1.5);
//  } else {
//    return 50.0;
//  }
//}
//
//double fractional_sca(const bool &do_capsnow, const double &h2osno, const double &snowmelt, const double &int_snow, const double &n_melt, 
//  const double &newsnow, const double &frac_sno) {
//
//  double fsno_new = frac_sno;
//  if (!do_capsnow) {
//    if (h2osno > 0.0) {
//      if (snowmelt > 0.0) { // first compute fsca change from melt during previous time step
//        fsno_new = sca_fraction_depletion(h2osno, int_snow, n_melt);
//      }
//      if (newsnow > 0.0) { // update fsca by new snow event, add to previous fsca
//        fsno_new = sca_fraction_accumulation_snow(newsnow, fsno_new);
//      }
//    } else { // h2osno == 0
//      // initialize fsca when no snow present initially
//      if (newsnow > 0.0) {
//        fsno_new = sca_fraction_accumulation_bare(newsnow);
//      } else {
//        fsno_new = 0.0;
//      }
//    }
//    return fsno_new;
//  }
//  return 1.0; // fully snow-covered if do_capsnow
//}
//
//double integrated_snow(const bool &do_capsnow, const double &h2osno, const double &newsnow, const double &frac_sno, 
//  const double &n_melt, const double &int_snow) {
//
//  return (do_capsnow) ? 5.0e2 : (newsnow > 0.0) ? std::min(1.e8, (h2osno + newsnow) / 
//    (0.5 * (cos(ELM_PI * pow((1.0 - std::max(frac_sno, 1.e-6)), 
//    (1.0 / n_melt))) + 1.0))) + newsnow : int_snow + newsnow;
//}
//
//
//double frost_deposition(const bool &do_capsnow, const double &h2osno, const double &int_snow) {
//  return (do_capsnow) ? int_snow : std::max(int_snow, h2osno);
//}
//
//double dz_snow(const bool &do_capsnow, const double &snow_depth, const double &snow_depth_old) {
//  return (do_capsnow) ? 0.0 : snow_depth - snow_depth_old; // update change in snow depth
//} 
//
//
//double update_h2osno(const bool &do_capsnow, const double &newsnow, const double &h2osno) {
//  return (do_capsnow) ? h2osno : h2osno + newsnow;
//}
//
//double depth_of_snow(const bool &do_capsnow, const bool &urbpoi, const double &h2osno, const double &newsnow, 
//  const double &forc_t, const double &frac_sno, const double &snow_depth) {
//
//  double snow_depth_new = snow_depth;
//  if (!do_capsnow) {
//    if (h2osno > 0.0) {
//      // for subgrid fluxes
//      if (subgridflag == 1 && !urbpoi) {
//        if (frac_sno > 0.0) {
//          snow_depth_new = snow_depth + newsnow / (snow_bulk_density(forc_t) * frac_sno);
//        } else {
//          snow_depth_new = 0.0;
//        }
//      } else {
//        // for uniform snow cover
//        snow_depth_new = snow_depth + newsnow / snow_bulk_density(forc_t);
//      }
//    } else { // h2osno == 0
//      // initialize snow_depth when no snow present initially
//      if (newsnow > 0.0) {
//        // for subgrid fluxes
//        if (subgridflag == 1 && !urbpoi) {
//          double z_avg = newsnow / snow_bulk_density(forc_t);
//          snow_depth_new = z_avg / frac_sno;
//        } else {
//          snow_depth_new = newsnow / snow_bulk_density(forc_t);
//        }
//      } else {
//        snow_depth_new = 0.0;
//      }
//    }
//  }
//  return snow_depth_new;
//}
//
//
//double oldfflag_sca(const bool &do_capsnow, const int &oldfflag, const double &h2osno, const double &snow_depth, const double &newsnow, 
//  const double &frac_sno) {
//
//  double fsno_new = frac_sno;
//  if (!do_capsnow && oldfflag == 1) {
//    if (h2osno > 0.0) {
//      // snow cover fraction in Niu et al. 2007
//      if (snow_depth > 0.0) {
//        fsno_new = sca_fraction_oldfflag(h2osno, newsnow, snow_depth);
//      }
//      if (h2osno < 1.0) {
//        fsno_new = std::min(fsno_new, h2osno);
//      }
//    } else { // h2osno == 0
//      if (newsnow > 0.0) {
//        // snow cover fraction in Niu et al. 2007
//        if (snow_depth > 0.0) {
//          fsno_new = sca_fraction_oldfflag(h2osno, newsnow, snow_depth);
//        }
//      }
//    }
//  }
//  return fsno_new;
//}
//
//
//
////in:
////urbpoi
////h2osno
////snowmelt
////n_melt
////newsnow
////bifall
////oldfflag
////
////out:
////int_snow - state
////frac_sno - state
////snow_depth - state
//
//
//// SCA is based on CLM 4.5 tech note eqs 7.14 (accumulation) and 7.15 (depletion) 
//void update_sca_and_snow_depth(const bool &urbpoi, const int &oldfflag, const double &h2osno, const double &snowmelt,
// const double &n_melt, const double &newsnow, const double &bifall, 
// double &int_snow, double &frac_sno, double &snow_depth) {
//
//  /*======================  FSCA PARAMETERIZATIONS  ======================
//  fsca parameterization based on *changes* in swe
//  first compute change from melt during previous time step */
//  if (h2osno > 0.0) {
//    if (snowmelt > 0.0) {
//      frac_sno = sca_fraction_depletion(h2osno, int_snow, n_melt);
//    }
//    // update fsca by new snow event, add to previous fsca
//    if (newsnow > 0.0) {
//      frac_sno = sca_fraction_accumulation_snow(newsnow, frac_sno);
//      // reset int_snow after accumulation events
//      int_snow = integrated_snow1(h2osno, newsnow, frac_sno, n_melt);
//    }
//    /*====================================================================*/
//    // for subgrid fluxes
//    if (subgridflag == 1 && !urbpoi) {
//      if (frac_sno > 0.0) {
//        snow_depth = snow_depth + newsnow / (bifall * frac_sno);
//      } else {
//        snow_depth = 0.0;
//      }
//    } else {
//      // for uniform snow cover
//      snow_depth = snow_depth + newsnow / bifall;
//    }
//    if (oldfflag == 1) {
//      // snow cover fraction in Niu et al. 2007
//      if (snow_depth > 0.0) {
//        frac_sno = sca_fraction_oldfflag(h2osno, newsnow, snow_depth);
//      }
//      if (h2osno < 1.0) {
//        frac_sno = std::min(frac_sno, h2osno);
//      }
//    }
//  } else { // h2osno == 0
//    // initialize frac_sno and snow_depth when no snow present initially
//    if (newsnow > 0.0) {
//      //frac_sno = 0.0;
//      frac_sno = sca_fraction_accumulation_bare(newsnow);
//      // make int_snow consistent w/ new fsno, h2osno
//      // reset prior to adding newsnow below
//      int_snow = integrated_snow1(h2osno, newsnow, frac_sno, n_melt);
//      // update snow_depth to be consistent with frac_sno, z_avg
//      if (subgridflag == 1 && !urbpoi) {
//        double z_avg = newsnow / bifall;
//        snow_depth = z_avg / frac_sno;
//      } else {
//        snow_depth = newsnow / bifall;
//      }
//      // use n&y07 formulation
//      if (oldfflag == 1) {
//        // snow cover fraction in Niu et al. 2007
//        if (snow_depth > 0.0) {
//          frac_sno = sca_fraction_oldfflag(h2osno, newsnow, snow_depth);
//        }
//      }
//    } else {
//      snow_depth = 0.0;
//      frac_sno = 0.0;
//    }
//  }
//}
//
//
//
//
//
//
//
//
//
//
//
//// eq 7.40
//int initialize_new_snow_layer(const int &ltype, const double &qflx_snow_grnd, const double &frac_sno, 
//  const double &forc_t, const double &t_grnd, double &h2osno, double &snow_depth, int &snl,
//  ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
//  ArrayD1 frac_iceold, ArrayD1 snw_rds) {
//
//  if (ltype == istwet && t_grnd > tfrz) {
//      h2osno = 0.0;
//      snow_depth = 0.0;
//  }
//  // When the snow accumulation exceeds 10 mm, initialize snow layer
//  // Currently, the water temperature for the precipitation is simply set
//  // as the surface air temperature
//  int newnode = 0; // flag for when snow node will be initialized
//  if (snl == 0 && qflx_snow_grnd > 0.0 && (frac_sno * snow_depth) >= 0.01) {
//    newnode = 1;
//    snl = 1;
//    dz[nlevsno-snl] = snow_depth; // meter
//    z[nlevsno-snl] = -0.5 * dz[nlevsno-snl];
//    zi[nlevsno-snl] = -dz[nlevsno-snl];
//    t_soisno[nlevsno-snl] = std::min(tfrz, forc_t); // K
//    h2osoi_ice[nlevsno-snl] = h2osno;               // kg/m2
//    h2osoi_liq[nlevsno-snl] = 0.0;                  // kg/m2
//    frac_iceold[nlevsno-snl] = 1.0;
//    snw_rds[nlevsno-snl] = snw_rds_min;
//  }
//  return newnode;
//}
//
//
//// eq 7.40
//int initialize_new_snow_layer1(const double &qflx_snow_grnd, const double &frac_sno, 
//  const double &snow_depth, const double &forc_t, const double &h2osno, int &snl,
//  ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
//  ArrayD1 frac_iceold, ArrayD1 snw_rds) {
//
//  // When the snow accumulation exceeds 10 mm, initialize snow layer
//  // Currently, the water temperature for the precipitation is simply set
//  // as the surface air temperature
//  int newnode = 0; // flag for when snow node will be initialized
//  if (snl == 0 && qflx_snow_grnd > 0.0 && (frac_sno * snow_depth) >= 0.01) {
//    newnode = 1;
//    snl = 1;
//    dz[nlevsno-snl] = snow_depth; // meter
//    z[nlevsno-snl] = -0.5 * dz[nlevsno-snl];
//    zi[nlevsno-snl] = -dz[nlevsno-snl];
//    t_soisno[nlevsno-snl] = std::min(tfrz, forc_t); // K
//    h2osoi_ice[nlevsno-snl] = h2osno;               // kg/m2
//    h2osoi_liq[nlevsno-snl] = 0.0;                  // kg/m2
//    frac_iceold[nlevsno-snl] = 1.0;
//    snw_rds[nlevsno-snl] = snw_rds_min;
//  }
//  return newnode;
//}
//
//
//void apply_new_snow_to_top_layer(const int &snl, const int &newnode, const double &newsnow, const double &dz_snowf, ArrayD1 h2osoi_ice, ArrayD1 dz) {
//  // The change of ice partial density of surface node due to precipitation.
//  // Only ice part of snowfall is added here, the liquid part will be added later.
//  if (snl > 0 && newnode == 0) {
//    h2osoi_ice[nlevsno-snl] += newsnow;
//    dz[nlevsno-snl] += dz_snowf;
//  }
//}
//
////in:
////Land
////qflx_snow_grnd
////dtime
////do_capsnow
////forc_t
////n_melt
////
////out:
////dz_snowf  - local
////newsnow  - local
////frac_sno - state
////frac_sno_eff - state
////int_snow - state
////h2osno - state
////snow_depth - state
//
//void update_snowcover(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, const double &qflx_snow_melt, const double &dtime, 
//const double &forc_t, const double &t_grnd, const double &n_melt, double &dz_snowf , double &newsnow , double &frac_sno, 
//double &frac_sno_eff, double &int_snow, double &h2osno, double &snow_depth) {
//  if (do_capsnow) {
//    dz_snowf = 0.0;
//    newsnow = qflx_snow_grnd * dtime;
//    frac_sno = 1.0;
//    int_snow = 5.e2;
//  } else {
//    double bifall;
//    if (forc_t > tfrz + 2.0) {
//      bifall = 50.0 + 1.7 * pow(17.0, 1.5);
//    } else if (forc_t > tfrz - 15.0) {
//      bifall = 50.0 + 1.7 * pow((forc_t - tfrz + 15.0), 1.5);
//    } else {
//      bifall = 50.0;
//    }
//    //double bifall = snow_bulk_density(forc_t);
//    // all snow falls on ground, no snow on h2osfc
//    newsnow = qflx_snow_grnd * dtime;
//    // update int_snow
//    int_snow = std::max(int_snow, h2osno); // h2osno could be larger due to frost
//    // snowmelt from previous time step * dtime
//    double snowmelt = qflx_snow_melt * dtime;
//    
//    // save initial snow
//    double temp_snow_depth = snow_depth;
//
//    // calculate frac_sno and snow_depth
//     update_sca_and_snow_depth(Land.urbpoi, oldfflag, h2osno, snowmelt,
//     n_melt, newsnow, bifall, 
//     int_snow, frac_sno, snow_depth);
//
//    //double snowmelt = qflx_snow_melt * dtime;
//    //int_snow = frost_deposition(do_capsnow, h2osno, int_snow);
//    //frac_sno = fractional_sca(do_capsnow, h2osno, snowmelt, int_snow, n_melt, newsnow, frac_sno);
//    //int_snow = integrated_snow(do_capsnow, h2osno, newsnow, frac_sno, n_melt, int_snow);
//    //snow_depth = depth_of_snow(do_capsnow, Land.urbpoi, h2osno, newsnow, forc_t, frac_sno, snow_depth);
//    //frac_sno = oldfflag_sca(do_capsnow, oldfflag, h2osno, snow_depth, newsnow, frac_sno);
//
//    //double dz_snowf = dz_snow(do_capsnow, snow_depth, snow_depth_old);
//
//    h2osno = h2osno + newsnow; // update h2osno for new snow
//    int_snow = int_snow + newsnow;
//    dz_snowf = (snow_depth - temp_snow_depth); // update change in snow depth
//  } // end else do_capsnow
//  if (Land.ltype == istwet && t_grnd > tfrz) {
//    h2osno = 0.0;
//    snow_depth = 0.0;
//  }
//}
//
//void update_snow1(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, const double &qflx_snow_melt, const double &dtime, 
//const double &forc_t, const double &t_grnd, const double &n_melt,
//
//double &frac_sno, double &frac_sno_eff, double &int_snow, double &h2osno, double &snow_depth,
// int &snl, ArrayD1 swe_old, ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
//  ArrayD1 frac_iceold, ArrayD1 snw_rds) {
//
//  if (!Land.lakpoi) {
//    double dz_snowf, newsnow; // local vars
//    save_initial_swe1(snl, h2osoi_liq, h2osoi_ice, swe_old);
//    
//    update_snowcover(Land, do_capsnow, oldfflag, qflx_snow_grnd, qflx_snow_melt, dtime, 
//    forc_t, t_grnd, n_melt, dz_snowf , newsnow , frac_sno, 
//    frac_sno_eff, int_snow, h2osno, snow_depth);
//    
//    frac_sno_eff = effective_snow_fraction(Land.ltype, frac_sno);
//
//    int newnode = initialize_new_snow_layer1(qflx_snow_grnd, frac_sno, 
//    snow_depth, forc_t, h2osno, 
//    snl, dz, z, zi, t_soisno, h2osoi_ice, h2osoi_liq, 
//    frac_iceold, snw_rds);
//
//    apply_new_snow_to_top_layer(snl, newnode, newsnow, dz_snowf, 
//    h2osoi_ice, dz);
//  }
//}
//
//
//void update_snow(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, 
//  const double &qflx_snow_melt, const double &dtime, const double &forc_t, const double &t_grnd, 
//  const double &n_melt, double &frac_sno, double &frac_sno_eff, double &int_snow, double &h2osno, 
//  double &snow_depth, int &snl, ArrayD1 swe_old, ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, 
//  ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, ArrayD1 frac_iceold, ArrayD1 snw_rds) {
//
//  if (!Land.lakpoi) {
//
//    double snow_depth_old = save_initial_swe(snl, snow_depth, h2osoi_liq, h2osoi_ice, snow_depth_old, swe_old);
//    
//    double newsnow = qflx_snow_grnd * dtime;
//
//    double snowmelt = qflx_snow_melt * dtime;
//
//    int_snow = frost_deposition(do_capsnow, h2osno, int_snow);
//
//    frac_sno = fractional_sca(do_capsnow, h2osno, snowmelt, int_snow, n_melt, newsnow, frac_sno);
//
//    int_snow = integrated_snow(do_capsnow, h2osno, newsnow, frac_sno, n_melt, int_snow);
//
//    snow_depth = depth_of_snow(do_capsnow, Land.urbpoi, h2osno, newsnow, forc_t, frac_sno, snow_depth);
//
//    frac_sno = oldfflag_sca(do_capsnow, oldfflag, h2osno, snow_depth, newsnow, frac_sno);
//
//    frac_sno_eff = effective_snow_fraction(Land.ltype, frac_sno);
//
//    h2osno = update_h2osno(do_capsnow, newsnow, h2osno);
//
//    double dz_snowf = dz_snow(do_capsnow, snow_depth, snow_depth_old);
//
//    int newnode = initialize_new_snow_layer(qflx_snow_grnd, frac_sno, snow_depth, forc_t, t_grnd, 
//      h2osno, snow_depth, snl, dz, z, zi, t_soisno, h2osoi_ice, h2osoi_liq, frac_iceold, snw_rds);
//
//    apply_new_snow_to_top_layer(snl, newnode, newsnow, dz_snowf, h2osoi_ice, dz);
//  }
//}



// solves eq 7.67 for frac_h2osfc using d determined via newton-raphson iteration of eq 7.66
double frac_h2osfc_newton_iteration(const double &sigma, const double &h2osfc) {
  double d = 0.0;
  for (int l = 0; l < 10; l++) {
    double fd = 0.5 * d * (1.0 + erf(d / (sigma * sqrt(2.0)))) +
         sigma / sqrt(2.0 * ELM_PI) * exp(-pow(d, 2) / (2.0 * pow(sigma, 2))) - h2osfc;
    double dfdd = 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))));
    d -= fd / dfdd;
  }
  // update the submerged areal fraction using the new d value
  return 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))));
}


void submerged_fraction(const double &micro_sigma, double &h2osfc, double &frac_h2osfc, ArrayD1 h2osoi_liq) {
  // a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
  const double min_h2osfc = 1.e-8; // arbitrary lower limit on h2osfc for safer numerics...
  if (h2osfc > min_h2osfc) {
    // Use newton-raphson method to iteratively determine frac_h20sfc
    // based on amount of surface water storage (h2osfc) and
    // microtopography variability (micro_sigma)
    double sigma = 1.0e3 * micro_sigma; // convert to mm
    frac_h2osfc = frac_h2osfc_newton_iteration(sigma, h2osfc);
  } else {
    frac_h2osfc = 0.0;
    h2osoi_liq[nlevsno] += h2osfc;
    h2osfc = 0.0;
  }
}


void adjust_h2osfc_and_snow_fractions(const double &h2osno, double &frac_h2osfc, double &frac_sno, double &frac_sno_eff) {
// adjust fh2o, fsno when sum is greater than zero
  if (frac_sno > (1.0 - frac_h2osfc) && h2osno > 0.0) {
    if (frac_h2osfc > 0.01) {
      frac_h2osfc = std::max((1.0 - frac_sno), 0.01);
      frac_sno = 1.0 - frac_h2osfc;
    } else {
      frac_sno = 1.0 - frac_h2osfc;
    }
    frac_sno_eff = frac_sno;
  }
}

template <class ArrayD1>
void innundation(const LandType &Land, const double &micro_sigma, const double &h2osno,
                double &h2osfc, ArrayD1 h2osoi_liq, double &frac_sno, double &frac_sno_eff, double &frac_h2osfc) {

  if (!Land.lakpoi) {
    // h2osfc only calculated for soil vegetated land units
    if (Land.ltype == istsoil || Land.ltype == istcrop) {
      submerged_fraction(micro_sigma, h2osfc, frac_h2osfc, h2osoi_liq);
      adjust_h2osfc_and_snow_fractions(h2osno, frac_h2osfc, frac_sno, frac_sno_eff);
    } else { // if landunit not istsoil/istcrop, set frac_h2osfc to zero
      frac_h2osfc = 0.0;
    }
  }
} // FracH2OSfc












} // namespace ELM

