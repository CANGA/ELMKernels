/*
DESCRIPTION:
Initialization snow layer(s) if the snow accumulation exceeds 10 mm.

INPUTS:


OUTPUTS:

*/

#include <algorithm>
#include <cmath>
#include "clm_constants.hh"

void SnowInit(
  const double& dtime,
  const int& ltype,
  const bool& urbpoi,
  const bool& do_capsnow,                            
  const int& oldfflag,
  const double& forc_t,
  const double& t_grnd,
  const double& qflx_snow_grnd,
  const double& qflx_snow_melt,
  const double& n_melt,
  
  double& snow_depth,
  double& h2osno,
  double& int_snow,
  double* swe_old,
  double* h2osoi_liq,
  double* h2osoi_ice,
  double* t_soisno,
  double* frac_iceold,
  int& snl,
  double* dz,
  double* z,
  double* zi,
  double* snw_rds,
  int& newnode,
  double& qflx_snow_h2osfc,
  double& frac_sno_eff,
  double& frac_sno)


{
  double rpi = SHR_CONST_PI;
  double temp_snow_depth, dz_snowf, newsnow, bifall, snowmelt, accum_factor, temp_intsnow, z_avg;
  
  // Determine snow height and snow water
  // Use Alta relationship, Anderson(1976); LaChapelle(1961),
  // U.S.Department of Agriculture Forest Service, Project F,
  // Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.
  qflx_snow_h2osfc = 0.0;
  // set temporary variables prior to updating
  temp_snow_depth = snow_depth;
  // save initial snow content

  for(int j = nlevsno - 1; j >= snl; j--) {
    swe_old[j] = 0.0;
  }
  for (int j = 0; j < snl; j++) {
    swe_old[j] = h2osoi_liq[j] + h2osoi_ice[j];
  }

  if (do_capsnow) {
    dz_snowf = 0.0;
    newsnow = qflx_snow_grnd * dtime;
    frac_sno = 1.0;
    int_snow = 5.e2;
  } else {
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
    int_snow = std::max(int_snow, h2osno); //h2osno could be larger due to frost

    // snowmelt from previous time step * dtime
    snowmelt = qflx_snow_melt * dtime;

    // set shape factor for accumulation of snow
    accum_factor = 0.1;

    /*======================  FSCA PARAMETERIZATIONS  ======================
    fsca parameterization based on *changes* in swe
    first compute change from melt during previous time step */

    if (h2osno > 0.0) {
      if(snowmelt > 0.0) {
        double smr = std::min(1.0, (h2osno / int_snow));
        frac_sno = 1.0 - pow((acos(std::min(1.0, (2.0 * smr - 1.0))) / rpi), n_melt);
      }
      // update fsca by new snow event, add to previous fsca
      if (newsnow > 0.0) {
        double fsno_new = 1.0 - (1.0 - tanh(accum_factor * newsnow)) * (1.0 - frac_sno);
        frac_sno = fsno_new;
        // reset int_snow after accumulation events
        temp_intsnow = (h2osno + newsnow) / (0.5 * (cos(rpi * pow( (1.0 - std::max(frac_sno, 1.e-6)), (1.0 / n_melt))) + 1.0));
        int_snow = std::min(1.e8, temp_intsnow);
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
          frac_sno = tanh(snow_depth / (2.5 * zlnd * pow((std::min(800.0, ((h2osno+ newsnow) / snow_depth / 100.0))), 1.0))); // why to the power of 1.0??
        }
        if (h2osno < 1.0) {
          frac_sno = std::min(frac_sno, h2osno);
        }
      }

    } else { //h2osno == 0
      // initialize frac_sno and snow_depth when no snow present initially
      if (newsnow > 0.0) {
        z_avg = newsnow / bifall;
        frac_sno = tanh(accum_factor * newsnow);
        // make int_snow consistent w/ new fsno, h2osno
        int_snow = 0.0; //reset prior to adding newsnow below
        temp_intsnow = (h2osno + newsnow) / (0.5 * (cos(rpi * pow( (1.0 - std::max(frac_sno, 1.e-6)), (1.0 / n_melt))) + 1.0));
        int_snow = std::min(1.e8, temp_intsnow);

        // update snow_depth and h2osno to be consistent with frac_sno, z_avg
        if (subgridflag == 1 && !urbpoi) {
          snow_depth = z_avg / frac_sno;
        } else {
          snow_depth = newsnow / bifall;
        }

        // use n&y07 formulation
        if (oldfflag == 1) {
          // snow cover fraction in Niu et al. 2007
          if(snow_depth > 0.0) {
            frac_sno = tanh(snow_depth / (2.5 * zlnd * pow((std::min(800.0, ((h2osno+ newsnow) / snow_depth / 100.0))), 1.0))); // why to the power of 1.0??
          }
        }
      } else {
        z_avg = 0.0;
        snow_depth = 0.0;
        frac_sno = 0.0;
      }
    }

    // no snow on surface water
    qflx_snow_h2osfc = 0.0;

    // update h2osno for new snow
    h2osno = h2osno + newsnow;
    int_snow = int_snow + newsnow;

    // update change in snow depth
    dz_snowf = (snow_depth - temp_snow_depth) / dtime;


  } // end else do_capsnow

  // set frac_sno_eff variable
  if (ltype == istsoil || ltype == istcrop) {
    if (subgridflag == 1) {
      frac_sno_eff = frac_sno;
    } else {
      frac_sno_eff = 1.0;
    }
  } else {
    frac_sno_eff = 1.0;
  }

  if (ltype == istwet && t_grnd > tfrz) {
    h2osno = 0.0;
    snow_depth = 0.0;
  }

  // When the snow accumulation exceeds 10 mm, initialize snow layer
  // Currently, the water temperature for the precipitation is simply set
  // as the surface air temperature
  newnode = 0; // flag for when snow node will be initialized
  if (snl == 0 && qflx_snow_grnd > 0.0 && (frac_sno * snow_depth) >= 0.01) {
    newnode = 1;
    snl = 1;
    dz[0] = snow_depth;                     // meter
    z[0] = -0.5 * dz[0];
    zi[1] = -dz[0];
    t_soisno[0] = std::min(tfrz, forc_t);   // K
    h2osoi_ice[0] = h2osno;                 // kg/m2
    h2osoi_liq[0] = 0.0;                    // kg/m2
    frac_iceold[0] = 1.0;
    // intitialize SNICAR variables for fresh snow: -- after aeorosol functionality added
    //call aerosol_vars%Reset(column=c)
    // call waterstate_vars%Reset(column=c)
    snw_rds[0] = snw_rds_min;
  }

  // The change of ice partial density of surface node due to precipitation.
  // Only ice part of snowfall is added here, the liquid part will be added
  // later.
  if (snl > 0 && newnode == 0) {
    h2osoi_ice[snl-1] = h2osoi_ice[snl-1] + newsnow;
    dz[snl-1] = dz[snl-1] + dz_snowf * dtime;
  }
} //SnowInit
