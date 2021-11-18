


#include "array.hh"
#include "ELMConstants.h"
#include "LandType.h"
#include <algorithm>
#include <cmath>

using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

namespace ELM {
namespace model_internal {



double sca_fraction_oldfflag(const double &h2osno, const double &newsnow, const double &snow_depth) {
  return tanh(snow_depth / (2.5 * zlnd *
              pow((std::min(800.0, ((h2osno + newsnow) / snow_depth / 100.0))), 1.0))); // why to the power of 1.0??
}

double sca_fraction_accumulation_bare(const double &newsnow) {
  return tanh(accum_factor * newsnow);
}

double sca_fraction_accumulation_snow(const double &newsnow, const double &frac_sno) {
  return 1.0 - (1.0 - sca_fraction_accumulation_bare(newsnow)) * (1.0 - frac_sno);
}

double sca_fraction_depletion(const double &h2osno, const double &int_snow, const double &n_melt) {
  double smr = std::min(1.0, (h2osno / int_snow));
  return 1.0 - pow((acos(std::min(1.0, (2.0 * smr - 1.0))) / ELM_PI), n_melt);
}


double save_initial_swe(const int& snl, const double &snow_depth, const ArrayD1 h2osoi_liq,
  const ArrayD1 h2osoi_ice, double &snow_depth_old, ArrayD1 swe_old) {
  // save initial snow content
  for (int j = 0; j < nlevsno - snl; j++) {
    swe_old[j] = 0.0;
  }
  for (int j = nlevsno - snl; j < nlevsno; j++) {
    swe_old[j] = h2osoi_liq[j] + h2osoi_ice[j];
  }
  return snow_depth;
} 

double effective_snow_fraction(const int &ltype, const double &frac_sno) { // set frac_sno_eff variable
  if (ltype == istsoil || ltype == istcrop) {
    if (subgridflag == 1) {
      return frac_sno;
    } else {
      return 1.0;
    }
  } else {
    return 1.0;
  }
}

double snow_bulk_density(const double &forc_t) {
  if (forc_t > tfrz + 2.0) {
    return 50.0 + 1.7 * pow(17.0, 1.5);
  } else if (forc_t > tfrz - 15.0) {
    return 50.0 + 1.7 * pow((forc_t - tfrz + 15.0), 1.5);
  } else {
    return 50.0;
  }
}

double integrated_snow(const bool &do_capsnow, const double &h2osno, const double &newsnow, const double &frac_sno, 
  const double &n_melt, const double &int_snow) {

  return (do_capsnow) ? 5.0e2 : (newsnow > 0.0) ? std::min(1.e8, (h2osno + newsnow) / 
    (0.5 * (cos(ELM_PI * pow((1.0 - std::max(frac_sno, 1.e-6)), 
    (1.0 / n_melt))) + 1.0))) + newsnow : int_snow + newsnow;
}


double frost_deposition(const bool &do_capsnow, const double &h2osno, const double &int_snow) {
  return (do_capsnow) ? int_snow : std::max(int_snow, h2osno);
}

double dz_snow(const bool &do_capsnow, const double &snow_depth, const double &snow_depth_old) {
  return (do_capsnow) ? 0.0 : snow_depth - snow_depth_old; // update change in snow depth
} 


double update_h2osno(const bool &do_capsnow, const double &newsnow, const double &h2osno) {
  return (do_capsnow) ? h2osno : h2osno + newsnow;
}

double fractional_sca(const bool &do_capsnow, const double &h2osno, const double &snowmelt, 
  const double &int_snow, const double &n_melt, const double &newsnow, const double &frac_sno) {

  double fsno_new = frac_sno;
  if (!do_capsnow) {
    if (h2osno > 0.0) {
      if (snowmelt > 0.0) { // first compute fsca change from melt during previous time step
        fsno_new = sca_fraction_depletion(h2osno, int_snow, n_melt);
      }
      if (newsnow > 0.0) { // update fsca by new snow event, add to previous fsca
        fsno_new = sca_fraction_accumulation_snow(newsnow, fsno_new);
      }
    } else { // h2osno == 0
      // initialize fsca when no snow present initially
      if (newsnow > 0.0) {
        fsno_new = sca_fraction_accumulation_bare(newsnow);
      } else {
        fsno_new = 0.0;
      }
    }
    return fsno_new;
  }
  return 1.0; // fully snow-covered if do_capsnow
}


double depth_of_snow(const bool &do_capsnow, const bool &urbpoi, const double &h2osno, const double &newsnow, 
  const double &forc_t, const double &frac_sno, const double &snow_depth) {

  double snow_depth_new = snow_depth;
  if (!do_capsnow) {
    if (h2osno > 0.0) {
      // for subgrid fluxes
      if (subgridflag == 1 && !urbpoi) {
        if (frac_sno > 0.0) {
          snow_depth_new = snow_depth + newsnow / (snow_bulk_density(forc_t) * frac_sno);
        } else {
          snow_depth_new = 0.0;
        }
      } else {
        // for uniform snow cover
        snow_depth_new = snow_depth + newsnow / snow_bulk_density(forc_t);
      }
    } else { // h2osno == 0
      // initialize snow_depth when no snow present initially
      if (newsnow > 0.0) {
        // for subgrid fluxes
        if (subgridflag == 1 && !urbpoi) {
          double z_avg = newsnow / snow_bulk_density(forc_t);
          snow_depth_new = z_avg / frac_sno;
        } else {
          snow_depth_new = newsnow / snow_bulk_density(forc_t);
        }
      } else {
        snow_depth_new = 0.0;
      }
    }
  }
  return snow_depth_new;
}


double oldfflag_sca(const bool &do_capsnow, const int &oldfflag, const double &h2osno, const double &snow_depth, const double &newsnow, 
  const double &frac_sno) {

  double fsno_new = frac_sno;
  if (!do_capsnow && oldfflag == 1) {
    if (h2osno > 0.0) {
      // snow cover fraction in Niu et al. 2007
      if (snow_depth > 0.0) {
        fsno_new = sca_fraction_oldfflag(h2osno, newsnow, snow_depth);
      }
      if (h2osno < 1.0) {
        fsno_new = std::min(fsno_new, h2osno);
      }
    } else { // h2osno == 0
      if (newsnow > 0.0) {
        // snow cover fraction in Niu et al. 2007
        if (snow_depth > 0.0) {
          fsno_new = sca_fraction_oldfflag(h2osno, newsnow, snow_depth);
        }
      }
    }
  }
  return fsno_new;
}



// eq 7.40
int initialize_new_snow_layer(const int &ltype, const double &qflx_snow_grnd, const double &frac_sno, 
  const double &forc_t, const double &t_grnd, double &h2osno, double &snow_depth, int &snl,
  ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
  ArrayD1 frac_iceold, ArrayD1 snw_rds) {

  if (ltype == istwet && t_grnd > tfrz) {
      h2osno = 0.0;
      snow_depth = 0.0;
  }
  // When the snow accumulation exceeds 10 mm, initialize snow layer
  // Currently, the water temperature for the precipitation is simply set
  // as the surface air temperature
  int newnode = 0; // flag for when snow node will be initialized
  if (snl == 0 && qflx_snow_grnd > 0.0 && (frac_sno * snow_depth) >= 0.01) {
    newnode = 1;
    snl = 1;
    dz[nlevsno-snl] = snow_depth; // meter
    z[nlevsno-snl] = -0.5 * dz[nlevsno-snl];
    zi[nlevsno-snl] = -dz[nlevsno-snl];
    t_soisno[nlevsno-snl] = std::min(tfrz, forc_t); // K
    h2osoi_ice[nlevsno-snl] = h2osno;               // kg/m2
    h2osoi_liq[nlevsno-snl] = 0.0;                  // kg/m2
    frac_iceold[nlevsno-snl] = 1.0;
    snw_rds[nlevsno-snl] = snw_rds_min;
  }
  return newnode;
}



void apply_new_snow_to_top_layer(const int &snl, const int &newnode, const double &newsnow, 
  const double &dz_snowf, ArrayD1 h2osoi_ice, ArrayD1 dz) {
  // The change of ice partial density of surface node due to precipitation.
  // Only ice part of snowfall is added here, the liquid part will be added later.
  if (snl > 0 && newnode == 0) {
    h2osoi_ice[nlevsno-snl] += newsnow;
    dz[nlevsno-snl] += dz_snowf;
  }
}



void update_snow(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, 
  const double &qflx_snow_melt, const double &dtime, const double &forc_t, const double &t_grnd, 
  const double &n_melt, double &frac_sno, double &frac_sno_eff, double &int_snow, double &h2osno, 
  double &snow_depth, int &snl, ArrayD1 swe_old, ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, 
  ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, ArrayD1 frac_iceold, ArrayD1 snw_rds) {

  if (!Land.lakpoi) {

    double snow_depth_old = save_initial_swe(snl, snow_depth, h2osoi_liq, h2osoi_ice, snow_depth_old, swe_old);
    
    double newsnow = qflx_snow_grnd * dtime;

    double snowmelt = qflx_snow_melt * dtime;

    int_snow = frost_deposition(do_capsnow, h2osno, int_snow);

    frac_sno = fractional_sca(do_capsnow, h2osno, snowmelt, int_snow, n_melt, newsnow, frac_sno);

    int_snow = integrated_snow(do_capsnow, h2osno, newsnow, frac_sno, n_melt, int_snow);

    snow_depth = depth_of_snow(do_capsnow, Land.urbpoi, h2osno, newsnow, forc_t, frac_sno, snow_depth);

    frac_sno = oldfflag_sca(do_capsnow, oldfflag, h2osno, snow_depth, newsnow, frac_sno);

    frac_sno_eff = effective_snow_fraction(Land.ltype, frac_sno);

    h2osno = update_h2osno(do_capsnow, newsnow, h2osno);

    double dz_snowf = dz_snow(do_capsnow, snow_depth, snow_depth_old);

    int newnode = initialize_new_snow_layer(qflx_snow_grnd, frac_sno, snow_depth, forc_t, t_grnd, 
      h2osno, snow_depth, snl, dz, z, zi, t_soisno, h2osoi_ice, h2osoi_liq, frac_iceold, snw_rds);

    apply_new_snow_to_top_layer(snl, newnode, newsnow, dz_snowf, h2osoi_ice, dz);
  }
}








} //model_internal
} // namespace ELM


