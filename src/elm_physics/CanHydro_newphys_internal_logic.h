
#pragma once

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



double sca_fraction_oldfflag(const double &h2osno, const double &newsnow, const double &snow_depth);

double sca_fraction_accumulation_bare(const double &newsnow);

double sca_fraction_accumulation_snow(const double &newsnow, const double &frac_sno);

double sca_fraction_depletion(const double &h2osno, const double &int_snow, const double &n_melt);


double save_initial_swe(const int& snl, const double &snow_depth, const ArrayD1 h2osoi_liq,
  const ArrayD1 h2osoi_ice, double &snow_depth_old, ArrayD1 swe_old);

double effective_snow_fraction(const int &ltype, const double &frac_sno);

double snow_bulk_density(const double &forc_t);

double integrated_snow(const bool &do_capsnow, const double &h2osno, const double &newsnow, const double &frac_sno, 
  const double &n_melt, const double &int_snow);


double frost_deposition(const bool &do_capsnow, const double &h2osno, const double &int_snow);

double dz_snow(const bool &do_capsnow, const double &snow_depth, const double &snow_depth_old);


double update_h2osno(const bool &do_capsnow, const double &newsnow, const double &h2osno);

double fractional_sca(const bool &do_capsnow, const double &h2osno, const double &snowmelt, 
  const double &int_snow, const double &n_melt, const double &newsnow, const double &frac_sno);


double depth_of_snow(const bool &do_capsnow, const bool &urbpoi, const double &h2osno, const double &newsnow, 
  const double &forc_t, const double &frac_sno, const double &snow_depth);

double oldfflag_sca(const bool &do_capsnow, const int &oldfflag, const double &h2osno, const double &snow_depth, const double &newsnow, 
  const double &frac_sno);



// eq 7.40
int initialize_new_snow_layer(const int &ltype, const double &qflx_snow_grnd, const double &frac_sno, 
  const double &forc_t, const double &t_grnd, double &h2osno, double &snow_depth, int &snl,
  ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
  ArrayD1 frac_iceold, ArrayD1 snw_rds);



void apply_new_snow_to_top_layer(const int &snl, const int &newnode, const double &newsnow, 
  const double &dz_snowf, ArrayD1 h2osoi_ice, ArrayD1 dz);



void update_snow(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, 
  const double &qflx_snow_melt, const double &dtime, const double &forc_t, const double &t_grnd, 
  const double &n_melt, double &frac_sno, double &frac_sno_eff, double &int_snow, double &h2osno, 
  double &snow_depth, int &snl, ArrayD1 swe_old, ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, 
  ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, ArrayD1 frac_iceold, ArrayD1 snw_rds);








} //model_internal
} // namespace ELM


