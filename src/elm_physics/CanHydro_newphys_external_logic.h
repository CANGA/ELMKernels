
#pragma once

#include "array.hh"
#include "elm_constants.h"
#include "LandType.h"
#include <algorithm>
#include <cmath>

using ArrayB1 = ELM::Array<bool, 1>;
using ArrayI1 = ELM::Array<int, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

namespace ELM {
namespace model_external {



double sca_fraction_oldfflag(const double &h2osno, const double &newsnow, const double &snow_depth);

double sca_fraction_accumulation_bare(const double &newsnow);

double sca_fraction_accumulation_snow(const double &newsnow, const double &frac_sno);

double sca_fraction_depletion(const double &h2osno, const double &int_snow, const double &n_melt);

// I'm not sure what this means - it matches ELM output, but the values are really large - gets used in SnowWater()
double integrated_snow1(const double &h2osno, const double &newsnow, const double &frac_sno, 
  const double &n_melt);


void save_initial_swe1(const int& snl, const ArrayD1 h2osoi_liq, const ArrayD1 h2osoi_ice, ArrayD1 swe_old);



double effective_snow_fraction(const int &ltype, const double &frac_sno);


// eq 7.40
int initialize_new_snow_layer1(const double &qflx_snow_grnd, const double &frac_sno, 
  const double &snow_depth, const double &forc_t, const double &h2osno, int &snl,
  ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
  ArrayD1 frac_iceold, ArrayD1 snw_rds);


void apply_new_snow_to_top_layer(const int &snl, const int &newnode, const double &newsnow, 
  const double &dz_snowf, ArrayD1 h2osoi_ice, ArrayD1 dz);



// SCA is based on CLM 4.5 tech note eqs 7.14 (accumulation) and 7.15 (depletion) 
void update_sca_and_snow_depth(const bool &urbpoi, const int &oldfflag, const double &h2osno, const double &snowmelt,
 const double &n_melt, const double &newsnow, const double &bifall, 
 double &int_snow, double &frac_sno, double &snow_depth);



void update_snowcover(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, const double &qflx_snow_melt, const double &dtime, 
const double &forc_t, const double &t_grnd, const double &n_melt, double &dz_snowf , double &newsnow , double &frac_sno, 
double &frac_sno_eff, double &int_snow, double &h2osno, double &snow_depth);




void update_snow1(const LandType &Land, const bool &do_capsnow, const int &oldfflag, const double &qflx_snow_grnd, const double &qflx_snow_melt, const double &dtime, 
const double &forc_t, const double &t_grnd, const double &n_melt,

double &frac_sno, double &frac_sno_eff, double &int_snow, double &h2osno, double &snow_depth,
 int &snl, ArrayD1 swe_old, ArrayD1 dz, ArrayD1 z, ArrayD1 zi, ArrayD1 t_soisno, ArrayD1 h2osoi_ice, ArrayD1 h2osoi_liq, 
  ArrayD1 frac_iceold, ArrayD1 snw_rds);









} //model_external
} // namespace ELM


