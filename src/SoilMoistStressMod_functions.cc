

void calc_effective_soilporosity(
  const double* watsat,
  const double* h2osoi_ice,
  const double* dz,

  double* eff_por)
{
  double vol_ice;

  for (int i = 0; i < nlevgrnd; i++) {
    //compute the volumetric ice content
    vol_ice = std::min(watsat[i], (h2osoi_ice[nlevsno + i] / (denice * dz[nlevsno + i])));
    // compute the maximum soil space to fill liquid water and air
    eff_por[i] = watsat[i] - vol_ice;
  }

}


void calc_volumetric_h2oliq(
  const double* eff_por,
  const double* h2osoi_liq,
  const double* dz,

  double* vol_liq)
{
  for (int i = 0; i < nlevgrnd; i++) {
    // volume of liquid is no greater than effective void space
    vol_liq[i] = std::min(eff_por[i], (h2osoi_liq[nlevsno + i] / (dz[h2osoi_liq] * denh2o)));
  }
}

// CLM call tree:
// calc_root_moist_stress() -> normalize_unfrozen_rootfr() -> array_normalization()
//                          \ 
//                            calc_root_moist_stress_clm45default() -> soil_water_retention_curve%soil_suction() -- this is located in SoilWaterRetentionCurveClappHornberg1978Mod.F90
void calc_root_moist_stress(

  )
{

}

void normalize_unfrozen_rootfr()
{
  
}

void array_normalization (double* arr_inout)
{
  double arr_sum = 0.0;

  for (int i = 0; i < nlevgrnd) {
    arr_sum += arr_inout[i];
  }

  for (int i = 0; i < nlevgrnd) {
    arr_inout[i] = arr_inout[i] / arr_sum;
  }
}