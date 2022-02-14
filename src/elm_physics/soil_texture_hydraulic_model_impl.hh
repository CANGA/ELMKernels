
#pragma once

namespace ELM::soil_hydraulics {

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

} // namespace ELM::soil_hydraulics

