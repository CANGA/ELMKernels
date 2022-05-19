
#pragma once

namespace ELM {

// from ColumnDataType.F90 and WaterStateType.F90
//-----------------------------------------------------------------------
// set cold-start initial values for select members of col_ws
//-----------------------------------------------------------------------
template <typename ArrayD1>
ACCELERATE
void init_snow_state(const bool& urbpoi, const int& snl, double& h2osno, double& int_snow, double& snow_depth,
                     double& h2osfc, double& h2ocan, double& frac_h2osfc, double& fwet, double& fdry, double& frac_sno,
                     ArrayD1 snw_rds)
{
  using ELMdims::nlevsno;

  // this should/will be intitialized from input data - 0.0 for now
  // if given swe (h2osno): snow_depth = h2osno / bdsno;
  // if given snow depth:  h2osno = snow_depth * bdsno;
  h2osno = 0.0;
  int_snow = 0.0;
  snow_depth = 0.0;
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
    if (snow_depth > 0.0) {
      const double snowbd{std::min(400.0, h2osno / snow_depth)}; // bulk density of snow (kg/m3)
      const double fmelt{pow(snowbd / 100.0, 1.0)};
      // 100 is the assumed fresh snow density; 1 is a melting factor that could be
      // reconsidered, optimal value of 1.5 in Niu et al., 2007
      frac_sno = tanh(snow_depth / (2.5 * ELMconst::ZLND * fmelt));
    }
  }

  // initial snow radius
  if (snl > 0) {
    for (int i = 0; i < nlevsno - snl; ++i) {
      snw_rds(i) = 0.0;
    }
    for (int i = nlevsno - snl; i < nlevsno; ++i) {
      snw_rds(i) = ELMconst::SNW_RDS_MIN;
    }
  } else if (h2osno > 0.0) {
    snw_rds(nlevsno - 1) = ELMconst::SNW_RDS_MIN;
    for (int i = 0; i < nlevsno - 1; ++i) {
      snw_rds(i) = 0.0;
    }
  } else {
    for (int i = 0; i < nlevsno; ++i) {
      snw_rds(i) = 0.0;
    }
  }
}

template <class ArrayD1>
ACCELERATE
void init_snow_layers(const double& snow_depth, const bool& lakpoi, int& snl, ArrayD1 dz, ArrayD1 z, ArrayD1 zi)
{
  using ELMdims::nlevsno;

  for (int i = 0; i < nlevsno; i++) {
    dz[i] = spval;
    z[i] = spval;
    zi[i] = spval;
  }

  if (!lakpoi) {
    if (snow_depth < 0.01) {
      snl = 0;
      for (int i = 0; i < nlevsno; i++) {
        dz[i] = 0.0;
        z[i] = 0.0;
        zi[i] = 0.0;
      }
      zi[nlevsno] = 0.0;
    } else {
      if ((snow_depth >= 0.01) && (snow_depth <= 0.03)) {
        snl = 1;
        dz[4] = snow_depth;
      } else if ((snow_depth > 0.03) && (snow_depth <= 0.04)) {
        snl = 2;
        dz[3] = snow_depth / 2.0;
        dz[4] = dz[3];
      } else if ((snow_depth > 0.04) && (snow_depth <= 0.07)) {
        snl = 2;
        dz[3] = 0.02;
        dz[4] = snow_depth - dz[3];
      } else if ((snow_depth > 0.07) && (snow_depth <= 0.12)) {
        snl = 3;
        dz[2] = 0.02;
        dz[3] = (snow_depth - 0.02) / 2.0;
        dz[4] = dz[3];
      } else if ((snow_depth > 0.12) && (snow_depth <= 0.18)) {
        snl = 3;
        dz[2] = 0.02;
        dz[3] = 0.05;
        dz[4] = snow_depth - dz[2] - dz[3];
      } else if ((snow_depth > 0.18) && (snow_depth <= 0.29)) {
        snl = 4;
        dz[1] = 0.02;
        dz[2] = 0.05;
        dz[3] = (snow_depth - dz[1] - dz[2]) / 2.0;
        dz[4] = dz[3];
      } else if ((snow_depth > 0.29) && (snow_depth <= 0.41)) {
        snl = 4;
        dz[1] = 0.02;
        dz[2] = 0.05;
        dz[3] = 0.11;
        dz[4] = snow_depth - dz[1] - dz[2] - dz[3];
      } else if ((snow_depth > 0.41) && (snow_depth <= 0.64)) {
        snl = 5;
        dz[0] = 0.02;
        dz[1] = 0.05;
        dz[2] = 0.11;
        dz[3] = (snow_depth - dz[0] - dz[1] - dz[2]) / 2.0;
        dz[4] = dz[3];
      } else if (snow_depth > 0.64) {
        snl = 5;
        dz[0] = 0.02;
        dz[1] = 0.05;
        dz[2] = 0.11;
        dz[3] = 0.23;
        dz[4] = snow_depth - dz[0] - dz[1] - dz[2] - dz[3];
      }
    }
    for (int j = nlevsno - 1; j >= nlevsno - snl; j--) {
      z[j] = zi[j + 1] - 0.5 * dz[j];
      zi[j] = zi[j + 1] - dz[j];
    }
  } else {
    snl = 0;
    for (int i = 0; i < nlevsno; i++) {
      dz[i] = 0.0;
      z[i] = 0.0;
      zi[i] = 0.0;
    }
    zi[nlevsno] = 0.0;
  }
}

} // namespace ELM
