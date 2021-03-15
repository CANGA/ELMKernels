// derived from initVerticalMod.F90

#pragma once

#include "ELMConstants.h"

namespace ELM {

template <class ArrayD1>
void InitSnowLayers(const double &snow_depth, const bool &lakpoi, int &snl, ArrayD1 dz, ArrayD1 z, ArrayD1 zi) {

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
