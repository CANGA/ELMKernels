
#include "elm_constants.h"
#include <cmath>

namespace ELM {

void array_normalization(double *arr_inout) {
  double arr_sum = 0.0;

  for (int i = 0; i < nlevgrnd; i++) {
    arr_sum += arr_inout[i];
  }
  for (int i = 0; i < nlevgrnd; i++) {
    if (arr_sum > 0.0) {
      arr_inout[i] /= arr_sum;
    }
  }
}


void soil_suction(const double &smpsat, const double &s, const double &bsw, double &smp) {
  smp = -smpsat * pow(s, (-bsw));
}

} // namespace ELM

