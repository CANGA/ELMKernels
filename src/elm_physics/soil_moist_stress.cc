
#include "soil_moist_stress.h"

namespace ns = ELM::soil_moist_stress;
  
void ns::array_normalization(double *arr_inout) {
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


double ns::soil_suction(const double &smpsat, const double &s, const double &bsw) {
  return -smpsat * pow(s, (-bsw));
}


double ns::dsuction_dsat(const double &bsw, const double &smp, const double &s) {
  return -bsw * smp / s;
}
