/*


*/
#include "clm_constants.hh"

void QSat(
  const double& T,
  const double& p,
  
  double& es,
  double& esdT,
  double& qs,
  double& qsdT)
{
  // For water vapor (temperature range 0C-100C)
  const double a0 =  6.11213476;
  const double a1 =  0.444007856;
  const double a2 =  0.143064234e-01;
  const double a3 =  0.264461437e-03;
  const double a4 =  0.305903558e-05;
  const double a5 =  0.196237241e-07;
  const double a6 =  0.892344772e-10;
  const double a7 = -0.373208410e-12;
  const double a8 =  0.209339997e-15;
  // For derivative:water vapor
  const double b0 =  0.444017302;
  const double b1 =  0.286064092e-01;
  const double b2 =  0.794683137e-03;
  const double b3 =  0.121211669e-04;
  const double b4 =  0.103354611e-06;
  const double b5 =  0.404125005e-09;
  const double b6 = -0.788037859e-12;
  const double b7 = -0.114596802e-13;
  const double b8 =  0.381294516e-16;
  // For ice (temperature range -75C-0C)
  const double c0 =  6.11123516;
  const double c1 =  0.503109514;
  const double c2 =  0.188369801e-01;
  const double c3 =  0.420547422e-03;
  const double c4 =  0.614396778e-05;
  const double c5 =  0.602780717e-07;
  const double c6 =  0.387940929e-09;
  const double c7 =  0.149436277e-11;
  const double c8 =  0.262655803e-14;
  // For derivative:ice
  const double d0 =  0.503277922;
  const double d1 =  0.377289173e-01;
  const double d2 =  0.126801703e-02;
  const double d3 =  0.249468427e-04;
  const double d4 =  0.313703411e-06;
  const double d5 =  0.257180651e-08;
  const double d6 =  0.133268878e-10;
  const double d7 =  0.394116744e-13;
  const double d8 =  0.498070196e-16;

  double td,vp,vp1,vp2,T_limit;
  T_limit = T - tfrz;

  if (T_limit > 100.0) { T_limit=100.0; }
  if (T_limit < -75.0) { T_limit=-75.0; }

  td = T_limit;
  if (td >= 0.0) {
    es = a0 + td *(a1 + td*(a2 + td*(a3 + td*(a4
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))));
            
    esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))));
  } else {
    es = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))));

    esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))));
  }
  es    = es    * 100.0;           // pa
  esdT  = esdT  * 100.0;           // pa/K
  vp    = 1.0 / (p - 0.378 * es);
  vp1   = 0.622 * vp;
  vp2   = vp1   * vp;
  qs    = es    * vp1;             // kg/kg
  qsdT  = esdT  * vp2 * p;         // 1 / K
}
