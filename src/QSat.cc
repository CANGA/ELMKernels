/*


*/
#include "clm_constants.h"
#include "qsat.h"

void QSat(
  const double& T,
  const double& p,
  
  double& es,
  double& esdT,
  double& qs,
  double& qsdT)
{
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
