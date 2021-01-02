#ifndef QSAT_H_
#define QSAT_H_

void QSat(const double& T, const double& p, double& es, double& esdT, double& qs, double& qsdT);

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

#endif
