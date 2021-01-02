#ifndef FRICTION_VELOCITY_H_
#define FRICTION_VELOCITY_H_

void MoninObukIni(const double& ur, const double& thv, const double& dthv, const double& zldis, const double& z0m, double& um,
double& obu);

void FrictionVelocityWind( const double& forc_hgt_u_patch, const double& displa, const double& um, const double& obu, 
const double& z0m, double& ustar);

void FrictionVelocityTemperature( const double& forc_hgt_t_patch, const double& displa, const double& obu, const double& z0h, 
double& temp1);

void FrictionVelocityHumidity( const double& forc_hgt_q_patch, const double& forc_hgt_t_patch, const double& displa, 
const double& obu, const double& z0h, const double& z0q, const double& temp1, double& temp2);

#endif

