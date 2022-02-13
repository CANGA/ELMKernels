/*! \file qsat.h
\brief Internal function derived from QSatMod.F90
*/
#pragma once

namespace ELM::qsat {

/*! Computes saturation mixing ratio and the change in saturation. (internal) */
void qsat(const double &T, const double &p, double &es, double &esdT, double &qs, double &qsdT);

} // namespace ELM::qsat

