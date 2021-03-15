/*! \file QSat.hh
\brief Internal function derived from QSatMod.F90
*/
#pragma once

namespace ELM {

/*! Computes saturation mixing ratio and the change in saturation. (internal) */
void QSat(const double &T, const double &p, double &es, double &esdT, double &qs, double &qsdT);

} // namespace ELM
