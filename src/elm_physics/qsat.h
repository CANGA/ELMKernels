/*! \file qsat.h
\brief Internal function derived from QSatMod.F90
*/
#pragma once

#include "elm_constants.h"

#include "kokkos_includes.hh"

namespace ELM {

/*! Computes saturation mixing ratio and the change in saturation. (internal) */
ACCELERATE
void qsat(const double& T, const double& p, double& es, double& esdT, double& qs, double& qsdT);

} // namespace ELM

#include "qsat_impl.hh"
