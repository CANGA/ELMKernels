/*
functions derived from SoilMoistStressMod.F90

These functions get called from CanopyFluxes.
ELM calc_root_moist_stress() call tree:
calc_root_moist_stress() -> normalize_unfrozen_rootfr() -> array_normalization()
                            >
                              calc_root_moist_stress_clm45default() -> soil_water_retention_curve%soil_suction()

I call calc_root_moist_stress() -> normalize_unfrozen_rootfr() -> array_normalization()
                                   >
                                     soil_suction()


-- note: rootfr_unf[nlevgrnd] needs to be initialized to

soil_suction is SoilWaterRetentionCurveClappHornberg1978Mod, currently without derivative
*/
#pragma once

#include "elm_constants.h"
#include <algorithm>
#include <cmath>

namespace ELM::soil_moist_stress {

/*
DESCRIPTION: normalize array elements with respect to array sum
IN/OUT:
arr_inout[>=nlevgrnd]  [double] double array to normalize
*/
void array_normalization(double *arr_inout);

/*
DESCRIPTION: compute soil suction potential, negative
INPUTS:
smpsat   [double] minimum soil suction, positive [mm]
s        [double] relative saturation [0-1]
bsw      [double] shape parameter

OUTPUTS:
[double] soil suction, negative, [mm]
*/
double soil_suction(const double& smpsat, const double& s, const double& bsw);

/*
bsw      [double] shape parameter
smp      [double] soil suction, negative, [mm]
s        [double] relative saturation [0-1]

OUTPUTS:
[double] d(smp)/ds [mm]
*/
double dsuction_dsat(const double& bsw, const double& smp, const double& s);

/*
DESCRIPTION: normalize root fraction for total unfrozen depth

INPUTS:
t_soisno[nlevgrnd+nlevsno] [double] soil temperature (Kelvin)
rootfr[nlevgrnd]           [double] fraction of roots in each soil layer
altmax_indx                [int] index corresponding to maximum active layer depth from current year
altmax_lastyear_indx       [int] index corresponding to maximum active layer depth from prior year

OUTPUTS:
rootfr_unf[nlevgrnd]       [double] root fraction defined for unfrozen layers only
*/
template <class ArrayD1>
void normalize_unfrozen_rootfr(const ArrayD1 t_soisno, const ArrayD1 rootfr, const int& altmax_indx,
                               const int& altmax_lastyear_indx, double *rootfr_unf);

/*
DESCRIPTION: compute the effective soil porosity

INPUTS:
watsat[nlevgrnd]             [double] volumetric soil water at saturation (porosity)
h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
eff_porosity[nlevgrnd]       [double] effective soil porosity
*/
template <class ArrayD1>
void calc_effective_soilporosity(const ArrayD1 watsat, const ArrayD1 h2osoi_ice, const ArrayD1 dz, ArrayD1 eff_por);

/*
DESCRIPTION: compute the volumetric liquid water content

INPUTS:
eff_porosity[nlevgrnd]       [double] effective soil porosity
h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
dz[nlevgrnd+nlevsno]         [double] layer thickness (m)

OUTPUTS:
vol_liq[nlevgrnd+nlevsno]    [double] volumetric liquid water content
*/
template <class ArrayD1>
void calc_volumetric_h2oliq(const ArrayD1 eff_por, const ArrayD1 h2osoi_liq, const ArrayD1 dz, double *vol_liq);

/*
DESCRIPTION: compute integrated soil water stress (btran), also effective root fraction
INPUTS:
h2osoi_liqvol[nlevgrnd+nlevsno] [double] liquid volumetric moisture
rootfr[nlevgrnd]                [double] fraction of roots in each soil layer
t_soisno[nlevgrnd+nlevsno]      [double] col soil temperature (Kelvin)
tc_stress                       [double] critical soil temperature for soil water stress (C)
sucsat[nlevgrnd]                [double] minimum soil suction (mm)
watsat[nlevgrnd]                [double] volumetric soil water at saturation (porosity)
bsw[nlevgrnd]                   [double] Clapp and Hornberger "b
smpso[numpft]                   [double] soil water potential at full stomatal opening (mm)
smpsc[numpft]                   [double] soil water potential at full stomatal closure (mm)
eff_porosity[nlevgrnd]          [double] effective soil porosity
altmax_indx                     [int] index corresponding to maximum active layer depth from current year
altmax_lastyear_indx            [int] index corresponding to maximum active layer depth from prior year

OUTPUTS:
rootr[nlevgrnd]                 [double] effective fraction of roots in each soil layer
btran                           [double] transpiration wetness factor (0 to 1) (integrated soil water stress)
*/
template <class ArrayD1>
void calc_root_moist_stress(const double *h2osoi_liqvol, const ArrayD1 rootfr, const ArrayD1 t_soisno,
                            const double& tc_stress, const ArrayD1 sucsat, const ArrayD1 watsat, const ArrayD1 bsw,
                            const double& smpso, const double& smpsc, const ArrayD1 eff_porosity,
                            const int& altmax_indx, const int& altmax_lastyear_indx, ArrayD1 rootr, double& btran);

} // namespace ELM::soil_moist_stress

#include "soil_moist_stress_impl.hh"
