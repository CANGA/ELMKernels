/*! \file InitSnowLayers.hh
\brief Function derived from initVerticalMod.F90
*/
#pragma once

namespace ELM {

/*! Choose number of snow layers (max 5) based on depth of snow. Initialize snow layer thickness, depth, interface
depth. CLM uses negative indexing for the snow model -- snl is the negative of the number of snow layers - i=0 is layer
next to soil, i = snl+1 is top layer. Here we use a positive snl -- i=0 corresponds to snl = nlevsno & i = nlevsno-1 is
the layer adjacent to the ground surface & i = [nlevsno-snl] is the top active snow layer.

\param[in]  snow_depth           [double] snow height (m)
\param[in]  lakpoi               [bool] true => landunit is a lake point
\param[out[ snl                  [integer] number of snow layers
\param[out[ dz[nlevgrnd+nlevsno] [double]  layer thickness (m)
\param[out[ z[nlevgrnd+nlevsno]  [double]  layer depth (m)
\param[out[ zi[nlevgrnd+nlevsno] [double]  interface level below a "z" level (m
*/
template <class ArrayD1>
void InitSnowLayers(const double &snow_depth, const bool &lakpoi, int &snl, ArrayD1 dz, ArrayD1 z, ArrayD1 zi);

} // namespace ELM

#include <InitSnowLayers_impl.hh>
