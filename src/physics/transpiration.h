
#pragma once

#include "elm_constants.h"

#include "compile_options.hh"

namespace ELM::trans {

template <typename ArrayD1>
ACCELERATE
void transpiration(const bool& veg_active, 
                   const double& qflx_tran_veg,
                   const ArrayD1 rootr,
                   ArrayD1 qflx_rootsoi);


} // namespace ELM::trans

#include "transpiration_impl.hh"
