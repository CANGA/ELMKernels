
#pragma once

#include "compile_options.hh"
#include "data_types.hh"

namespace ELM {

  void kokkos_evaluate_conservation(ELMStateType& S,
                                    const double& dtime);

}
