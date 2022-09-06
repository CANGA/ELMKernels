
#include "invoke_kernel.hh"
#include "surface_fluxes.h"

#include "init_timestep_kokkos.hh"
#include "init_timestep.h"

void ELM::kokkos_init_timestep(ELMStateType& S)
{
  auto init_step_kernel = ELM_LAMBDA (const int& idx) {
    ELM::init_timestep(S.Land.lakpoi, S.veg_active(idx),
                     S.frac_veg_nosno_alb(idx),
                     S.snl(idx), S.h2osno(idx),
                     Kokkos::subview(S.h2osoi_ice, idx, Kokkos::ALL),
                     Kokkos::subview(S.h2osoi_liq, idx, Kokkos::ALL),
                     S.do_capsnow(idx),
                     S.frac_veg_nosno(idx),
                     Kokkos::subview(S.frac_iceold, idx, Kokkos::ALL));
  }; // end init_step lambda
  invoke_kernel(init_step_kernel, std::make_tuple(S.snl.extent(0)), "kokkos_init_timestep");
}
