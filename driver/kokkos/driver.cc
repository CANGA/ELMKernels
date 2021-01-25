#include "clm_constants.h"
#include "landtype.h"
#include "BareGroundFluxes.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_functors.hh"


using ArrayD1 = Kokkos::View<double*>;
using ArrayI1 = Kokkos::View<int*>;
using ArrayD2 = Kokkos::View<double**>;

template<typename Array_t>
Array_t create(const std::string& name, int D0) {
  return Array_t(name, D0);
}

template<typename Array_t>
Array_t create(const std::string& name, int D0, int D1) {
  return Array_t(name, D0, D1);
}

template<typename Array_t>
Array_t create(const std::string& name, int D0, int D1, int D2) {
  return Array_t(name, D0, D1, D2);
}

template<class Array_t, typename Scalar_t>
void assign(Array_t& arr, Scalar_t val) {
  Kokkos::deep_copy(arr, val);
}





int main(int argc, char** argv)
{
  Kokkos::initialize( argc, argv);

  const int ncells = 300000;
  const int ntimes = 2000000;

  { // scope to make Kokkos happy
  // instantiate data
  ELM::LandType land;
  auto frac_vec_nosno_in = create<ArrayI1>("frac_vec_nosno_in", ncells);
  auto forc_u = create<ArrayD1>("forc_u", ncells);
  auto forc_v = create<ArrayD1>("forc_v", ncells);
  auto forc_q_in = create<ArrayD1>("forc_q_in", ncells);
  auto forc_th_in = create<ArrayD1>("forc_th_in", ncells);
  auto forc_hgt_u_patch_in = create<ArrayD1>("forc_hgt_u_patch_in", ncells);
  auto thm_in = create<ArrayD1>("thm_in", ncells);
  auto thv_in = create<ArrayD1>("thv_in", ncells);
  auto t_grnd = create<ArrayD1>("t_grnd", ncells);
  auto qg = create<ArrayD1>("qg", ncells);
  auto z0mg_in = create<ArrayD1>("z0mg_in", ncells);
  auto dlrad = create<ArrayD1>("dlrad", ncells);
  auto ulrad = create<ArrayD1>("ulrad", ncells);

  auto forc_hgt_t_patch = create<ArrayD1>("forc_hgt_t_patch", ncells);
  auto forc_hgt_q_patch = create<ArrayD1>("forc_hgt_q_patch", ncells);
  auto z0mg = create<ArrayD1>("z0mg", ncells);
  auto zii = create<ArrayD1>("zii", ncells);
  auto beta = create<ArrayD1>("beta", ncells);
  auto z0hg = create<ArrayD1>("z0hg", ncells);
  auto z0qg = create<ArrayD1>("z0qg", ncells);

  auto snl = create<ArrayI1>("snl", ncells); assign(snl, 1); // fix me!
  auto forc_rho = create<ArrayD1>("forc_rho", ncells);
  auto soilbeta = create<ArrayD1>("soilbeta", ncells);
  auto dqgdT = create<ArrayD1>("dqgdT", ncells);
  auto htvp = create<ArrayD1>("htvp", ncells);
  auto t_h2osfc = create<ArrayD1>("t_h2osfc", ncells);
  auto qg_snow = create<ArrayD1>("qg_snow", ncells);
  auto qg_soil = create<ArrayD1>("qg_soil", ncells);
  auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", ncells);

  auto t_soisno = create<ArrayD2>("t_soisno", ncells, ELM::nlevsno + ELM::nlevgrnd);
  auto forc_pbot = create<ArrayD1>("forc_pbot", ncells);

  auto cgrnds = create<ArrayD1>("cgrnds", ncells);
  auto cgrndl = create<ArrayD1>("cgrndl", ncells);
  auto cgrnd = create<ArrayD1>("cgrnd", ncells);
  auto eflx_sh_grnd = create<ArrayD1>("eflx_sh_grnd", ncells);
  auto eflx_sh_tot = create<ArrayD1>("eflx_sh_tot", ncells);
  auto eflx_sh_snow = create<ArrayD1>("eflx_sh_snow", ncells);
  auto eflx_sh_soil = create<ArrayD1>("eflx_sh_soil", ncells);
  auto eflx_sh_h2osfc = create<ArrayD1>("eflx_sh_h2osfc", ncells);
  auto qflx_evap_soi = create<ArrayD1>("qflx_evap_soi", ncells);
  auto qflx_evap_tot = create<ArrayD1>("qflx_evap_tot", ncells);
  auto qflx_ev_snow = create<ArrayD1>("qflx_ev_snow", ncells);
  auto qflx_ev_soil = create<ArrayD1>("qflx_ev_soil", ncells);
  auto qflx_ev_h2osfc = create<ArrayD1>("qflx_ev_h2osfc", ncells);
  auto t_ref2m = create<ArrayD1>("t_ref2m", ncells);
  auto t_ref2m_r = create<ArrayD1>("t_ref2m_r", ncells);
  auto q_ref2m = create<ArrayD1>("q_ref2m", ncells);
  auto rh_ref2m = create<ArrayD1>("rh_ref2m", ncells);
  auto rh_ref2m_r = create<ArrayD1>("rh_ref2m_r", ncells);

  // initialize
  Kokkos::View<ELM::BareGroundFluxes*> bg_fluxes("bare ground fluxes", ncells);

  BGFCaller bgf_caller(bg_fluxes, land, frac_vec_nosno_in, forc_u, forc_v, forc_q_in, forc_th_in, forc_hgt_u_patch_in, 
    thm_in, thv_in, t_grnd, qg, z0mg_in, dlrad, ulrad, forc_hgt_t_patch, forc_hgt_q_patch, z0mg, zii, beta, 
    z0hg, z0qg, snl, forc_rho, soilbeta, dqgdT, htvp, t_h2osfc, qg_snow, qg_soil, qg_h2osfc, t_soisno, forc_pbot, 
    cgrnds, cgrndl, cgrnd, eflx_sh_grnd, eflx_sh_tot, eflx_sh_snow, eflx_sh_soil, eflx_sh_h2osfc, qflx_evap_soi, 
    qflx_evap_tot, qflx_ev_snow, qflx_ev_soil, qflx_ev_h2osfc, t_ref2m, t_ref2m_r, q_ref2m, rh_ref2m, rh_ref2m_r);

  for (int time=0; time!=ntimes; ++time) {
   
   // parallel functor call
    Kokkos::parallel_for(ncells, bgf_caller);
    std::cout << "running, t = " << time << std::endl;
  }

  }
  Kokkos::finalize();
  return 0;
}

