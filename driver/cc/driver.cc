#include <string>
#include "clm_constants.h"
#include "landtype.h"
#include "BareGroundFluxes.h"
#include "array.hh"


using ArrayD1 = ELM::Array<double,1>;
using ArrayI1 = ELM::Array<int,1>;
using ArrayD2 = ELM::Array<double,2>;

template<class Array_t>
Array_t create(const std::string& name, int D0) {
  return Array_t(D0);
}

template<class Array_t>
Array_t create(const std::string& name, int D0, int D1) {
  return Array_t(D0, D1);
}

template<class Array_t>
Array_t create(const std::string& name, int D0, int D1, int D2) {
  return Array_t(D0, D1, D2);
}


template<class Array_t, typename Scalar_t>
void assign(Array_t& arr, Scalar_t val) {
  ELM::deep_copy(arr, val);
}

int main(int argc, char** argv)
{
  const int ncells = 3;
  const int ntimes = 2;


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

  auto t_soisno = create<ArrayD2>("t_soisno", ncells, ELM::nlevsno+1); // is this correct?  This is required by BareGroundFluxes_impl.hh:108
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
  auto bg_fluxes = create<ELM::Array<ELM::BareGroundFluxes,1> >("bg_fluxes", ncells);

  // iterate in time
  for (int time=0; time!=ntimes; ++time) {
    for (int c=0; c!=ncells; ++c) {
      bg_fluxes[c].InitializeFlux(land,
              frac_vec_nosno_in[c],
              forc_u[c],
              forc_v[c],
              forc_q_in[c],
              forc_th_in[c],
              forc_hgt_u_patch_in[c],
              thm_in[c],
              thv_in[c],
              t_grnd[c],
              qg[c],
              z0mg_in[c],
              dlrad[c],
              ulrad[c]);

      bg_fluxes[c].StabilityIteration(land,
              forc_hgt_t_patch[c],
              forc_hgt_q_patch[c],
              z0mg[c],
              zii[c],
              beta[c],
              z0hg[c],
              z0qg[c]);

      bg_fluxes[c].ComputeFlux(land,
              snl[c],
              forc_rho[c],
              soilbeta[c],
              dqgdT[c],
              htvp[c],
              t_h2osfc[c],
              qg_snow[c],
              qg_soil[c],
              qg_h2osfc[c],
              t_soisno[c],
              forc_pbot[c],
              cgrnds[c],
              cgrndl[c],
              cgrnd[c],
              eflx_sh_grnd[c],
              eflx_sh_tot[c],
              eflx_sh_snow[c],
              eflx_sh_soil[c],
              eflx_sh_h2osfc[c],
              qflx_evap_soi[c],
              qflx_evap_tot[c],
              qflx_ev_snow[c],
              qflx_ev_soil[c],
              qflx_ev_h2osfc[c],
              t_ref2m[c],
              t_ref2m_r[c],
              q_ref2m[c],
              rh_ref2m[c],
              rh_ref2m_r[c]);
    }
  }

  return 0;
}

