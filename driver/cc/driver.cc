#include "clm_constants.h"
#include "landtype.h"
#include "BareGroundFluxes.h"
#include "array.hh"

int main(int argc, char** argv)
{
  const int ncells = 3;
  const int ntimes = 2;

  // instantiate data
  ELM::LandType land;
  ELM::Array<int,1> frac_vec_nosno_in(ncells);
  ELM::Array<double,1> forc_u(ncells);
  ELM::Array<double,1> forc_v(ncells);
  ELM::Array<double,1> forc_q_in(ncells);
  ELM::Array<double,1> forc_th_in(ncells);
  ELM::Array<double,1> forc_hgt_u_patch_in(ncells);
  ELM::Array<double,1> thm_in(ncells);
  ELM::Array<double,1> thv_in(ncells);
  ELM::Array<double,1> t_grnd(ncells);
  ELM::Array<double,1> qg(ncells);
  ELM::Array<double,1> z0mg_in(ncells);
  ELM::Array<double,1> dlrad(ncells);
  ELM::Array<double,1> ulrad(ncells);

  ELM::Array<double,1> forc_hgt_t_patch(ncells);
  ELM::Array<double,1> forc_hgt_q_patch(ncells);
  ELM::Array<double,1> z0mg(ncells);
  ELM::Array<double,1> zii(ncells);
  ELM::Array<double,1> beta(ncells);
  ELM::Array<double,1> z0hg(ncells);
  ELM::Array<double,1> z0qg(ncells);

  ELM::Array<int,1> snl(ncells, 1); // fix me!
  ELM::Array<double,1> forc_rho(ncells);
  ELM::Array<double,1> soilbeta(ncells);
  ELM::Array<double,1> dqgdT(ncells);
  ELM::Array<double,1> htvp(ncells);
  ELM::Array<double,1> t_h2osfc(ncells);
  ELM::Array<double,1> qg_snow(ncells);
  ELM::Array<double,1> qg_soil(ncells);
  ELM::Array<double,1> qg_h2osfc(ncells);

  ELM::Array<double,2> t_soisno(ncells, ELM::nlevsno+1); // is this correct?  This is required by BareGroundFluxes_impl.hh:108
  ELM::Array<double,1> forc_pbot(ncells);

  ELM::Array<double,1> cgrnds(ncells);
  ELM::Array<double,1> cgrndl(ncells);
  ELM::Array<double,1> cgrnd(ncells);
  ELM::Array<double,1> eflx_sh_grnd(ncells);
  ELM::Array<double,1> eflx_sh_tot(ncells);
  ELM::Array<double,1> eflx_sh_snow(ncells);
  ELM::Array<double,1> eflx_sh_soil(ncells);
  ELM::Array<double,1> eflx_sh_h2osfc(ncells);
  ELM::Array<double,1> qflx_evap_soi(ncells);
  ELM::Array<double,1> qflx_evap_tot(ncells);
  ELM::Array<double,1> qflx_ev_snow(ncells);
  ELM::Array<double,1> qflx_ev_soil(ncells);
  ELM::Array<double,1> qflx_ev_h2osfc(ncells);
  ELM::Array<double,1> t_ref2m(ncells);
  ELM::Array<double,1> t_ref2m_r(ncells);
  ELM::Array<double,1> q_ref2m(ncells);
  ELM::Array<double,1> rh_ref2m(ncells);
  ELM::Array<double,1> rh_ref2m_r(ncells);

  // initialize
  ELM::Array<ELM::BareGroundFluxes, 1> bg_fluxes(ncells);
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
  }

  // iterate in time
  for (int time=0; time!=ntimes; ++time) {
    for (int c=0; c!=ncells; ++c) {
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

