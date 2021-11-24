#include "bareground_fluxes.h"
#include "elm_constants.h"
#include "landtype.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>

/*

tests BareGroundFluxes kernels: 
InitializeFlux_BG()
StabilityIteration_BG()
ComputeFlux_BG()

implicitly tests FrictionVelocity kernels

the following data comes from the files BareGroundFluxes_IN.txt and BareGroundFluxes_OUT.txt located in test/data
int frac_veg_nosno
int snl
double forc_u
double forc_v
double forc_q
double forc_th
double thm
double thv
double t_grnd
double qg
double z0mg
double dlrad
double ulrad
double forc_hgt_t_patch
double forc_hgt_u_patch
double forc_hgt_q_patch
double z0hg
double z0qg
double forc_rho
double soilbeta
double dqgdT
double htvp
double t_h2osfc
double qg_snow
double qg_soil
double qg_h2osfc
double forc_pbot
double cgrnds
double cgrndl
double cgrnd
double eflx_sh_grnd
double eflx_sh_tot
double eflx_sh_snow
double eflx_sh_soil
double eflx_sh_h2osfc
double qflx_evap_soi
double qflx_evap_tot
double qflx_ev_snow
double qflx_ev_soil
double qflx_ev_h2osfc
double t_ref2m
double t_ref2m_r
double q_ref2m
double rh_ref2m
double rh_ref2m_r
ArrayD1 t_soisno[nlevsno+nlevgrnd]

*/


using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files 
  const std::string data_dir("/Users/80x/Software/elm_kernels/test/data/");
  const std::string input_file = data_dir + "BareGroundFluxes_IN.txt";
  const std::string output_file = data_dir + "BareGroundFluxes_OUT.txt";

  // hardwired
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;


  // temporary data to pass between functions
  double zldis;   // reference height "minus" zero displacement height [m]
  double displa;  // displacement height [m]
  double dth;     // diff of virtual temp. between ref. height and surface
  double dqh;     // diff of humidity between ref. height and surface
  double obu;     // Monin-Obukhov length (m)
  double ur;      // wind speed at reference height [m/s]
  double um;      // wind speed including the stablity effect [m/s]
  double temp1;   // relation for potential temperature profile
  double temp2;   // relation for specific humidity profile
  double temp12m; // relation for potential temperature profile applied at 2-m
  double temp22m; // relation for specific humidity profile applied at 2-m
  double ustar;   // friction velocity [m/s]



  // ELM state variables
  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto snl = create<ArrayI1>("snl", n_grid_cells);
  
  auto forc_u = create<ArrayD1>("forc_u", n_grid_cells);
  auto forc_v = create<ArrayD1>("forc_v", n_grid_cells);
  auto forc_q = create<ArrayD1>("forc_q", n_grid_cells);
  auto forc_th = create<ArrayD1>("forc_th", n_grid_cells);
  auto thm = create<ArrayD1>("thm", n_grid_cells);
  auto thv = create<ArrayD1>("thv", n_grid_cells);
  auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
  auto qg = create<ArrayD1>("qg", n_grid_cells);
  auto z0mg = create<ArrayD1>("z0mg", n_grid_cells);
  auto dlrad = create<ArrayD1>("dlrad", n_grid_cells);
  auto ulrad = create<ArrayD1>("ulrad", n_grid_cells);
  auto forc_hgt_t_patch = create<ArrayD1>("forc_hgt_t_patch", n_grid_cells);
  auto forc_hgt_u_patch = create<ArrayD1>("forc_hgt_u_patch", n_grid_cells);
  auto forc_hgt_q_patch = create<ArrayD1>("forc_hgt_q_patch", n_grid_cells);
  auto z0hg = create<ArrayD1>("z0hg", n_grid_cells);
  auto z0qg = create<ArrayD1>("z0qg", n_grid_cells);
  auto forc_rho = create<ArrayD1>("forc_rho", n_grid_cells);
  auto soilbeta = create<ArrayD1>("soilbeta", n_grid_cells);
  auto dqgdT = create<ArrayD1>("dqgdT", n_grid_cells);
  auto htvp = create<ArrayD1>("htvp", n_grid_cells);
  auto t_h2osfc = create<ArrayD1>("t_h2osfc", n_grid_cells);
  auto qg_snow = create<ArrayD1>("qg_snow", n_grid_cells);
  auto qg_soil = create<ArrayD1>("qg_soil", n_grid_cells);
  auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", n_grid_cells);
  auto forc_pbot = create<ArrayD1>("forc_pbot", n_grid_cells);
  auto cgrnds = create<ArrayD1>("cgrnds", n_grid_cells);
  auto cgrndl = create<ArrayD1>("cgrndl", n_grid_cells);
  auto cgrnd = create<ArrayD1>("cgrnd", n_grid_cells);
  auto eflx_sh_grnd = create<ArrayD1>("eflx_sh_grnd", n_grid_cells);
  auto eflx_sh_tot = create<ArrayD1>("eflx_sh_tot", n_grid_cells);
  auto eflx_sh_snow = create<ArrayD1>("eflx_sh_snow", n_grid_cells);
  auto eflx_sh_soil = create<ArrayD1>("eflx_sh_soil", n_grid_cells);
  auto eflx_sh_h2osfc = create<ArrayD1>("eflx_sh_h2osfc", n_grid_cells);
  auto qflx_evap_soi = create<ArrayD1>("qflx_evap_soi", n_grid_cells);
  auto qflx_evap_tot = create<ArrayD1>("qflx_evap_tot", n_grid_cells);
  auto qflx_ev_snow = create<ArrayD1>("qflx_ev_snow", n_grid_cells);
  auto qflx_ev_soil = create<ArrayD1>("qflx_ev_soil", n_grid_cells);
  auto qflx_ev_h2osfc = create<ArrayD1>("qflx_ev_h2osfc", n_grid_cells);
  auto t_ref2m = create<ArrayD1>("t_ref2m", n_grid_cells);
  auto t_ref2m_r = create<ArrayD1>("t_ref2m_r", n_grid_cells);
  auto q_ref2m = create<ArrayD1>("q_ref2m", n_grid_cells);
  auto rh_ref2m = create<ArrayD1>("rh_ref2m", n_grid_cells);
  auto rh_ref2m_r = create<ArrayD1>("rh_ref2m_r", n_grid_cells);
  
  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);


  // input and output utility class objects
  ELM::IO::ELMtestinput in(input_file);
  ELM::IO::ELMtestinput out(output_file);

  for (std::size_t t = 1; t < 49; ++t) {

    // get input and output state for time t
    in.getState(t);
    out.getState(t);

    // parse input state and assign to variables
    //in.parseState(frac_veg_nosno); hardwired
    in.parseState(snl);
    in.parseState(forc_u);
    in.parseState(forc_v);
    in.parseState(forc_q);
    in.parseState(forc_th);
    in.parseState(thm);
    in.parseState(thv);
    in.parseState(t_grnd);
    in.parseState(qg);
    in.parseState(z0mg);
    in.parseState(dlrad);
    in.parseState(ulrad);
    in.parseState(forc_hgt_t_patch);
    in.parseState(forc_hgt_u_patch);
    in.parseState(forc_hgt_q_patch);
    in.parseState(z0hg);
    in.parseState(z0qg);
    in.parseState(forc_rho);
    in.parseState(soilbeta);
    in.parseState(dqgdT);
    in.parseState(htvp);
    in.parseState(t_h2osfc);
    in.parseState(qg_snow);
    in.parseState(qg_soil);
    in.parseState(qg_h2osfc);
    in.parseState(forc_pbot);
    in.parseState(cgrnds);
    in.parseState(cgrndl);
    in.parseState(cgrnd);
    in.parseState(eflx_sh_grnd);
    in.parseState(eflx_sh_tot);
    in.parseState(eflx_sh_snow);
    in.parseState(eflx_sh_soil);
    in.parseState(eflx_sh_h2osfc);
    in.parseState(qflx_evap_soi);
    in.parseState(qflx_evap_tot);
    in.parseState(qflx_ev_snow);
    in.parseState(qflx_ev_soil);
    in.parseState(qflx_ev_h2osfc);
    in.parseState(t_ref2m);
    in.parseState(t_ref2m_r);
    in.parseState(q_ref2m);
    in.parseState(rh_ref2m);
    in.parseState(rh_ref2m_r);
    in.parseState(t_soisno[idx]);
    
    assign(frac_veg_nosno, 0); // hardwire to make it run
    
    // call BareGroundFluxes kernels
    ELM::InitializeFlux_BG(Land, frac_veg_nosno[idx], forc_u[idx], forc_v[idx], forc_q[idx], forc_th[idx], forc_hgt_u_patch[idx], thm[idx],
                      thv[idx], t_grnd[idx], qg[idx], z0mg[idx], dlrad[idx], ulrad[idx], zldis, displa, dth, dqh, obu, ur, um);
    
    ELM::StabilityIteration_BG(Land, frac_veg_nosno[idx], forc_hgt_t_patch[idx], forc_hgt_u_patch[idx], forc_hgt_q_patch[idx],
                          z0mg[idx], zldis, displa, dth, dqh, ur, forc_q[idx], forc_th[idx], thv[idx], z0hg[idx],
                          z0qg[idx], obu, um, temp1, temp2, temp12m, temp22m, ustar);
    
    ELM::ComputeFlux_BG(Land, frac_veg_nosno[idx], snl[idx], forc_rho[idx], soilbeta[idx], dqgdT[idx], htvp[idx], t_h2osfc[idx],
                   qg_snow[idx], qg_soil[idx], qg_h2osfc[idx], t_soisno[idx], forc_pbot[idx], dth,
                   dqh, temp1, temp2, temp12m, temp22m, ustar, forc_q[idx], thm[idx], cgrnds[idx], cgrndl[idx], cgrnd[idx],
                   eflx_sh_grnd[idx], eflx_sh_tot[idx], eflx_sh_snow[idx], eflx_sh_soil[idx], eflx_sh_h2osfc[idx],
                   qflx_evap_soi[idx], qflx_evap_tot[idx], qflx_ev_snow[idx], qflx_ev_soil[idx], qflx_ev_h2osfc[idx], t_ref2m[idx],
                   t_ref2m_r[idx], q_ref2m[idx], rh_ref2m[idx], rh_ref2m_r[idx]);
    
    
    // compare kernel output to ELM output state
    //out.compareOutput(frac_veg_nosno); hardwired
    out.compareOutput(snl);
    out.compareOutput(forc_u);
    out.compareOutput(forc_v);
    out.compareOutput(forc_q);
    out.compareOutput(forc_th);
    out.compareOutput(thm);
    out.compareOutput(thv);
    out.compareOutput(t_grnd);
    out.compareOutput(qg);
    out.compareOutput(z0mg);
    out.compareOutput(dlrad);
    out.compareOutput(ulrad);
    out.compareOutput(forc_hgt_t_patch);
    out.compareOutput(forc_hgt_u_patch);
    out.compareOutput(forc_hgt_q_patch);
    out.compareOutput(z0hg);
    out.compareOutput(z0qg);
    out.compareOutput(forc_rho);
    out.compareOutput(soilbeta);
    out.compareOutput(dqgdT);
    out.compareOutput(htvp);
    out.compareOutput(t_h2osfc);
    out.compareOutput(qg_snow);
    out.compareOutput(qg_soil);
    out.compareOutput(qg_h2osfc);
    out.compareOutput(forc_pbot);
    out.compareOutput(cgrnds);
    out.compareOutput(cgrndl);
    out.compareOutput(cgrnd);
    out.compareOutput(eflx_sh_grnd);
    out.compareOutput(eflx_sh_tot);
    out.compareOutput(eflx_sh_snow);
    out.compareOutput(eflx_sh_soil);
    out.compareOutput(eflx_sh_h2osfc);
    out.compareOutput(qflx_evap_soi);
    out.compareOutput(qflx_evap_tot);
    out.compareOutput(qflx_ev_snow);
    out.compareOutput(qflx_ev_soil);
    out.compareOutput(qflx_ev_h2osfc);
    out.compareOutput(t_ref2m);
    out.compareOutput(t_ref2m_r);
    out.compareOutput(q_ref2m);
    out.compareOutput(rh_ref2m);
    out.compareOutput(rh_ref2m_r);
    out.compareOutput(t_soisno[idx]);
  }
  return 0;
}
