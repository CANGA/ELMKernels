
#include "SnowSNICAR.h"
#include "ELMConstants.h"
#include "LandType.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>

/*    

tests SnowSNICAR kernels: 
SurfRadZeroFluxes()
SurfRadAbsorbed()
SurfRadLayers()
SurfRadReflected()

the following data comes from the files SnowSNICAR_IN.txt and SnowSNICAR_OUT.txt located in test/data


urbpoi
flg_slr_in
coszen
h2osno
snl
h2osoi_liq
h2osoi_ice
snw_rds

snl_top
snl_btm
flx_abs_lcl
flx_abs
flg_nosnl
h2osoi_ice_lcl
h2osoi_liq_lcl
snw_rds_lcl
rds_bcint_lcl
rds_bcext_lcl
mu_not
flx_slrd_lcl
flx_slri_lcl

ss_alb_oc1
asm_prm_oc1
ext_cff_mss_oc1
ss_alb_oc2
asm_prm_oc2
ext_cff_mss_oc2
ss_alb_dst1
asm_prm_dst1
ext_cff_mss_dst1
ss_alb_dst2
asm_prm_dst2
ext_cff_mss_dst2
ss_alb_dst3
asm_prm_dst3
ext_cff_mss_dst3
ss_alb_dst4
asm_prm_dst4
ext_cff_mss_dst4

ss_alb_snw_drc
asm_prm_snw_drc
ext_cff_mss_snw_drc
ss_alb_snw_dfs
asm_prm_snw_dfs
ext_cff_mss_snw_dfs

ss_alb_bc1
asm_prm_bc1
ext_cff_mss_bc1
ss_alb_bc2
asm_prm_bc2
ext_cff_mss_bc2


mss_cnc_aer_in
bcenh
g_star
omega_star
tau_star

*/

using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;
using ArrayD3 = ELM::Array<double, 3>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(name, D0, D1, D2); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files 
  const std::string data_dir("/Users/80x/Software/elm_kernels/test/data/");
  const std::string input_file = data_dir + "SnowSNICAR_IN.txt";
  const std::string output_file = data_dir + "SnowSNICAR_OUT.txt";

  // hardwired
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;


// ELM state variables



auto urbpoi = create<ArrayI1>("urbpoi", n_grid_cells);
auto flg_slr_in = create<ArrayI1>("flg_slr_in", n_grid_cells);
auto coszen = create<ArrayD1>("coszen", n_grid_cells);
auto h2osno = create<ArrayD1>("h2osno", n_grid_cells);
auto snl = create<ArrayI1>("snl", n_grid_cells);
auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto snw_rds = create<ArrayD2>("snw_rds", n_grid_cells, ELM::nlevsno);


auto snl_top = create<ArrayI1>("snl_top", n_grid_cells);
auto snl_btm = create<ArrayI1>("snl_btm", n_grid_cells);
auto flx_abs_lcl = create<ArrayD3>("flx_abs_lcl", n_grid_cells, ELM::nlevsno+1, ELM::SNICAR::numrad_snw);
auto flx_abs = create<ArrayD3>("flx_abs", n_grid_cells, ELM::nlevsno+1, ELM::numrad);
auto flg_nosnl = create<ArrayI1>("flg_nosnl", n_grid_cells);
auto h2osoi_ice_lcl = create<ArrayD2>("h2osoi_ice_lcl", n_grid_cells, ELM::nlevsno);
auto h2osoi_liq_lcl = create<ArrayD2>("h2osoi_liq_lcl", n_grid_cells, ELM::nlevsno);
auto snw_rds_lcl = create<ArrayD2>("snw_rds_lcl", n_grid_cells, ELM::nlevsno);
auto rds_bcint_lcl = create<ArrayD2>("rds_bcint_lcl", n_grid_cells, ELM::nlevsno);
auto rds_bcext_lcl = create<ArrayD2>("rds_bcext_lcl", n_grid_cells, ELM::nlevsno);
auto mu_not = create<ArrayD1>("mu_not", n_grid_cells);
auto flx_slrd_lcl = create<ArrayD2>("flx_slrd_lcl", n_grid_cells, ELM::SNICAR::numrad_snw);
auto flx_slri_lcl = create<ArrayD2>("flx_slri_lcl", n_grid_cells, ELM::SNICAR::numrad_snw);




auto ss_alb_oc1 = create<ArrayD1>("ss_alb_oc1", ELM::SNICAR::numrad_snw);
auto asm_prm_oc1 = create<ArrayD1>("asm_prm_oc1", ELM::SNICAR::numrad_snw);
auto ext_cff_mss_oc1 = create<ArrayD1>("ext_cff_mss_oc1", ELM::SNICAR::numrad_snw);
auto ss_alb_oc2 = create<ArrayD1>("ss_alb_oc2", ELM::SNICAR::numrad_snw);
auto asm_prm_oc2 = create<ArrayD1>("asm_prm_oc2", ELM::SNICAR::numrad_snw);
auto ext_cff_mss_oc2 = create<ArrayD1>("ext_cff_mss_oc2", ELM::SNICAR::numrad_snw);
auto ss_alb_dst1 = create<ArrayD1>("ss_alb_dst1", ELM::SNICAR::numrad_snw);
auto asm_prm_dst1 = create<ArrayD1>("asm_prm_dst1", ELM::SNICAR::numrad_snw);
auto ext_cff_mss_dst1 = create<ArrayD1>("ext_cff_mss_dst1", ELM::SNICAR::numrad_snw);
auto ss_alb_dst2 = create<ArrayD1>("ss_alb_dst2", ELM::SNICAR::numrad_snw);

auto asm_prm_dst2 = create<ArrayD1>("asm_prm_dst2", ELM::SNICAR::numrad_snw);
auto ext_cff_mss_dst2 = create<ArrayD1>("ext_cff_mss_dst2", ELM::SNICAR::numrad_snw);
auto ss_alb_dst3 = create<ArrayD1>("ss_alb_dst3", ELM::SNICAR::numrad_snw);
auto asm_prm_dst3 = create<ArrayD1>("asm_prm_dst3", ELM::SNICAR::numrad_snw);
auto ext_cff_mss_dst3 = create<ArrayD1>("ext_cff_mss_dst3", ELM::SNICAR::numrad_snw);
auto ss_alb_dst4 = create<ArrayD1>("ss_alb_dst4", ELM::SNICAR::numrad_snw);
auto asm_prm_dst4 = create<ArrayD1>("asm_prm_dst4", ELM::SNICAR::numrad_snw);
auto ext_cff_mss_dst4 = create<ArrayD1>("ext_cff_mss_dst4", ELM::SNICAR::numrad_snw);




auto ss_alb_snw_drc = create<ArrayD2>("ss_alb_snw_drc", ELM::SNICAR::numrad_snw, ELM::SNICAR::idx_Mie_snw_mx);
auto asm_prm_snw_drc = create<ArrayD2>("asm_prm_snw_drc", ELM::SNICAR::numrad_snw, ELM::SNICAR::idx_Mie_snw_mx);
auto ext_cff_mss_snw_drc = create<ArrayD2>("ext_cff_mss_snw_drc", ELM::SNICAR::numrad_snw, ELM::SNICAR::idx_Mie_snw_mx);
auto ss_alb_snw_dfs = create<ArrayD2>("ss_alb_snw_dfs", ELM::SNICAR::numrad_snw, ELM::SNICAR::idx_Mie_snw_mx);
auto asm_prm_snw_dfs = create<ArrayD2>("asm_prm_snw_dfs", ELM::SNICAR::numrad_snw, ELM::SNICAR::idx_Mie_snw_mx);
auto ext_cff_mss_snw_dfs = create<ArrayD2>("ext_cff_mss_snw_dfs", ELM::SNICAR::numrad_snw, ELM::SNICAR::idx_Mie_snw_mx);

auto ss_alb_bc1 = create<ArrayD2>("ss_alb_bc1", E::SNICAR::idx_bc_nclrds_max, ELM::SNICAR::numrad_snw);
auto asm_prm_bc1 = create<ArrayD2>("asm_prm_bc1", ELM::SNICAR::idx_bc_nclrds_max, ELM::SNICAR::numrad_snw);
auto ext_cff_mss_bc1 = create<ArrayD2>("ext_cff_mss_bc1", ELM::SNICAR::idx_bc_nclrds_max, ELM::SNICAR::numrad_snw);
auto ss_alb_bc2 = create<ArrayD2>("ss_alb_bc2", E::SNICAR::idx_bc_nclrds_max, ELM::SNICAR::numrad_snw);
auto asm_prm_bc2 = create<ArrayD2>("asm_prm_bc2", ELM::SNICAR::idx_bc_nclrds_max, ELM::SNICAR::numrad_snw);
auto ext_cff_mss_bc2 = create<ArrayD2>("ext_cff_mss_bc2", ELM::SNICAR::idx_bc_nclrds_max, ELM::SNICAR::numrad_snw);


mss_cnc_aer_in[nlevsno][sno_nbr_aer]
bcenh[idx_bcint_icerds_max][idx_bc_nclrds_max][numrad_snw]
g_star[numrad_snw][nlevsno]
omega_star[numrad_snw][nlevsno]
tau_star[numrad_snw][nlevsno]

auto mss_cnc_aer_in
auto bcenh
auto g_star
auto omega_star
auto tau_star

mss_cnc_aer_in
bcenh
g_star
omega_star
tau_star

}






