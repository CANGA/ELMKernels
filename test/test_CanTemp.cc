
#include "canopy_temperature.h"
#include "elm_constants.h"
#include "land_data.h"
#include "read_test_input.hh"
#include "array.hh"

#include <iostream>
#include <string>

/*

tests canopy_temperature kernels: 
old_ground_temp()
ground_temp()
calc_soilalpha()
calc_soilbeta()
humidities()
ground_properties()
forcing_height()
init_energy_fluxes()

implicitly tests QSat() and calc_soilevap_stress()

the following data comes from the files CanopyTemperature_IN.txt and CanopyTemperature_OUT.txt located in test/data

veg_active
snl
frac_veg_nosno
t_h2osfc
t_h2osfc_bef
frac_sno_eff
frac_h2osfc
t_grnd
frac_sno
smpmin
soilalpha
soilalpha_u
soilbeta
forc_q
forc_pbot
qg_snow
qg_soil
qg
qg_h2osfc
dqgdT
forc_th
elai
esai
htop
emg
emv
htvp
z0mg
z0hg
z0qg
z0mv
z0hv
z0qv
thv
z0m
displa
forc_hgt_u
forc_hgt_t
forc_hgt_q
z_0_town
z_d_town
forc_t
forc_hgt_u_patch
forc_hgt_t_patch
forc_hgt_q_patch
thm
eflx_sh_tot
eflx_sh_tot_u
eflx_sh_tot_r
eflx_lh_tot
eflx_lh_tot_u
eflx_lh_tot_r
eflx_sh_veg
qflx_evap_tot
qflx_evap_veg
qflx_tran_veg
t_soisno[nlevgrnd+nlevsno]
tssbef[nlevgrnd+nlevsno]
h2osoi_liq[nlevgrnd+nlevsno]
h2osoi_ice[nlevgrnd+nlevsno]
dz[nlevgrnd+nlevsno]
watsat[nlevgrnd]
sucsat[nlevgrnd]
bsw[nlevgrnd]
watdry[nlevgrnd]
watopt[nlevgrnd]
rootfr_road_perv[nlevgrnd]
rootr_road_perv[nlevgrnd]
watfc[nlevgrnd]
displar[numpft]
z0mr[numpft]
*/


using ArrayI1 = ELM::Array<int, 1>;
using ArrayB1 = ELM::Array<bool, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files 
  const std::string data_dir("/Users/80x/Software/elm_kernels/test/data/");
  const std::string input_file = data_dir + "CanopyTemperature_IN.txt";
  const std::string output_file = data_dir + "CanopyTemperature_OUT.txt";

  // hardwired
  ELM::LandType Land;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;

  // temporary data to pass between functions
  double qred; // soil surface relative humidity
  double hr;   // relative humidity

  // ELM state variables
  auto veg_active = create<ArrayB1>("veg_active", n_grid_cells);
  auto snl = create<ArrayI1>("snl", n_grid_cells);
  auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
  auto t_h2osfc = create<ArrayD1>("t_h2osfc", n_grid_cells);
  auto t_h2osfc_bef = create<ArrayD1>("t_h2osfc_bef", n_grid_cells);
  auto frac_sno_eff = create<ArrayD1>("frac_sno_eff", n_grid_cells);
  auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
  auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
  auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);
  auto smpmin = create<ArrayD1>("smpmin", n_grid_cells);
  auto soilalpha = create<ArrayD1>("soilalpha", n_grid_cells);
  auto soilalpha_u = create<ArrayD1>("soilalpha_u", n_grid_cells);
  auto soilbeta = create<ArrayD1>("soilbeta", n_grid_cells);
  auto forc_q = create<ArrayD1>("forc_q", n_grid_cells);
  auto forc_pbot = create<ArrayD1>("forc_pbot", n_grid_cells);
  auto qg_snow = create<ArrayD1>("qg_snow", n_grid_cells);
  auto qg_soil = create<ArrayD1>("qg_soil", n_grid_cells);
  auto qg = create<ArrayD1>("qg", n_grid_cells);
  auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", n_grid_cells);
  auto dqgdT = create<ArrayD1>("dqgdT", n_grid_cells);
  auto forc_th = create<ArrayD1>("forc_th", n_grid_cells);
  auto elai = create<ArrayD1>("elai", n_grid_cells);
  auto esai = create<ArrayD1>("esai", n_grid_cells);
  auto htop = create<ArrayD1>("htop", n_grid_cells);
  auto emg = create<ArrayD1>("emg", n_grid_cells);
  auto emv = create<ArrayD1>("emv", n_grid_cells);
  auto htvp = create<ArrayD1>("htvp", n_grid_cells);
  auto z0mg = create<ArrayD1>("z0mg", n_grid_cells);
  auto z0hg = create<ArrayD1>("z0hg", n_grid_cells);
  auto z0qg = create<ArrayD1>("z0qg", n_grid_cells);
  auto z0mv = create<ArrayD1>("z0mv", n_grid_cells);
  auto z0hv = create<ArrayD1>("z0hv", n_grid_cells);
  auto z0qv = create<ArrayD1>("z0qv", n_grid_cells);
  auto thv = create<ArrayD1>("thv", n_grid_cells);
  auto z0m = create<ArrayD1>("z0m", n_grid_cells);
  auto displa = create<ArrayD1>("displa", n_grid_cells);
  auto forc_hgt_u = create<ArrayD1>("forc_hgt_u", n_grid_cells);
  auto forc_hgt_t = create<ArrayD1>("forc_hgt_t", n_grid_cells);
  auto forc_hgt_q = create<ArrayD1>("forc_hgt_q", n_grid_cells);
  auto z_0_town = create<ArrayD1>("z_0_town", n_grid_cells);
  auto z_d_town = create<ArrayD1>("z_d_town", n_grid_cells);
  auto forc_t = create<ArrayD1>("forc_t", n_grid_cells);
  auto forc_hgt_u_patch = create<ArrayD1>("forc_hgt_u_patch", n_grid_cells);
  auto forc_hgt_t_patch = create<ArrayD1>("forc_hgt_t_patch", n_grid_cells);
  auto forc_hgt_q_patch = create<ArrayD1>("forc_hgt_q_patch", n_grid_cells);
  auto thm = create<ArrayD1>("thm", n_grid_cells);
  auto eflx_sh_tot = create<ArrayD1>("eflx_sh_tot", n_grid_cells);
  auto eflx_sh_tot_u = create<ArrayD1>("eflx_sh_tot_u", n_grid_cells);
  auto eflx_sh_tot_r = create<ArrayD1>("eflx_sh_tot_r", n_grid_cells);
  auto eflx_lh_tot = create<ArrayD1>("eflx_lh_tot", n_grid_cells);
  auto eflx_lh_tot_u = create<ArrayD1>("eflx_lh_tot_u", n_grid_cells);
  auto eflx_lh_tot_r = create<ArrayD1>("eflx_lh_tot_r", n_grid_cells);
  auto eflx_sh_veg = create<ArrayD1>("eflx_sh_veg", n_grid_cells);
  auto qflx_evap_tot = create<ArrayD1>("qflx_evap_tot", n_grid_cells);
  auto qflx_evap_veg = create<ArrayD1>("qflx_evap_veg", n_grid_cells);
  auto qflx_tran_veg = create<ArrayD1>("qflx_tran_veg", n_grid_cells);

  auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevgrnd+ELM::nlevsno);
  auto tssbef = create<ArrayD2>("tssbef", n_grid_cells, ELM::nlevgrnd+ELM::nlevsno);
  auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevgrnd+ELM::nlevsno);
  auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevgrnd+ELM::nlevsno);
  auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevgrnd+ELM::nlevsno);
  auto watsat = create<ArrayD2>("watsat", n_grid_cells, ELM::nlevgrnd);
  auto sucsat = create<ArrayD2>("sucsat", n_grid_cells, ELM::nlevgrnd);
  auto bsw = create<ArrayD2>("bsw", n_grid_cells, ELM::nlevgrnd);
  auto watdry = create<ArrayD2>("watdry", n_grid_cells, ELM::nlevgrnd);
  auto watopt = create<ArrayD2>("watopt", n_grid_cells, ELM::nlevgrnd);
  auto rootfr_road_perv = create<ArrayD2>("rootfr_road_perv", n_grid_cells, ELM::nlevgrnd);
  auto rootr_road_perv = create<ArrayD2>("rootr_road_perv", n_grid_cells, ELM::nlevgrnd);
  auto watfc = create<ArrayD2>("watfc", n_grid_cells, ELM::nlevgrnd);
  auto displar = create<ArrayD2>("displar", n_grid_cells, ELM::numpft);
  auto z0mr = create<ArrayD2>("z0mr", n_grid_cells, ELM::numpft);


  // input and output utility class objects
  ELM::IO::ELMtestinput in(input_file);
  ELM::IO::ELMtestinput out(output_file);

  for (std::size_t t = 1; t < 49; ++t) {

    // get input and output state for time t
    in.getState(t);
    out.getState(t);

    // parse input state and assign to variables
    in.parseState(veg_active);
    in.parseState(snl);
    in.parseState(frac_veg_nosno);
    in.parseState(t_h2osfc);
    in.parseState(t_h2osfc_bef);
    in.parseState(frac_sno_eff);
    in.parseState(frac_h2osfc);
    in.parseState(t_grnd);
    in.parseState(frac_sno);
    in.parseState(smpmin);
    in.parseState(soilalpha);
    in.parseState(soilalpha_u);
    in.parseState(soilbeta);
    in.parseState(forc_q);
    in.parseState(forc_pbot);
    in.parseState(qg_snow);
    in.parseState(qg_soil);
    in.parseState(qg);
    in.parseState(qg_h2osfc);
    in.parseState(dqgdT);
    in.parseState(forc_th);
    in.parseState(elai);
    in.parseState(esai);
    in.parseState(htop);
    in.parseState(emg);
    in.parseState(emv);
    in.parseState(htvp);
    in.parseState(z0mg);
    in.parseState(z0hg);
    in.parseState(z0qg);
    in.parseState(z0mv);
    in.parseState(z0hv);
    in.parseState(z0qv);
    in.parseState(thv);
    in.parseState(z0m);
    in.parseState(displa);
    in.parseState(forc_hgt_u);
    in.parseState(forc_hgt_t);
    in.parseState(forc_hgt_q);
    //in.parseState(z_0_town);
    //in.parseState(z_d_town);
    in.parseState(forc_t);
    in.parseState(forc_hgt_u_patch);
    in.parseState(forc_hgt_t_patch);
    in.parseState(forc_hgt_q_patch);
    in.parseState(thm);
    in.parseState(eflx_sh_tot);
    in.parseState(eflx_sh_tot_u);
    in.parseState(eflx_sh_tot_r);
    in.parseState(eflx_lh_tot);
    in.parseState(eflx_lh_tot_u);
    in.parseState(eflx_lh_tot_r);
    in.parseState(eflx_sh_veg);
    in.parseState(qflx_evap_tot);
    in.parseState(qflx_evap_veg);
    in.parseState(qflx_tran_veg);
    in.parseState(t_soisno[idx]);
    in.parseState(tssbef[idx]);
    in.parseState(h2osoi_liq[idx]);
    in.parseState(h2osoi_ice[idx]);
    in.parseState(dz[idx]);
    in.parseState(watsat[idx]);
    in.parseState(sucsat[idx]);
    in.parseState(bsw[idx]);
    in.parseState(watdry[idx]);
    in.parseState(watopt[idx]);
    //in.parseState(rootfr_road_perv[idx]);
    //in.parseState(rootr_road_perv[idx]);
    in.parseState(watfc[idx]);
    in.parseState(displar[idx]);
    in.parseState(z0mr[idx]);

    // call CanopyTemperature kernels
    ELM::canopy_temperature::old_ground_temp(Land, t_h2osfc[idx], t_soisno[idx], t_h2osfc_bef[idx], tssbef[idx]);

    ELM::canopy_temperature::ground_temp(Land, snl[idx], frac_sno_eff[idx], frac_h2osfc[idx], t_h2osfc[idx], t_soisno[idx],
                             t_grnd[idx]);

    ELM::canopy_temperature::calc_soilalpha(Land, frac_sno[idx], frac_h2osfc[idx], smpmin[idx], h2osoi_liq[idx], h2osoi_ice[idx],
                            dz[idx], t_soisno[idx], watsat[idx], sucsat[idx], bsw[idx], watdry[idx], watopt[idx],
                            rootfr_road_perv[idx], rootr_road_perv[idx], qred, hr, soilalpha[idx], soilalpha_u[idx]);

    ELM::canopy_temperature::calc_soilbeta(Land, frac_sno[idx], frac_h2osfc[idx], watsat[idx], watfc[idx], h2osoi_liq[idx],
                           h2osoi_ice[idx], dz[idx], soilbeta[idx]);

    ELM::canopy_temperature::humidities(Land, snl[idx], forc_q[idx], forc_pbot[idx], t_h2osfc[idx], t_grnd[idx], frac_sno[idx],
                             frac_sno_eff[idx], frac_h2osfc[idx], qred, hr, t_soisno[idx], qg_snow[idx], qg_soil[idx],
                             qg[idx], qg_h2osfc[idx], dqgdT[idx]);

    ELM::canopy_temperature::ground_properties(Land, snl[idx], frac_sno[idx], forc_th[idx], forc_q[idx], elai[idx], esai[idx], htop[idx],
                           displar[idx], z0mr[idx],
                           h2osoi_liq[idx], h2osoi_ice[idx],
                           emg[idx], emv[idx], htvp[idx], z0mg[idx], z0hg[idx], z0qg[idx], z0mv[idx], z0hv[idx], z0qv[idx],
                           thv[idx], z0m[idx], displa[idx]);

    ELM::canopy_temperature::forcing_height(Land, veg_active[idx], frac_veg_nosno[idx], forc_hgt_u[idx], forc_hgt_t[idx], forc_hgt_q[idx],
                                 z0m[idx], z0mg[idx], z_0_town[idx], z_d_town[idx], forc_t[idx], displa[idx], forc_hgt_u_patch[idx],
                                 forc_hgt_t_patch[idx], forc_hgt_q_patch[idx], thm[idx]);

    ELM::canopy_temperature::init_energy_fluxes(Land, eflx_sh_tot[idx], eflx_sh_tot_u[idx], eflx_sh_tot_r[idx], eflx_lh_tot[idx],
                                eflx_lh_tot_u[idx], eflx_lh_tot_r[idx], eflx_sh_veg[idx], qflx_evap_tot[idx], qflx_evap_veg[idx],
                                qflx_tran_veg[idx]);



    // compare kernel output to ELM output state
    out.compareOutput(veg_active);
    out.compareOutput(snl);
    out.compareOutput(frac_veg_nosno);
    out.compareOutput(t_h2osfc);
    out.compareOutput(t_h2osfc_bef);
    out.compareOutput(frac_sno_eff);
    out.compareOutput(frac_h2osfc);
    out.compareOutput(t_grnd);
    out.compareOutput(frac_sno);
    out.compareOutput(smpmin);
    out.compareOutput(soilalpha);
    out.compareOutput(soilalpha_u);
    out.compareOutput(soilbeta);
    out.compareOutput(forc_q);
    out.compareOutput(forc_pbot);
    out.compareOutput(qg_snow);
    out.compareOutput(qg_soil);
    out.compareOutput(qg);
    out.compareOutput(qg_h2osfc);
    out.compareOutput(dqgdT);
    out.compareOutput(forc_th);
    out.compareOutput(elai);
    out.compareOutput(esai);
    out.compareOutput(htop);
    out.compareOutput(emg);
    out.compareOutput(emv);
    out.compareOutput(htvp);
    out.compareOutput(z0mg);
    out.compareOutput(z0hg);
    out.compareOutput(z0qg);
    out.compareOutput(z0mv);
    out.compareOutput(z0hv);
    out.compareOutput(z0qv);
    out.compareOutput(thv);
    out.compareOutput(z0m);
    out.compareOutput(displa);
    out.compareOutput(forc_hgt_u);
    out.compareOutput(forc_hgt_t);
    out.compareOutput(forc_hgt_q);
    //out.compareOutput(z_0_town);
    //out.compareOutput(z_d_town);
    out.compareOutput(forc_t);
    out.compareOutput(forc_hgt_u_patch);
    out.compareOutput(forc_hgt_t_patch);
    out.compareOutput(forc_hgt_q_patch);
    out.compareOutput(thm);
    out.compareOutput(eflx_sh_tot);
    out.compareOutput(eflx_sh_tot_u);
    out.compareOutput(eflx_sh_tot_r);
    out.compareOutput(eflx_lh_tot);
    out.compareOutput(eflx_lh_tot_u);
    out.compareOutput(eflx_lh_tot_r);
    out.compareOutput(eflx_sh_veg);
    out.compareOutput(qflx_evap_tot);
    out.compareOutput(qflx_evap_veg);
    out.compareOutput(qflx_tran_veg);
    out.compareOutput(t_soisno[idx]);
    out.compareOutput(tssbef[idx]);
    out.compareOutput(h2osoi_liq[idx]);
    out.compareOutput(h2osoi_ice[idx]);
    out.compareOutput(dz[idx]);
    out.compareOutput(watsat[idx]);
    out.compareOutput(sucsat[idx]);
    out.compareOutput(bsw[idx]);
    out.compareOutput(watdry[idx]);
    out.compareOutput(watopt[idx]);
    //out.compareOutput(rootfr_road_perv[idx]);
    //out.compareOutput(rootr_road_perv[idx]);
    out.compareOutput(watfc[idx]);
    out.compareOutput(displar[idx]);
    out.compareOutput(z0mr[idx]);
  }
  return 0;
}
