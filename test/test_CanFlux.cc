#include "CanopyFluxes.h"
#include "ELMConstants.h"
#include "LandType.h"
#include "read_test_input.hh"
#include "read_input.hh"
#include "ReadPFTConstants.hh"
#include "array.hh"
#include "vegproperties.h"

#include <iostream>
#include <string>

/*

tests CanopyFluxes kernels: 
InitializeFlux_Can()
StabilityIteration_Can()
ComputeFlux_Can()


VegProperties veg

int snl
int frac_veg_nosno
int nrad
int altmax_indx
int altmax_lastyear_indx

double frac_sno
double forc_hgt_u_patch
double thm
double thv
double max_dayl
double dayl
double tc_stress
double elai
double esai
double emv
double emg
double qg
double t_grnd
double forc_t
double forc_pbot
double forc_lwrad
double forc_u
double forc_v
double forc_q
double forc_th
double z0mg
double btran
double displa
double z0mv
double z0hv
double z0qv
double t_veg
double forc_hgt_t_patch
double forc_hgt_q_patch
double fwet
double fdry
double laisun
double laisha
double forc_rho
double snow_depth
double soilbeta
double frac_h2osfc
double t_h2osfc
double sabv
double h2ocan
double htop
double t10
double vcmaxcintsha
double vcmaxcintsun
double forc_pco2
double forc_po2
double qflx_tran_veg
double qflx_evap_veg
double eflx_sh_veg
double qg_snow
double qg_soil
double qg_h2osfc
double dqgdT
double htvp
double eflx_sh_grnd
double eflx_sh_snow
double eflx_sh_soil
double eflx_sh_h2osfc
double qflx_evap_soi
double qflx_ev_snow
double qflx_ev_soil
double qflx_ev_h2osfc
double dlrad
double ulrad
double cgrnds
double cgrndl
double cgrnd
double t_ref2m
double t_ref2m_r
double q_ref2m
double rh_ref2m
double rh_ref2m_r

ArrayD1 rootr
ArrayD1 eff_porosity
ArrayD1 tlai_z
ArrayD1 parsha_z
ArrayD1 parsun_z
ArrayD1 laisha_z
ArrayD1 laisun_z
ArrayD1 t_soisno
ArrayD1 h2osoi_ice
ArrayD1 h2osoi_liq
ArrayD1 dz
ArrayD1 rootfr
ArrayD1 sucsat
ArrayD1 watsat
ArrayD1 bsw



*/

void copy_to_vegstruct(const ELM::Array<double, 1> array, double *vegvar) {
  for (std::size_t i = 0; i < array.extent(0); ++i)
    vegvar[i] = array[i];
}


using ArrayI1 = ELM::Array<int, 1>;
using ArrayD1 = ELM::Array<double, 1>;
using ArrayS1 = ELM::Array<std::string, 1>;
using ArrayD2 = ELM::Array<double, 2>;

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t, typename Scalar_t> void assign(Array_t &arr, Scalar_t val) { ELM::deep_copy(arr, val); }

int main(int argc, char **argv) {

  // data files 
  const std::string data_dir("/Users/80x/Software/elm_kernels/test/data/");
  const std::string input_file = data_dir + "CanopyFluxes_IN.txt";
  const std::string output_file = data_dir + "CanopyFluxes_OUT.txt";
  const std::string pft_file = "clm_params_c180524.nc";


  // hardwired
  ELM::LandType Land;
  ELM::VegProperties Veg;
  Land.ltype = 1;
  Land.ctype = 1;
  Land.vtype = 12;
  int n_grid_cells = 1;
  double dtime = 1800.0;
  int idx = 0;
  int MPI_COMM_WORLD;




  // temporary data to pass between functions
  double wtg = 0.0;         // heat conductance for ground [m/s]
  double wtgq = 0.0;        // latent heat conductance for ground [m/s]
  double wtalq = 0.0;       // normalized latent heat cond. for air and leaf [-]
  double wtlq0 = 0.0;       // normalized latent heat conductance for leaf [-]
  double wtaq0 = 0.0;       // normalized latent heat conductance for air [-]
  double wtl0 = 0.0;        // normalized heat conductance for leaf [-]
  double wta0 = 0.0;        // normalized heat conductance for air [-]
  double wtal = 0.0;        // normalized heat conductance for air and leaf [-]
  double dayl_factor = 0.0; // scalar (0-1) for daylength effect on Vcmax
  double air = 0.0;         // atmos. radiation temporay set
  double bir = 0.0;         // atmos. radiation temporay set
  double cir = 0.0;         // atmos. radiation temporay set
  double el = 0.0;          // vapor pressure on leaf surface [pa]
  double qsatl = 0.0;       // leaf specific humidity [kg/kg]
  double qsatldT = 0.0;     // derivative of "qsatl" on "t_veg"
  double taf = 0.0;         // air temperature within canopy space [K]
  double qaf = 0.0;         // humidity of canopy air [kg/kg]
  double um = 0.0;          // wind speed including the stablity effect [m/s]
  double ur = 0.0;          // wind speed at reference height [m/s]
  double dth = 0.0;         // diff of virtual temp. between ref. height and surface
  double dqh = 0.0;         // diff of humidity between ref. height and surface
  double obu = 0.0;         // Monin-Obukhov length (m)
  double zldis = 0.0;       // reference height "minus" zero displacement height [m]
  double temp1 = 0.0;       // relation for potential temperature profile
  double temp2 = 0.0;       // relation for specific humidity profile
  double temp12m = 0.0;     // relation for potential temperature profile applied at 2-m
  double temp22m = 0.0;     // relation for specific humidity profile applied at 2-m
  double tlbef = 0.0;       // leaf temperature from previous iteration [K]
  double delq = 0.0;        // temporary
  double dt_veg = 0.0;      // change in t_veg, last iteration (Kelvin)

  //if (read_pft_file) {
    // PFT variables
    //const auto maxpfts = ELM::IO::get_maxpfts(MPI_COMM_WORLD, data_dir, pft_file, "z0mr"); // I should write a better method to get the 'pft' dimension var from the params.nc files
    //std::cout << " maxpfts: " << 25 << std::endl;
    auto pftnames = create<ArrayS1>("pftnames", 25);
    auto z0mr = create<ArrayD1>("z0mr", 25);
    auto displar = create<ArrayD1>("displar", 25);
    auto dleaf = create<ArrayD1>("dleaf", 25);
    auto c3psn = create<ArrayD1>("c3psn", 25);
    auto xl = create<ArrayD1>("xl", 25);
    auto roota_par = create<ArrayD1>("roota_par", 25);
    auto rootb_par = create<ArrayD1>("rootb_par", 25);
    auto slatop = create<ArrayD1>("slatop", 25);
    auto leafcn = create<ArrayD1>("leafcn", 25);
    auto flnr = create<ArrayD1>("flnr", 25);
    auto smpso = create<ArrayD1>("smpso", 25);
    auto smpsc = create<ArrayD1>("smpsc", 25);
    auto fnitr = create<ArrayD1>("fnitr", 25);
    auto fnr = create<ArrayD1>("fnr", 25);
    auto act25 = create<ArrayD1>("act25", 25);
    auto kcha = create<ArrayD1>("kcha", 25);
    auto koha = create<ArrayD1>("koha", 25);
    auto cpha = create<ArrayD1>("cpha", 25);
    auto vcmaxha = create<ArrayD1>("vcmaxha", 25);
    auto jmaxha = create<ArrayD1>("jmaxha", 25);
    auto tpuha = create<ArrayD1>("tpuha", 25);
    auto lmrha = create<ArrayD1>("lmrha", 25);
    auto vcmaxhd = create<ArrayD1>("vcmaxhd", 25);
    auto jmaxhd = create<ArrayD1>("jmaxhd", 25);
    auto tpuhd = create<ArrayD1>("tpuhd", 25);
    auto lmrhd = create<ArrayD1>("lmrhd", 25);
    auto lmrse = create<ArrayD1>("lmrse", 25);
    auto qe = create<ArrayD1>("qe", 25);
    auto theta_cj = create<ArrayD1>("theta_cj", 25);
    auto bbbopt = create<ArrayD1>("bbbopt", 25);
    auto mbbopt = create<ArrayD1>("mbbopt", 25);
    auto nstor = create<ArrayD1>("nstor", 25);
    auto br_xr = create<ArrayD1>("br_xr", 25);
    auto tc_stress = create<ArrayD1>("tc_stress", 1); // only one value - keep in container for compatibility with NetCDF reader
    auto rholvis = create<ArrayD1>("rholvis", 25); // numrad
    auto rholnir = create<ArrayD1>("rholnir", 25); // numrad
    auto rhosvis = create<ArrayD1>("rhosvis", 25); // numrad
    auto rhosnir = create<ArrayD1>("rhosnir", 25); // numrad
    auto taulvis = create<ArrayD1>("taulvis", 25); // numrad
    auto taulnir = create<ArrayD1>("taulnir", 25); // numrad
    auto tausvis = create<ArrayD1>("tausvis", 25); // numrad
    auto tausnir = create<ArrayD1>("tausnir", 25); // numrad

    // read pft data from clm_params.nc
    ELM::ReadPFTConstants(data_dir, pft_file, pftnames, z0mr, displar, dleaf, c3psn, xl, roota_par, rootb_par,
                          slatop, leafcn, flnr, smpso, smpsc, fnitr, fnr, act25, kcha, koha, cpha, vcmaxha, jmaxha, tpuha,
                          lmrha, vcmaxhd, jmaxhd, tpuhd, lmrhd, lmrse, qe, theta_cj, bbbopt, mbbopt, nstor, br_xr,
                          tc_stress, rholvis, rholnir, rhosvis, rhosnir, taulvis, taulnir, tausvis, tausnir);
  //}

  // ELM state variables
auto snl = create<ArrayI1>("snl", n_grid_cells);
auto frac_veg_nosno = create<ArrayI1>("frac_veg_nosno", n_grid_cells);
auto nrad = create<ArrayI1>("nrad", n_grid_cells);
auto altmax_indx = create<ArrayI1>("altmax_indx", n_grid_cells);
auto altmax_lastyear_indx = create<ArrayI1>("altmax_lastyear_indx", n_grid_cells);


auto frac_sno = create<ArrayD1>("frac_sno", n_grid_cells);
auto forc_hgt_u_patch = create<ArrayD1>("forc_hgt_u_patch", n_grid_cells);
auto thm = create<ArrayD1>("thm", n_grid_cells);
auto thv = create<ArrayD1>("thv", n_grid_cells);
auto max_dayl = create<ArrayD1>("max_dayl", n_grid_cells);
auto dayl = create<ArrayD1>("dayl", n_grid_cells);
//auto tc_stress = create<ArrayD1>("tc_stress", n_grid_cells);
auto elai = create<ArrayD1>("elai", n_grid_cells);
auto esai = create<ArrayD1>("esai", n_grid_cells);
auto emv = create<ArrayD1>("emv", n_grid_cells);
auto emg = create<ArrayD1>("emg", n_grid_cells);
auto qg = create<ArrayD1>("qg", n_grid_cells);
auto t_grnd = create<ArrayD1>("t_grnd", n_grid_cells);
auto forc_t = create<ArrayD1>("forc_t", n_grid_cells);
auto forc_pbot = create<ArrayD1>("forc_pbot", n_grid_cells);
auto forc_lwrad = create<ArrayD1>("forc_lwrad", n_grid_cells);
auto forc_u = create<ArrayD1>("forc_u", n_grid_cells);
auto forc_v = create<ArrayD1>("forc_v", n_grid_cells);
auto forc_q = create<ArrayD1>("forc_q", n_grid_cells);
auto forc_th = create<ArrayD1>("forc_th", n_grid_cells);
auto z0mg = create<ArrayD1>("z0mg", n_grid_cells);
auto btran = create<ArrayD1>("btran", n_grid_cells);
auto displa = create<ArrayD1>("displa", n_grid_cells);
auto z0mv = create<ArrayD1>("z0mv", n_grid_cells);
auto z0hv = create<ArrayD1>("z0hv", n_grid_cells);
auto z0qv = create<ArrayD1>("z0qv", n_grid_cells);
auto t_veg = create<ArrayD1>("t_veg", n_grid_cells);
auto forc_hgt_t_patch = create<ArrayD1>("forc_hgt_t_patch", n_grid_cells);
auto forc_hgt_q_patch = create<ArrayD1>("forc_hgt_q_patch", n_grid_cells);
auto fwet = create<ArrayD1>("fwet", n_grid_cells);
auto fdry = create<ArrayD1>("fdry", n_grid_cells);
auto laisun = create<ArrayD1>("laisun", n_grid_cells);
auto laisha = create<ArrayD1>("laisha", n_grid_cells);
auto forc_rho = create<ArrayD1>("forc_rho", n_grid_cells);
auto snow_depth = create<ArrayD1>("snow_depth", n_grid_cells);
auto soilbeta = create<ArrayD1>("soilbeta", n_grid_cells);
auto frac_h2osfc = create<ArrayD1>("frac_h2osfc", n_grid_cells);
auto t_h2osfc = create<ArrayD1>("t_h2osfc", n_grid_cells);
auto sabv = create<ArrayD1>("sabv", n_grid_cells);
auto h2ocan = create<ArrayD1>("h2ocan", n_grid_cells);
auto htop = create<ArrayD1>("htop", n_grid_cells);
auto t10 = create<ArrayD1>("t10", n_grid_cells);
auto vcmaxcintsha = create<ArrayD1>("vcmaxcintsha", n_grid_cells);
auto vcmaxcintsun = create<ArrayD1>("vcmaxcintsun", n_grid_cells);
auto forc_pco2 = create<ArrayD1>("forc_pco2", n_grid_cells);
auto forc_po2 = create<ArrayD1>("forc_po2", n_grid_cells);
auto qflx_tran_veg = create<ArrayD1>("qflx_tran_veg", n_grid_cells);
auto qflx_evap_veg = create<ArrayD1>("qflx_evap_veg", n_grid_cells);
auto eflx_sh_veg = create<ArrayD1>("eflx_sh_veg", n_grid_cells);
auto qg_snow = create<ArrayD1>("qg_snow", n_grid_cells);
auto qg_soil = create<ArrayD1>("qg_soil", n_grid_cells);
auto qg_h2osfc = create<ArrayD1>("qg_h2osfc", n_grid_cells);
auto dqgdT = create<ArrayD1>("dqgdT", n_grid_cells);
auto htvp = create<ArrayD1>("htvp", n_grid_cells);
auto eflx_sh_grnd = create<ArrayD1>("eflx_sh_grnd", n_grid_cells);
auto eflx_sh_snow = create<ArrayD1>("eflx_sh_snow", n_grid_cells);
auto eflx_sh_soil = create<ArrayD1>("eflx_sh_soil", n_grid_cells);
auto eflx_sh_h2osfc = create<ArrayD1>("eflx_sh_h2osfc", n_grid_cells);
auto qflx_evap_soi = create<ArrayD1>("qflx_evap_soi", n_grid_cells);
auto qflx_ev_snow = create<ArrayD1>("qflx_ev_snow", n_grid_cells);
auto qflx_ev_soil = create<ArrayD1>("qflx_ev_soil", n_grid_cells);
auto qflx_ev_h2osfc = create<ArrayD1>("qflx_ev_h2osfc", n_grid_cells);
auto dlrad = create<ArrayD1>("dlrad", n_grid_cells);
auto ulrad = create<ArrayD1>("ulrad", n_grid_cells);
auto cgrnds = create<ArrayD1>("cgrnds", n_grid_cells);
auto cgrndl = create<ArrayD1>("cgrndl", n_grid_cells);
auto cgrnd = create<ArrayD1>("cgrnd", n_grid_cells);
auto t_ref2m = create<ArrayD1>("t_ref2m", n_grid_cells);
auto t_ref2m_r = create<ArrayD1>("t_ref2m_r", n_grid_cells);
auto q_ref2m = create<ArrayD1>("q_ref2m", n_grid_cells);
auto rh_ref2m = create<ArrayD1>("rh_ref2m", n_grid_cells);
auto rh_ref2m_r = create<ArrayD1>("rh_ref2m_r", n_grid_cells);




auto rootr = create<ArrayD2>("rootr", n_grid_cells, ELM::nlevgrnd);
auto eff_porosity = create<ArrayD2>("eff_porosity", n_grid_cells, ELM::nlevgrnd);
auto tlai_z = create<ArrayD2>("tlai_z", n_grid_cells, ELM::nlevcan);
auto parsha_z = create<ArrayD2>("parsha_z", n_grid_cells, ELM::nlevcan);
auto parsun_z = create<ArrayD2>("parsun_z", n_grid_cells, ELM::nlevcan);
auto laisha_z = create<ArrayD2>("laisha_z", n_grid_cells, ELM::nlevcan);
auto laisun_z = create<ArrayD2>("laisun_z", n_grid_cells, ELM::nlevcan);
auto t_soisno = create<ArrayD2>("t_soisno", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto h2osoi_ice = create<ArrayD2>("h2osoi_ice", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto h2osoi_liq = create<ArrayD2>("h2osoi_liq", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto dz = create<ArrayD2>("dz", n_grid_cells, ELM::nlevsno + ELM::nlevgrnd);
auto rootfr = create<ArrayD2>("rootfr", n_grid_cells, ELM::nlevgrnd);
auto sucsat = create<ArrayD2>("sucsat", n_grid_cells, ELM::nlevgrnd);
auto watsat = create<ArrayD2>("watsat", n_grid_cells, ELM::nlevgrnd);
auto bsw = create<ArrayD2>("bsw", n_grid_cells, ELM::nlevgrnd);


copy_to_vegstruct(fnr,Veg.fnr);
copy_to_vegstruct(act25,Veg.act25);
copy_to_vegstruct(kcha,Veg.kcha);
copy_to_vegstruct(koha,Veg.koha);
copy_to_vegstruct(cpha,Veg.cpha);
copy_to_vegstruct(vcmaxha,Veg.vcmaxha);
copy_to_vegstruct(jmaxha,Veg.jmaxha);
copy_to_vegstruct(tpuha,Veg.tpuha);
copy_to_vegstruct(lmrha,Veg.lmrha);
copy_to_vegstruct(vcmaxhd,Veg.vcmaxhd);
copy_to_vegstruct(jmaxhd,Veg.jmaxhd);
copy_to_vegstruct(tpuhd,Veg.tpuhd);
copy_to_vegstruct(lmrhd,Veg.lmrhd);
copy_to_vegstruct(lmrse,Veg.lmrse);
copy_to_vegstruct(qe,Veg.qe);
copy_to_vegstruct(theta_cj,Veg.theta_cj);
copy_to_vegstruct(bbbopt,Veg.bbbopt);
copy_to_vegstruct(mbbopt,Veg.mbbopt);
copy_to_vegstruct(c3psn,Veg.c3psn);
copy_to_vegstruct(slatop,Veg.slatop);
copy_to_vegstruct(leafcn,Veg.leafcn);
copy_to_vegstruct(flnr,Veg.flnr);
copy_to_vegstruct(fnitr,Veg.fnitr);


// input and output utility class objects
ELM::IO::ELMtestinput in(input_file);
ELM::IO::ELMtestinput out(output_file);
for (std::size_t t = 1; t < 49; ++t) {
  // get input and output state for time t
  in.getState(t);
  out.getState(t);


  // parse input state and assign to variables
in.parseState(snl);
in.parseState(frac_veg_nosno);
in.parseState(nrad);
in.parseState(altmax_indx);
in.parseState(altmax_lastyear_indx);
in.parseState(frac_sno);
in.parseState(forc_hgt_u_patch);
in.parseState(thm);
in.parseState(thv);
in.parseState(max_dayl);
in.parseState(dayl);
//in.parseState(tc_stress);
in.parseState(elai);
in.parseState(esai);
in.parseState(emv);
in.parseState(emg);
in.parseState(qg);
in.parseState(t_grnd);
in.parseState(forc_t);
in.parseState(forc_pbot);
in.parseState(forc_lwrad);
in.parseState(forc_u);
in.parseState(forc_v);
in.parseState(forc_q);
in.parseState(forc_th);
in.parseState(z0mg);
in.parseState(btran);
in.parseState(displa);
in.parseState(z0mv);
in.parseState(z0hv);
in.parseState(z0qv);
in.parseState(t_veg);
in.parseState(forc_hgt_t_patch);
in.parseState(forc_hgt_q_patch);
in.parseState(fwet);
in.parseState(fdry);
in.parseState(laisun);
in.parseState(laisha);
in.parseState(forc_rho);
in.parseState(snow_depth);
in.parseState(soilbeta);
in.parseState(frac_h2osfc);
in.parseState(t_h2osfc);
in.parseState(sabv);
in.parseState(h2ocan);
in.parseState(htop);
in.parseState(t10);
in.parseState(vcmaxcintsha);
in.parseState(vcmaxcintsun);
in.parseState(forc_pco2);
in.parseState(forc_po2);
in.parseState(qflx_tran_veg);
in.parseState(qflx_evap_veg);
in.parseState(eflx_sh_veg);
in.parseState(qg_snow);
in.parseState(qg_soil);
in.parseState(qg_h2osfc);
in.parseState(dqgdT);
in.parseState(htvp);
in.parseState(eflx_sh_grnd);
in.parseState(eflx_sh_snow);
in.parseState(eflx_sh_soil);
in.parseState(eflx_sh_h2osfc);
in.parseState(qflx_evap_soi);
in.parseState(qflx_ev_snow);
in.parseState(qflx_ev_soil);
in.parseState(qflx_ev_h2osfc);
in.parseState(dlrad);
in.parseState(ulrad);
in.parseState(cgrnds);
in.parseState(cgrndl);
in.parseState(cgrnd);
in.parseState(t_ref2m);
in.parseState(t_ref2m_r);
in.parseState(q_ref2m);
in.parseState(rh_ref2m);
in.parseState(rh_ref2m_r);
in.parseState(rootr[idx]);
in.parseState(eff_porosity[idx]);
in.parseState(tlai_z[idx]);
in.parseState(parsha_z[idx]);
in.parseState(parsun_z[idx]);
in.parseState(laisha_z[idx]);
in.parseState(laisun_z[idx]);
in.parseState(t_soisno[idx]);
in.parseState(h2osoi_ice[idx]);
in.parseState(h2osoi_liq[idx]);
in.parseState(dz[idx]);
in.parseState(rootfr[idx]);
in.parseState(sucsat[idx]);
in.parseState(watsat[idx]);
in.parseState(bsw[idx]);


    // call CanopyFluxes kernels
    ELM::InitializeFlux_Can(
        Land, snl[idx], frac_veg_nosno[idx], frac_sno[idx], forc_hgt_u_patch[idx], thm[idx], thv[idx], max_dayl[idx], dayl[idx],
        altmax_indx[idx], altmax_lastyear_indx[idx], t_soisno[idx],
        h2osoi_ice[idx], h2osoi_liq[idx],
        dz[idx], rootfr[idx], tc_stress[idx],
        sucsat[idx], watsat[idx],
        bsw[idx], smpso,
        smpsc, elai[idx], esai[idx], emv[idx], emg[idx], qg[idx], t_grnd[idx], forc_t[idx],
        forc_pbot[idx], forc_lwrad[idx], forc_u[idx], forc_v[idx], forc_q[idx], forc_th[idx], z0mg[idx], btran[idx], displa[idx], z0mv[idx],
        z0hv[idx], z0qv[idx], rootr[idx], eff_porosity[idx],
        dayl_factor, air, bir, cir, el, qsatl, qsatldT, taf, qaf, um, ur, obu, zldis, delq, t_veg[idx]);

    ELM::StabilityIteration_Can(
        Land, dtime, snl[idx], frac_veg_nosno[idx], frac_sno[idx], forc_hgt_u_patch[idx], forc_hgt_t_patch[idx],
        forc_hgt_q_patch[idx], dleaf, fwet[idx], fdry[idx], laisun[idx], laisha[idx],
        forc_rho[idx], snow_depth[idx], soilbeta[idx], frac_h2osfc[idx], t_h2osfc[idx], sabv[idx], h2ocan[idx], htop[idx],
        t_soisno[idx], air, bir, cir, ur, zldis, displa[idx], elai[idx], esai[idx], t_grnd[idx],
        forc_pbot[idx], forc_q[idx], forc_th[idx], z0mg[idx], z0mv[idx], z0hv[idx], z0qv[idx], thm[idx], thv[idx], qg[idx], Veg, nrad[idx],
        t10[idx], tlai_z[idx], vcmaxcintsha[idx], vcmaxcintsun[idx],
        parsha_z[idx], parsun_z[idx],
        laisha_z[idx], laisun_z[idx], forc_pco2[idx], forc_po2[idx],
        dayl_factor, btran[idx], qflx_tran_veg[idx], qflx_evap_veg[idx], eflx_sh_veg[idx], wtg, wtl0, wta0, wtal, el, qsatl,
        qsatldT, taf, qaf, um, dth, dqh, obu, temp1, temp2, temp12m, temp22m, tlbef, delq, dt_veg, t_veg[idx], wtgq,
        wtalq, wtlq0, wtaq0);

    ELM::ComputeFlux_Can(
        Land, dtime, snl[idx], frac_veg_nosno[idx], frac_sno[idx], t_soisno[idx], frac_h2osfc[idx],
        t_h2osfc[idx], sabv[idx], qg_snow[idx], qg_soil[idx], qg_h2osfc[idx], dqgdT[idx], htvp[idx], wtg, wtl0, wta0, wtal, air, bir,
        cir, qsatl, qsatldT, dth, dqh, temp1, temp2, temp12m, temp22m, tlbef, delq, dt_veg, t_veg[idx], t_grnd[idx],
        forc_pbot[idx], qflx_tran_veg[idx], qflx_evap_veg[idx], eflx_sh_veg[idx], forc_q[idx], forc_rho[idx], thm[idx], emv[idx],
        emg[idx], forc_lwrad[idx], wtgq, wtalq, wtlq0, wtaq0, h2ocan[idx], eflx_sh_grnd[idx], eflx_sh_snow[idx], eflx_sh_soil[idx],
        eflx_sh_h2osfc[idx], qflx_evap_soi[idx], qflx_ev_snow[idx], qflx_ev_soil[idx], qflx_ev_h2osfc[idx], dlrad[idx], ulrad[idx],
        cgrnds[idx], cgrndl[idx], cgrnd[idx], t_ref2m[idx], t_ref2m_r[idx], q_ref2m[idx], rh_ref2m[idx], rh_ref2m_r[idx]);


// compare kernel output to ELM output state
out.compareOutput(snl);
out.compareOutput(frac_veg_nosno);
out.compareOutput(nrad);
out.compareOutput(altmax_indx);
out.compareOutput(altmax_lastyear_indx);
out.compareOutput(frac_sno);
out.compareOutput(forc_hgt_u_patch);
out.compareOutput(thm);
out.compareOutput(thv);
out.compareOutput(max_dayl);
out.compareOutput(dayl);
//out.compareOutput(tc_stress);
out.compareOutput(elai);
out.compareOutput(esai);
out.compareOutput(emv);
out.compareOutput(emg);
out.compareOutput(qg);
out.compareOutput(t_grnd);
out.compareOutput(forc_t);
out.compareOutput(forc_pbot);
out.compareOutput(forc_lwrad);
out.compareOutput(forc_u);
out.compareOutput(forc_v);
out.compareOutput(forc_q);
out.compareOutput(forc_th);
out.compareOutput(z0mg);
out.compareOutput(btran);
out.compareOutput(displa);
out.compareOutput(z0mv);
out.compareOutput(z0hv);
out.compareOutput(z0qv);
out.compareOutput(t_veg);
out.compareOutput(forc_hgt_t_patch);
out.compareOutput(forc_hgt_q_patch);
out.compareOutput(fwet);
out.compareOutput(fdry);
out.compareOutput(laisun);
out.compareOutput(laisha);
out.compareOutput(forc_rho);
out.compareOutput(snow_depth);
out.compareOutput(soilbeta);
out.compareOutput(frac_h2osfc);
out.compareOutput(t_h2osfc);
out.compareOutput(sabv);
out.compareOutput(h2ocan);
out.compareOutput(htop);
out.compareOutput(t10);
out.compareOutput(vcmaxcintsha);
out.compareOutput(vcmaxcintsun);
out.compareOutput(forc_pco2);
out.compareOutput(forc_po2);
out.compareOutput(qflx_tran_veg);
out.compareOutput(qflx_evap_veg);
out.compareOutput(eflx_sh_veg);
out.compareOutput(qg_snow);
out.compareOutput(qg_soil);
out.compareOutput(qg_h2osfc);
out.compareOutput(dqgdT);
out.compareOutput(htvp);
out.compareOutput(eflx_sh_grnd);
out.compareOutput(eflx_sh_snow);
out.compareOutput(eflx_sh_soil);
out.compareOutput(eflx_sh_h2osfc);
out.compareOutput(qflx_evap_soi);
out.compareOutput(qflx_ev_snow);
out.compareOutput(qflx_ev_soil);
out.compareOutput(qflx_ev_h2osfc);
out.compareOutput(dlrad);
out.compareOutput(ulrad);
out.compareOutput(cgrnds);
out.compareOutput(cgrndl);
out.compareOutput(cgrnd);
out.compareOutput(t_ref2m);
out.compareOutput(t_ref2m_r);
out.compareOutput(q_ref2m);
out.compareOutput(rh_ref2m);
out.compareOutput(rh_ref2m_r);
out.compareOutput(rootr[idx]);
out.compareOutput(eff_porosity[idx]);
out.compareOutput(tlai_z[idx]);
out.compareOutput(parsha_z[idx]);
out.compareOutput(parsun_z[idx]);
out.compareOutput(laisha_z[idx]);
out.compareOutput(laisun_z[idx]);
out.compareOutput(t_soisno[idx]);
out.compareOutput(h2osoi_ice[idx]);
out.compareOutput(h2osoi_liq[idx]);
out.compareOutput(dz[idx]);
out.compareOutput(rootfr[idx]);
out.compareOutput(sucsat[idx]);
out.compareOutput(watsat[idx]);
out.compareOutput(bsw[idx]);

}

}


