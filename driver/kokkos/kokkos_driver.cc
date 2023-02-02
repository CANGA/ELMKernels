
#include <iostream>
#include <string>
#include <unordered_map>

// utilities
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"
#include "helper_functions.hh"

// constants
#include "elm_constants.h"

// initialization routines
#include "initialize_elm_kokkos.hh"
#include "init_timestep_kokkos.hh"

// parallel physics
#include "aerosol_kokkos.hh"
#include "albedo_kokkos.hh"
#include "canopy_hydrology_kokkos.hh"
#include "surface_radiation_kokkos.hh"
#include "canopy_temperature_kokkos.hh"
#include "bareground_fluxes_kokkos.hh"
#include "canopy_fluxes_kokkos.hh"
#include "soil_temperature_kokkos.hh"
#include "snow_hydrology_kokkos.hh"
#include "surface_fluxes_kokkos.hh"
#include "conserved_quantity_kokkos.hh"

// conditional compilation options
#include "compile_options.hh"
#include "data_types.hh"

#include "elm_kokkos_interface.hh"

using ELM::Utils::assign;
using AtmForcType = ELM::AtmForcType;

int main(int argc, char **argv) {

{ // enclosing scope

  using namespace ELM::ELMdims;

  Kokkos::initialize(argc, argv);

  { // inner scope

    std::string input_dir = INPUT_DATA_DIR;

    std::string fname_surfdata(
     input_dir+"E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc");
    std::string fname_snicar(
      input_dir+"pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
    std::string fname_forc(
     input_dir+"pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc");
    std::string fname_param(
      input_dir+"E3SM/components/elm/test_submodules/inputdata/lnd/clm2/paramdata/clm_params_c180524.nc");
    std::string fname_aerosol(
      input_dir+"pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc");
    std::string fname_snowage(
      input_dir+"pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc");

    const int n_procs = 1;
    const int ncols = 1;
    int idx = 0; // hardwire for ncols = 1
    const int ntimes = 100;
    const int myrank = 0;
    const double dtime = 1800.0;
    const double dtime_d = 1800.0 / 86400.0;
    const auto start = ELM::Utils::Date(2014, 1, 1);

    auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
    // domain and processor topology info
    auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { 1, 1 },
          { 0, 0 });

    // number of forcing timesteps to store in device memory
    int atm_nsteps = 101;
    // starting time of forcing file
    const auto fstart = ELM::Utils::Date(1985, 1, 1);
    // ELM State
    auto S = std::make_shared<ELMStateType>(ncols, dd, fname_forc, fstart, atm_nsteps);

    // fix this soon
    //S.get()->Land.ltype = 1;
    S.get()->Land.ltype = 1;
    S.get()->Land.ctype = 1;
    S.get()->Land.vtype = 12;
    S.get()->Land.lakpoi = false;
    S.get()->Land.urbpoi = false;

    assign(S.get()->vtype, 12);

    // hardwired params
    S.get()->lat = 71.323;
    S.get()->lon = 203.3886;
    S.get()->lat_r = S.get()->lat * ELM::ELMconst::ELM_PI() / 180.0;
    S.get()->lon_r = S.get()->lon * ELM::ELMconst::ELM_PI() / 180.0;
    const double dewmx = 0.1;
    const double irrig_rate = 0.0;
    const int n_irrig_steps_left = 0;
    const int oldfflag = 1;
    //auto veg_active = create<ViewB1>("veg_active", ncols); // need value
    assign(S.get()->veg_active, true);                               // hardwired
    //auto do_capsnow = create<ViewB1>("do_capsnow", ncols); // need value
    assign(S.get()->do_capsnow, false);                               // hardwired
    assign(S.get()->topo_slope, 0.070044865858546);
    assign(S.get()->topo_std, 3.96141847422387);
    assign(S.get()->snl, 0);
    assign(S.get()->snow_depth, 0.0);
    assign(S.get()->frac_sno, 0.0);
    assign(S.get()->int_snow, 0.0);
    assign(S.get()->h2osoi_liq, 0.0);
    assign(S.get()->h2osoi_ice, 0.0);
    assign(S.get()->t_h2osfc, 274.0);

    assign(S.get()->eflx_sh_grnd, 0.0);
    assign(S.get()->eflx_sh_snow, 0.0);
    assign(S.get()->eflx_sh_soil, 0.0);
    assign(S.get()->eflx_sh_h2osfc, 0.0);
    assign(S.get()->qflx_evap_soi, 0.0);
    assign(S.get()->qflx_ev_snow, 0.0);
    assign(S.get()->qflx_ev_soil, 0.0);
    assign(S.get()->qflx_ev_h2osfc, 0.0);

    assign(S.get()->altmax_indx, 5);
    assign(S.get()->altmax_lastyear_indx, 0);
    assign(S.get()->t10, 276.0);
    assign(S.get()->t_veg, 283.0);

    assign(S.get()->xmf, 0.0);
    assign(S.get()->xmf_h2osfc, 0.0);
    assign(S.get()->eflx_h2osfc_snow, 0.0);

    // hardwired grid info
    // this comes from ELM, but is wrong?
    // doesn't matter for now, but definitely wrong
    {
      double dz_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.017512817916255204,
      0.02757896925967625, 0.0454700332424132, 0.07496741098620856,
      0.12360036510228053, 0.20378255101043175, 0.33598062644843263,
      0.5539384053686849, 0.9132900315890611, 1.5057607013992766,
      2.482579696981332, 4.0930819526214, 6.7483512780057175,
      11.12615029420442, 13.851152141963599 };
      auto h_dz = Kokkos::create_mirror_view(S.get()->dz);
      for (int n = 0; n < ncols; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_dz(n, i) = dz_hardwire[i];
        }
      }
      Kokkos::deep_copy(S.get()->dz, h_dz);

      double zsoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.007100635417193535,
      0.02792500041531687, 0.06225857393654604, 0.11886506690014327,
      0.21219339590896316, 0.3660657971047043, 0.6197584979298266,
      1.0380270500015696, 1.7276353086671965, 2.8646071131796917,
      4.73915671146575, 7.829766507142356, 12.92532061670855,
      21.32646906315379, 35.17762120511739 };
      auto h_zsoi = Kokkos::create_mirror_view(S.get()->zsoi);
      for (int n = 0; n < ncols; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_zsoi(n, i) = zsoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(S.get()->zsoi, h_zsoi);

      double zisoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.017512817916255204, 0.04509178717593146, 0.09056182041834465, 
      0.16552923140455322, 0.28912959650683373, 0.4929121475172655,
      0.8288927739656982, 1.382831179334383, 2.2961212109234443,
      3.8018819123227208, 6.284461609304053, 10.377543561925453,
      17.12589483993117, 28.252045134135592, 42.10319727609919 };
      auto h_zisoi = Kokkos::create_mirror_view(S.get()->zisoi);
      for (int n = 0; n < ncols; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd() + 1; ++i) {
          h_zisoi(n, i) = zisoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(S.get()->zisoi, h_zisoi);
    }

    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // initialize data containers and read time-invariant data from files
    // these calls only need to occur once @ beginning of simulation
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

    ELM::initialize_kokkos_elm(*S.get(), fname_surfdata,
      fname_param, fname_snicar, fname_snowage, fname_aerosol);

    // hardwired soil/snow water state
    {
      double h2osoi_ice_hardwire[] = {
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 51.095355179469955, 131.99213225849098,
        17.829256395227745, 95.72899575304584, 155.31526899797177,
        0.01, 0.01, 0.01,
        0.01, 0.01 };

      double h2osoi_liq_hardwire[] = {
        0.0, 0.0, 0.0,
        0.0, 0.0, 7.045411435071487,
        14.353496179256807, 36.308518784697064, 62.46145027256513,
        97.14000248023912, 97.47148319510016, 78.52160092062527,
        65.63904088905001, 41.25305599181871, 70.8566046019581,
        0.01, 0.01, 0.01,
        0.01, 0.01 };

      auto h_soi_ice = Kokkos::create_mirror_view(S.get()->h2osoi_ice);
      auto h_soi_liq = Kokkos::create_mirror_view(S.get()->h2osoi_liq);
      for (int n = 0; n < ncols; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_soi_ice(n, i) = h2osoi_ice_hardwire[i];
          h_soi_liq(n, i) = h2osoi_liq_hardwire[i];
        }
      }
      Kokkos::deep_copy(S.get()->h2osoi_ice, h_soi_ice);
      Kokkos::deep_copy(S.get()->h2osoi_liq, h_soi_liq);


      double h2osoi_vol_hardwire[] = {
        0.4016484663460637, 0.5196481455614503, 0.7967166638201649,
        0.8331813710901114, 0.7859200286330449, 0.7517405589446893,
        0.6621235242027332, 0.1535948180493002, 0.15947477948341815,
        0.15954052527228618, 8.420726808634413e-06, 5.107428986500891e-06,
        3.0978122726178113e-06, 1.8789181213767733e-06, 1.5092697845407248e-06 };
      auto h_soi_vol = Kokkos::create_mirror_view(S.get()->h2osoi_vol);
      for (int n = 0; n < ncols; ++n) {
        for (int i = 0; i < nlevgrnd(); ++i) {
          h_soi_vol(n, i) = h2osoi_vol_hardwire[i];
        }
      }
      Kokkos::deep_copy(S.get()->h2osoi_vol, h_soi_vol);

    }


    // hardwired soil/snow temp info
    {
     double tsoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 278.3081064745931,
      276.1568781897738, 275.55803480737063, 275.2677090940866,
      274.7286996980052, 273.15, 272.4187794248787,
      270.65049816473027, 267.8224112387398, 265.7450135695632,
      264.49481140089864, 264.14163363048056, 264.3351872934207,
      264.1163763444719, 263.88852987294865 };
      auto h_tsoi = Kokkos::create_mirror_view(S.get()->t_soisno);
      for (int n = 0; n < ncols; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_tsoi(n, i) = tsoi_hardwire[i];
        }
      }
      auto h_tgrnd = Kokkos::create_mirror_view(S.get()->t_grnd);
      h_tgrnd(idx) = h_tsoi(idx, nlevsno() - S.get()->snl(idx));
      Kokkos::deep_copy(S.get()->t_soisno, h_tsoi);
      Kokkos::deep_copy(S.get()->t_grnd, h_tgrnd);
    }




    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /*                          TIME LOOP                                                                  */
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    ELM::Utils::Date current(start);

    for (int t = 0; t < ntimes; ++t) {

      ELM::Utils::Date time_plus_half_dt(current);
      time_plus_half_dt.increment_seconds(static_cast<size_t>((dtime + 1.0)/2)); // round to nearest second

      // get coszen, day length,
      // phenology data, atmospheric forcing,
      // aerosol forcing and snowpack state,
      // and a handful of variables that need to be reset every timestep
      ELM::kokkos_init_timestep(S.get(), dtime, current, fname_surfdata);


      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // Main physics calls
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

      {
        // canhydro::fraction_wet
        ELM::kokkos_frac_wet(*S.get());

        // call surface albedo and SNICAR kernels
        ELM::kokkos_albedo_snicar(*S.get());

        // call canopy_hydrology kernels
        ELM::kokkos_canopy_hydrology(*S.get(), dtime);

        // call surface_radiation kernels
        ELM::kokkos_surface_radiation(*S.get());

        // call canopy_temperature kernels
        ELM::kokkos_canopy_temperature(*S.get());

        // call bareground_fluxes kernels
        ELM::kokkos_bareground_fluxes(*S.get());

        // call canopy_fluxes kernels
        ELM::kokkos_canopy_fluxes(*S.get(), dtime);

        // call soil_temperature kernels
        ELM::kokkos_soil_temperature(*S.get(), dtime);

        // call snow_hydrology kernels
        ELM::kokkos_snow_hydrology(*S.get(), dtime, time_plus_half_dt);

        // call surface_fluxes kernels
        ELM::kokkos_surface_fluxes(*S.get(), dtime);

        ELM::kokkos_evaluate_conservation(*S.get(), dtime);
      }

      // print diagnostics

  //    std::cout << "soil temp fluxes:  "
  //        << hs_soil(0) << "  "
  //        << hs_h2osfc(0) << "  "
  //        << hs_top(0) << "  "
  //        << hs_top_snow(0) << "  "
  //        << dhsdT(0) << "  "
  //        << soil_e_balance(0) << "  "
  //        << sabg_chk(0) << std::endl;

      std::cout << "lwrad_out:  "
          << S.get()->eflx_lwrad_out(0) << "  "
          << S.get()->eflx_lwrad_net(0) << std::endl;

      for (int i = 0; i < nlevsno() + nlevgrnd(); ++i)
        std::cout << "column vars:  " << i <<
         "  t_soisno:  " << S.get()->t_soisno(0, i) <<
         "  h2osoi_ice:  " << S.get()->h2osoi_ice(0, i) <<
         "  h2osoi_liq:  " << S.get()->h2osoi_liq(0, i) <<
         "  dz:  " << S.get()->dz(0, i) <<
         "  zsoi:  " << S.get()->zsoi(0, i) <<
         "  zisoi:  " << S.get()->zisoi(0, i) << std::endl;
         std::cout << "last   zisoi:  " << S.get()->zisoi(0, nlevsno()+nlevgrnd()) << std::endl;




         for (int i = 0; i < nlevgrnd(); ++i)
         std::cout << i <<"  qflx_rootsoi:  " << S.get()->qflx_rootsoi(0, i) << "\n";


      for (int i = 0; i < ncols; ++i) {
        std::cout << "h2osno: " << S.get()->h2osno(i) << std::endl;
        std::cout << "t_grnd: " << S.get()->t_grnd(i) << std::endl;
        std::cout << "t_h2osfc: " << S.get()->t_h2osfc(i) << std::endl;
        std::cout << "t10: " << S.get()->t10(i) << std::endl;
        std::cout << "t_veg: " << S.get()->t_veg(i) << std::endl;
        std::cout << "snow_depth: " << S.get()->snow_depth(i) << std::endl;
        std::cout << "frac_sno: " << S.get()->frac_sno(i) << std::endl;
        std::cout << "frac_sno_eff: " << S.get()->frac_sno_eff(i) << std::endl;
        std::cout << "qflx_tran_veg: " << S.get()->qflx_tran_veg(i) << std::endl;
        std::cout << "frac_veg_nosno_alb: " << S.get()->frac_veg_nosno_alb(i) << std::endl;
      }

      current.increment_seconds(dtime);

    } // time loop


    // test interface
    ELM::ELMInterface test(ncols);
    test.setup();
    auto state_ptr = test.getPrimaryVars();
    std::cout << "test vars  !!!!!!!!!!!!! " << std::endl;
    std::cout << "snl " << state_ptr.get()->snl(0) << std::endl;
    std::cout << "snow_depth " << state_ptr.get()->snow_depth(0) << std::endl;
    std::cout << "frac_sno " << state_ptr.get()->frac_sno(0) << std::endl;
    std::cout << "int_snow " << state_ptr.get()->int_snow(0) << std::endl;
    std::cout << "h2ocan " << state_ptr.get()->h2ocan(0) << std::endl;
    std::cout << "h2osno " << state_ptr.get()->h2osno(0) << std::endl;
    std::cout << "h2osfc " << state_ptr.get()->h2osfc(0) << std::endl;
    std::cout << "t_grnd " << state_ptr.get()->t_grnd(0) << std::endl;
    std::cout << "t_h2osfc " << state_ptr.get()->t_h2osfc(0) << std::endl;
    std::cout << "t_h2osfc_bef " << state_ptr.get()->t_h2osfc_bef(0) << std::endl;
  } // inner scope

  Kokkos::finalize();
} // enclosing scope
return 0;
}
