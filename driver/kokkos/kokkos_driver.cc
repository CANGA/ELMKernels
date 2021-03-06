
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

// input data readers and structs
#include "pft_data.h"
#include "snicar_data.h"
#include "aerosol_physics.h"
#include "aerosol_data.h"
#include "phenology_data.h"

// atm data helper functions
// contains definition of atm data struct 
// that should be moved to atm_data.h 
#include "atm_data_helpers.hh"

// initialization routines
#include "initialize_elm_kokkos.hh"
#include "init_timestep_kokkos.hh"

// serial physics
#include "day_length.h"
#include "incident_shortwave.h"

// parallel physics
#include "albedo_kokkos.hh"
#include "canopy_hydrology_kokkos.hh"
#include "surface_radiation_kokkos.hh"
#include "canopy_temperature_kokkos.hh"
#include "bareground_fluxes_kokkos.hh"
#include "canopy_fluxes_kokkos.hh"
#include "soil_temperature_kokkos.hh"
#include "snow_hydrology_kokkos.hh"
#include "surface_fluxes_kokkos.hh"

// conditional compilation options
#include "invoke_kernel.hh"
#include "compile_options.hh"
#include "data_types.hh"

// elm state struct
#include "elm_state.h"

using ELM::Utils::create;
using ELM::Utils::assign;
using AtmForcType = ELM::AtmForcType;

// these phenology data functions will stay here until I decide what to do about host-device transfer
std::unordered_map<std::string, h_ViewD2> get_phen_host_views(const ELM::PhenologyDataManager<ViewD2>& phen_data)
{
  std::unordered_map<std::string, h_ViewD2> phen_host_views;
  phen_host_views["MONTHLY_LAI"] = Kokkos::create_mirror_view(phen_data.mlai);
  phen_host_views["MONTHLY_SAI"] = Kokkos::create_mirror_view(phen_data.msai);
  phen_host_views["MONTHLY_HEIGHT_TOP"] = Kokkos::create_mirror_view(phen_data.mhtop);
  phen_host_views["MONTHLY_HEIGHT_BOT"] = Kokkos::create_mirror_view(phen_data.mhbot);
  return phen_host_views;
}

void update_phenology(ELM::PhenologyDataManager<ViewD2>& phen_data,
                      std::unordered_map<std::string, h_ViewD2>& host_phen_views,
                      const std::shared_ptr<ELMStateType>& S,
                      const h_ViewI1 vtype,
                      const ELM::Utils::Date& current,
                      const std::string& fname_surfdata)
{
  // copy device data to host
  // copying entire views is likely inefficient, but it's currently necessary
  // could be eliminated by shifting older month indices in parallel kernel
  // and reading new data into a mirror of a subview (or a subview of a mirror?)
  // then we would only need one copy from host view into the device view
  // instead of the two we currently have
  // will fix later - too infrequently run (once per month) to cause concern
  if (phen_data.need_data()) {
    Kokkos::deep_copy(host_phen_views["MONTHLY_LAI"], phen_data.mlai);
    Kokkos::deep_copy(host_phen_views["MONTHLY_SAI"], phen_data.msai);
    Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_TOP"], phen_data.mhtop);
    Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_BOT"], phen_data.mhbot);
  }
  // reads three months of data on first call
  // after first call, read new data if phen_data.need_new_data_ == true
  auto phen_updated = phen_data.read_data(host_phen_views, fname_surfdata, current, vtype); // if needed
  // copy host views to device
  // could be made more efficient, see above
  if (phen_updated) {
    Kokkos::deep_copy(phen_data.mlai, host_phen_views["MONTHLY_LAI"]);
    Kokkos::deep_copy(phen_data.msai, host_phen_views["MONTHLY_SAI"]);
    Kokkos::deep_copy(phen_data.mhtop, host_phen_views["MONTHLY_HEIGHT_TOP"]);
    Kokkos::deep_copy(phen_data.mhbot, host_phen_views["MONTHLY_HEIGHT_BOT"]);
  }
  // run parallel kernel to process phenology data
  phen_data.get_data(current, S->snow_depth,
                     S->frac_sno, S->vtype, S->elai, S->esai,
                     S->htop, S->hbot, S->tlai, S->tsai,
                     S->frac_veg_nosno_alb);
}





int main(int argc, char **argv) {

{ // enclosing scope

  using namespace ELM::ELMdims;

  Kokkos::initialize(argc, argv);

  { // inner scope

    std::string fname_surfdata(
    "/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc");
    std::string fname_snicar(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
    std::string fname_forc(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc");
    std::string fname_param(
      "/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata/lnd/clm2/paramdata/clm_params_c180524.nc");
    std::string fname_aerosol(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc");
    std::string fname_snowage(
      "/Users/80x/Software/kernel_test_E3SM/pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc");

    const int n_procs = 1;
    const int ncells = 1;
    int idx = 0; // hardwire for ncells = 1
    const int ntimes = 3000;
    const int myrank = 0;
    const double dtime = 1800.0;
    const double dtime_d = 1800.0 / 86400.0;
    const auto start = ELM::Utils::Date(2014, 1, 1);

    auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
    auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { 1, 1 },
          { 0, 0 });
    
    // ELM State
    auto S = std::make_shared<ELMStateType>(ncells);
    
    S->Land.ltype = 1;
    S->Land.ctype = 1;
    S->Land.vtype = 12;
    S->Land.lakpoi = false;
    S->Land.urbpoi = false;

    assign(S->vtype, 12);
    auto host_vtype = Kokkos::create_mirror_view(S->vtype);

    // hardwired params
    S->lat = 71.323;
    S->lon = 203.3886;
    S->lat_r = S->lat * ELM::ELMconst::ELM_PI / 180.0;
    S->lon_r = S->lon * ELM::ELMconst::ELM_PI / 180.0;
    const double dewmx = 0.1;
    const double irrig_rate = 0.0;
    const int n_irrig_steps_left = 0;
    const int oldfflag = 1;
    //auto veg_active = create<ViewB1>("veg_active", ncells); // need value
    assign(S->veg_active, true);                               // hardwired
    //auto do_capsnow = create<ViewB1>("do_capsnow", ncells); // need value
    assign(S->do_capsnow, false);                               // hardwired
    assign(S->topo_slope, 0.070044865858546);
    assign(S->topo_std, 3.96141847422387);
    assign(S->snl, 0);
    assign(S->snow_depth, 0.0);
    assign(S->frac_sno, 0.0);
    assign(S->int_snow, 0.0);
    assign(S->h2osoi_liq, 0.0);
    assign(S->h2osoi_ice, 0.0);
    assign(S->t_h2osfc, 274.0);

    assign(S->eflx_sh_grnd, 0.0);
    assign(S->eflx_sh_snow, 0.0);
    assign(S->eflx_sh_soil, 0.0);
    assign(S->eflx_sh_h2osfc, 0.0);
    assign(S->qflx_evap_soi, 0.0);
    assign(S->qflx_ev_snow, 0.0);
    assign(S->qflx_ev_soil, 0.0);
    assign(S->qflx_ev_h2osfc, 0.0);

    assign(S->altmax_indx, 5);
    assign(S->altmax_lastyear_indx, 0);
    assign(S->t10, 276.0);
    assign(S->t_veg, 283.0);

    assign(S->xmf_dummy, 0.0);
    assign(S->xmf_h2osfc_dummy, 0.0);
    assign(S->eflx_h2osfc_snow_dummy, 0.0);
    assign(S->eflx_building_heat_dummy, 0.0);

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
      auto h_dz = Kokkos::create_mirror_view(S->dz);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
          h_dz(n, i) = dz_hardwire[i];
        }
      }
      Kokkos::deep_copy(S->dz, h_dz);

      double zsoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.007100635417193535,
      0.02792500041531687, 0.06225857393654604, 0.11886506690014327,
      0.21219339590896316, 0.3660657971047043, 0.6197584979298266,
      1.0380270500015696, 1.7276353086671965, 2.8646071131796917,
      4.73915671146575, 7.829766507142356, 12.92532061670855,
      21.32646906315379, 35.17762120511739 };
      auto h_zsoi = Kokkos::create_mirror_view(S->zsoi);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
          h_zsoi(n, i) = zsoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(S->zsoi, h_zsoi);

      double zisoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.017512817916255204, 0.04509178717593146, 0.09056182041834465, 
      0.16552923140455322, 0.28912959650683373, 0.4929121475172655,
      0.8288927739656982, 1.382831179334383, 2.2961212109234443,
      3.8018819123227208, 6.284461609304053, 10.377543561925453,
      17.12589483993117, 28.252045134135592, 42.10319727609919 };
      auto h_zisoi = Kokkos::create_mirror_view(S->zisoi);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd + 1; ++i) {
          h_zisoi(n, i) = zisoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(S->zisoi, h_zisoi);
    }

    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    // init data containers and read time-invariant data from files
    // these call only need to occur once @ beginning of simulation
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

    auto snicar_data = std::make_shared<ELM::SnicarData<ViewD1, ViewD2, ViewD3>>();
    auto snw_rds_table = std::make_shared<ELM::SnwRdsTable<ViewD3>>();
    auto pft_data = std::make_shared<ELM::PFTData<ViewD1, ViewD2>>();
    auto aerosol_data = std::make_shared<ELM::AerosolDataManager<ViewD1>>();
    

    ELM::initialize_kokkos_elm(S, snicar_data, snw_rds_table, pft_data,
                                aerosol_data, dd, fname_surfdata, fname_param,
                                fname_snicar, fname_snowage, fname_aerosol);

    int atm_nsteps = 101;
    const auto fstart = ELM::Utils::Date(1985, 1, 1);
    auto atm_forcing = std::make_shared<ELM::AtmForcObjects>(fname_forc, fstart, atm_nsteps, ncells);

    // phenology data manager
    // make host mirrors - need to be persistent
    ELM::PhenologyDataManager<ViewD2> phen_data(dd, ncells, 17);
    auto host_phen_views = get_phen_host_views(phen_data);

    // containers for aerosol deposition and concentration within snowpack layers
    auto aerosol_masses = std::make_shared<ELM::AerosolMasses<ViewD2>>(ncells);
    auto aerosol_concentrations = std::make_shared<ELM::AerosolConcentrations<ViewD2>>(ncells);

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

      auto h_soi_ice = Kokkos::create_mirror_view(S->h2osoi_ice);
      auto h_soi_liq = Kokkos::create_mirror_view(S->h2osoi_liq);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
          h_soi_ice(n, i) = h2osoi_ice_hardwire[i];
          h_soi_liq(n, i) = h2osoi_liq_hardwire[i];
        }
      }
      Kokkos::deep_copy(S->h2osoi_ice, h_soi_ice);
      Kokkos::deep_copy(S->h2osoi_liq, h_soi_liq);


      double h2osoi_vol_hardwire[] = {
        0.4016484663460637, 0.5196481455614503, 0.7967166638201649,
        0.8331813710901114, 0.7859200286330449, 0.7517405589446893,
        0.6621235242027332, 0.1535948180493002, 0.15947477948341815,
        0.15954052527228618, 8.420726808634413e-06, 5.107428986500891e-06,
        3.0978122726178113e-06, 1.8789181213767733e-06, 1.5092697845407248e-06 };
      auto h_soi_vol = Kokkos::create_mirror_view(S->h2osoi_vol);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevgrnd; ++i) {
          h_soi_vol(n, i) = h2osoi_vol_hardwire[i];
        }
      }
      Kokkos::deep_copy(S->h2osoi_vol, h_soi_vol);

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
      auto h_tsoi = Kokkos::create_mirror_view(S->t_soisno);
      for (int n = 0; n < ncells; ++n) {
        for (int i = 0; i < nlevsno + nlevgrnd; ++i) {
          h_tsoi(n, i) = tsoi_hardwire[i];
        }
      }
      auto h_tgrnd = Kokkos::create_mirror_view(S->t_grnd);
      h_tgrnd(idx) = h_tsoi(idx, nlevsno - S->snl(idx));
      Kokkos::deep_copy(S->t_soisno, h_tsoi);
      Kokkos::deep_copy(S->t_grnd, h_tgrnd);
    }




    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /*                          TIME LOOP                                                                  */
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
    ELM::Utils::Date current(start);

    for (int t = 0; t < ntimes; ++t) {

      ELM::Utils::Date time_plus_half_dt(current);
      time_plus_half_dt.increment_seconds(dtime/2);

      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // get coszen
      // only one value currently
      // will change when slope aspect modifier is completed
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      {
        // there are three methods to calculate zenith angle
        // they all produce similar results for the lat/lon tested here.

        // for now a single value for coszen is appropriate
        // but a slope based factor would necessitate per-cell values

        // first method - average cosz for dt_start to dt_end
        auto decday = ELM::Utils::decimal_doy(current) + 1.0;
        assign(S->coszen, ELM::incident_shortwave::average_cosz(S->lat_r, S->lon_r, dtime, decday));

        // second method - point cosz at dt_start + dt/2
        //auto thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, decday + dtime / 86400.0 /2.0);

        // third method - calc avg cosz over forcing dt (larger than model dt)
        // then calc point dt at start + dt/2
        // and use to calculate cosz_factor
        //ELM::Utils::Date forc_dt_start{forc_FSDS.get_data_start_time()};
        //forc_dt_start.increment_seconds(round(forc_FSDS.forc_t_idx(time_plus_half_dt, forc_FSDS.get_data_start_time()) * forc_FSDS.get_forc_dt_secs()));
        //double cosz_forc_decday = ELM::Utils::decimal_doy(forc_dt_start) + 1.0;
        //auto cosz_forcdt_avg = ELM::incident_shortwave::average_cosz(lat_r, lon_r, forc_FSDS.get_forc_dt_secs(), cosz_forc_decday);
        //auto thiscosz = ELM::incident_shortwave::coszen(lat_r, lon_r, decday + dtime / 86400.0 /2.0);
        //cosz_factor = (thiscosz > 0.001) ? std::min(thiscosz/cosz_forcdt_avg, 10.0) : 0.0;



        S->max_dayl = ELM::max_daylength(S->lat_r);
        S->dayl = ELM::daylength(S->lat_r, ELM::incident_shortwave::declination_angle2(current.doy + 1));
      }


      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // timestep init functions
      // read time-variable data
      // these are all self-invoking parallel kernels
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

      // read phenology data if required
      // reader will read 3 months of data on first call
      // subsequent calls only read the newest months (when phen_data.need_data() == true)
      // and shift the index of the two remaining older months
      update_phenology(phen_data, host_phen_views, S, host_vtype, current, fname_surfdata);

      // read new atm data if needed
      ELM::read_forcing(atm_forcing, dd, current, atm_nsteps);
      // get current time forcing values
      ELM::get_forcing(atm_forcing, S, dtime_d, time_plus_half_dt);

      // get aerosol mass (mss) and concentration in snowpack (cnc)
      ELM::aerosols::invoke_aerosol_source(time_plus_half_dt, dtime, S->snl, *aerosol_data, *aerosol_masses);
      ELM::aerosols::invoke_aerosol_concen_and_mass(dtime, S->do_capsnow, S->snl, S->h2osoi_liq,
        S->h2osoi_ice, S->snw_rds, S->qflx_snwcp_ice, *aerosol_masses, *aerosol_concentrations);

      // initialize a few variables at each dt
      // this function is pretty anemic
      // more should be incorporated into this call
      ELM::kokkos_init_timestep(S);


      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      // Main physics calls
      // single large loop to call all physics kernels
      // will move to less naive approach soon
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

      {
        // call surface albedo and SNICAR kernels
        ELM::kokkos_albedo_snicar(S, aerosol_concentrations, snicar_data, pft_data);

        // call canopy_hydrology kernels
        ELM::kokkos_canopy_hydrology(S, dtime);

        // call surface_radiation kernels
        ELM::kokkos_surface_radiation(S);

        // call canopy_temperature kernels
        ELM::kokkos_canopy_temperature(S, pft_data);

        // call bareground_fluxes kernels
        ELM::kokkos_bareground_fluxes(S);

        // call canopy_fluxes kernels
        ELM::kokkos_canopy_fluxes(S, dtime);

        // call soil_temperature kernels
        ELM::kokkos_soil_temperature(S, dtime);

        // call snow_hydrology kernels
        ELM::kokkos_snow_hydrology(S, aerosol_masses, aerosol_data, snw_rds_table, time_plus_half_dt, dtime);

        // call surface_fluxes kernels
        ELM::kokkos_surface_fluxes(S, dtime);
      }

  //    std::cout << "soil temp fluxes:  "
  //        << hs_soil(0) << "  "
  //        << hs_h2osfc(0) << "  "
  //        << hs_top(0) << "  "
  //        << hs_top_snow(0) << "  "
  //        << dhsdT(0) << "  "
  //        << soil_e_balance(0) << "  "
  //        << sabg_chk(0) << std::endl;

      std::cout << "lwrad_out:  "
          << S->eflx_lwrad_out(0) << "  "
          << S->eflx_lwrad_net(0) << std::endl;

      for (int i = 0; i < nlevsno + nlevgrnd; ++i)
        std::cout << "column vars:  " << i <<
         "  t_soisno:  " << S->t_soisno(0, i) <<
         "  h2osoi_ice:  " << S->h2osoi_ice(0, i) <<
         "  h2osoi_liq:  " << S->h2osoi_liq(0, i) <<
         "  dz:  " << S->dz(0, i) <<
         "  zsoi:  " << S->zsoi(0, i) <<
         "  zisoi:  " << S->zisoi(0, i) << std::endl;
         std::cout << "last   zisoi:  " << S->zisoi(0, nlevsno+nlevgrnd) << std::endl;




         for (int i = 0; i < nlevgrnd; ++i)
         std::cout << i <<"  qflx_rootsoi:  " << S->qflx_rootsoi(0, i) << "\n";


      for (int i = 0; i < ncells; ++i) {
        std::cout << "h2osno: " << S->h2osno(i) << std::endl;
        std::cout << "t_grnd: " << S->t_grnd(i) << std::endl;
        std::cout << "snow_depth: " << S->snow_depth(i) << std::endl;
        std::cout << "frac_sno: " << S->frac_sno(i) << std::endl;
        std::cout << "frac_sno_eff: " << S->frac_sno_eff(i) << std::endl;
        std::cout << "qflx_tran_veg: " << S->qflx_tran_veg(i) << std::endl;
        std::cout << "frac_veg_nosno_alb: " << S->frac_veg_nosno_alb(i) << std::endl;
      }

      current.increment_seconds(dtime);

    } // time loop

  } // inner scope

  Kokkos::finalize();
} // enclosing scope
return 0;
}
