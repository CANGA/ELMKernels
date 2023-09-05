
#include "elm_kokkos_interface.hh"

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
#include "conserved_quantity_kokkos.hh"

// constants
#include "elm_constants.h"

// conditional compilation options
#include "compile_options.hh"

// utilities
#include "read_input.hh"
#include "read_netcdf.hh"
#include "helper_functions.hh"



// initialization routines
#include "initialize_elm_kokkos.hh"
#include "init_timestep_kokkos.hh"

using namespace ELM::ELMdims;

namespace ELM {

// cell_per_col_ fixed as nlevgrnd() for now
ELMInterface::ELMInterface(size_t ncols) : ncols_(ncols), cell_per_col_(nlevgrnd())
{
  // temp dummy filepath
  std::string input_dir = INPUT_DATA_DIR;
  fname_surfdata_ = input_dir +
    "E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc";
  fname_forc_ = input_dir +
    "pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc";

  // allocate ELM State
  // TODO decouple notions of domain decomp, input file reads,
  // dates from state - remove everything aside from ncols_
  S_ = std::make_shared<ELMStateType>(ncols_,
    ELM::Utils::create_domain_decomposition_2D(
    ELM::Utils::square_numprocs(1), { 1, 1 },{ 0, 0 }), // domain decomp
    fname_forc_,
    ELM::Utils::Date(1985, 1, 1),
    101);
}

void ELMInterface::setup()
{
  using ELM::Utils::assign;

  // temporary - ATS will provide these 
  // the path/to/file for phenology and atm forcing inputs
  std::string input_dir = INPUT_DATA_DIR;
  // forcing files
  std::string fname_snicar(
    input_dir+"pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc");
  std::string fname_param(
    input_dir+"E3SM/components/elm/test_submodules/inputdata/lnd/clm2/paramdata/clm_params_c180524.nc");
  std::string fname_aerosol(
    input_dir+"pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc");
  std::string fname_snowage(
    input_dir+"pt-e3sm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc");



/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/*                          TEMPORARY - Provide some initial values here                               */
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/


      // fix this soon
    //S_->Land.ltype = 1;
    S_->Land.ltype = 1;
    S_->Land.ctype = 1;
    S_->Land.vtype = 12;
    S_->Land.lakpoi = false;
    S_->Land.urbpoi = false;

    assign(S_->vtype, 12);

    // hardwired params
    S_->lat = 71.323;
    S_->lon = 203.3886;
    S_->lat_r = S_->lat * ELM::ELMconst::ELM_PI() / 180.0;
    S_->lon_r = S_->lon * ELM::ELMconst::ELM_PI() / 180.0;
    const double dewmx = 0.1;
    const double irrig_rate = 0.0;
    const int n_irrig_steps_left = 0;
    const int oldfflag = 1;
    //auto veg_active = create<ViewB1>("veg_active", ncols); // need value
    assign(S_->veg_active, true);                               // hardwired
    //auto do_capsnow = create<ViewB1>("do_capsnow", ncols); // need value
    assign(S_->do_capsnow, false);                               // hardwired
    assign(S_->topo_slope, 0.070044865858546);
    assign(S_->topo_std, 3.96141847422387);
    assign(S_->snl, 0);
    assign(S_->snow_depth, 0.0);
    assign(S_->frac_sno, 0.0);
    assign(S_->int_snow, 0.0);
    assign(S_->h2osoi_liq, 0.0);
    assign(S_->h2osoi_ice, 0.0);
    assign(S_->t_h2osfc, 274.0);

    assign(S_->eflx_sh_grnd, 0.0);
    assign(S_->eflx_sh_snow, 0.0);
    assign(S_->eflx_sh_soil, 0.0);
    assign(S_->eflx_sh_h2osfc, 0.0);
    assign(S_->qflx_evap_soi, 0.0);
    assign(S_->qflx_ev_snow, 0.0);
    assign(S_->qflx_ev_soil, 0.0);
    assign(S_->qflx_ev_h2osfc, 0.0);

    assign(S_->altmax_indx, 5);
    assign(S_->altmax_lastyear_indx, 0);
    assign(S_->t10, 276.0);
    assign(S_->t_veg, 283.0);

    assign(S_->xmf, 0.0);
    assign(S_->xmf_h2osfc, 0.0);
    assign(S_->eflx_h2osfc_snow, 0.0);

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
      auto h_dz = Kokkos::create_mirror_view(S_->dz);
      for (int n = 0; n < ncols_; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_dz(n, i) = dz_hardwire[i];
        }
      }
      Kokkos::deep_copy(S_->dz, h_dz);

      double zsoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.007100635417193535,
      0.02792500041531687, 0.06225857393654604, 0.11886506690014327,
      0.21219339590896316, 0.3660657971047043, 0.6197584979298266,
      1.0380270500015696, 1.7276353086671965, 2.8646071131796917,
      4.73915671146575, 7.829766507142356, 12.92532061670855,
      21.32646906315379, 35.17762120511739 };
      auto h_zsoi = Kokkos::create_mirror_view(S_->zsoi);
      for (int n = 0; n < ncols_; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_zsoi(n, i) = zsoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(S_->zsoi, h_zsoi);

      double zisoi_hardwire[] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.017512817916255204, 0.04509178717593146, 0.09056182041834465, 
      0.16552923140455322, 0.28912959650683373, 0.4929121475172655,
      0.8288927739656982, 1.382831179334383, 2.2961212109234443,
      3.8018819123227208, 6.284461609304053, 10.377543561925453,
      17.12589483993117, 28.252045134135592, 42.10319727609919 };
      auto h_zisoi = Kokkos::create_mirror_view(S_->zisoi);
      for (int n = 0; n < ncols_; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd() + 1; ++i) {
          h_zisoi(n, i) = zisoi_hardwire[i];
        }
      }
      Kokkos::deep_copy(S_->zisoi, h_zisoi);
    }

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

      auto h_soi_ice = Kokkos::create_mirror_view(S_->h2osoi_ice);
      auto h_soi_liq = Kokkos::create_mirror_view(S_->h2osoi_liq);
      for (int n = 0; n < ncols_; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_soi_ice(n, i) = h2osoi_ice_hardwire[i];
          h_soi_liq(n, i) = h2osoi_liq_hardwire[i];
        }
      }
      Kokkos::deep_copy(S_->h2osoi_ice, h_soi_ice);
      Kokkos::deep_copy(S_->h2osoi_liq, h_soi_liq);


      double h2osoi_vol_hardwire[] = {
        0.4016484663460637, 0.5196481455614503, 0.7967166638201649,
        0.8331813710901114, 0.7859200286330449, 0.7517405589446893,
        0.6621235242027332, 0.1535948180493002, 0.15947477948341815,
        0.15954052527228618, 8.420726808634413e-06, 5.107428986500891e-06,
        3.0978122726178113e-06, 1.8789181213767733e-06, 1.5092697845407248e-06 };
      auto h_soi_vol = Kokkos::create_mirror_view(S_->h2osoi_vol);
      for (int n = 0; n < ncols_; ++n) {
        for (int i = 0; i < nlevgrnd(); ++i) {
          h_soi_vol(n, i) = h2osoi_vol_hardwire[i];
        }
      }
      Kokkos::deep_copy(S_->h2osoi_vol, h_soi_vol);

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
      auto h_tsoi = Kokkos::create_mirror_view(S_->t_soisno);
      for (int n = 0; n < ncols_; ++n) {
        for (int i = 0; i < nlevsno() + nlevgrnd(); ++i) {
          h_tsoi(n, i) = tsoi_hardwire[i];
        }
      }
      int idx = 0; // hardwire for ncols = 1
      auto h_tgrnd = Kokkos::create_mirror_view(S_->t_grnd);
      h_tgrnd(idx) = h_tsoi(idx, nlevsno() - S_->snl(idx));
      Kokkos::deep_copy(S_->t_soisno, h_tsoi);
      Kokkos::deep_copy(S_->t_grnd, h_tgrnd);
    }


/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/*    END TEMPORARY                         */
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

    ELM::initialize_kokkos_elm(*S_, fname_surfdata_,
      fname_param, fname_snicar, fname_snowage, fname_aerosol);

}


bool ELMInterface::advance(const Utils::Date& dt_start_date, const double& dt_seconds)
{
  ELM::Utils::Date time_plus_half_dt(dt_start_date);
  time_plus_half_dt.increment_seconds(static_cast<size_t>((dt_seconds + 1.0)/2)); // round to nearest second

  // get coszen, day length,
  // phenology data, atmospheric forcing,
  // aerosol forcing and snowpack state,
  // and a handful of variables that need to be reset every timestep
  ELM::kokkos_init_timestep(*S_, dt_seconds, dt_start_date, fname_surfdata_);


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // Main physics calls
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

  {
    // canhydro::fraction_wet
    ELM::kokkos_frac_wet(*S_);

    // call surface albedo and SNICAR kernels
    ELM::kokkos_albedo_snicar(*S_);

    // call canopy_hydrology kernels
    ELM::kokkos_canopy_hydrology(*S_, dt_seconds);

    // call surface_radiation kernels
    ELM::kokkos_surface_radiation(*S_);

    // call canopy_temperature kernels
    ELM::kokkos_canopy_temperature(*S_);

    // call bareground_fluxes kernels
    ELM::kokkos_bareground_fluxes(*S_);

    // call canopy_fluxes kernels
    ELM::kokkos_canopy_fluxes(*S_, dt_seconds);

    // call soil_temperature kernels
    ELM::kokkos_soil_temperature(*S_, dt_seconds);

    // call snow_hydrology kernels
    ELM::kokkos_snow_hydrology(*S_, dt_seconds, time_plus_half_dt);

    // call surface_fluxes kernels
    ELM::kokkos_surface_fluxes(*S_, dt_seconds);

    ELM::kokkos_evaluate_conservation(*S_, dt_seconds);
  }

  return false;
}

void
ELMInterface::copyPrimaryVars (
PrimaryVars<ViewI1, ViewD1, ViewD2>& primary_vars)
{
  Kokkos::deep_copy (primary_vars.snl, S_->snl );
  Kokkos::deep_copy (primary_vars.snow_depth, S_->snow_depth );
  Kokkos::deep_copy (primary_vars.frac_sno, S_->frac_sno );
  Kokkos::deep_copy (primary_vars.int_snow, S_->int_snow );
  Kokkos::deep_copy (primary_vars.snw_rds, S_->snw_rds );
  Kokkos::deep_copy (primary_vars.h2osoi_liq, S_->h2osoi_liq );
  Kokkos::deep_copy (primary_vars.h2osoi_ice, S_->h2osoi_ice );
  Kokkos::deep_copy (primary_vars.h2osoi_vol, S_->h2osoi_vol );
  Kokkos::deep_copy (primary_vars.h2ocan, S_->h2ocan );
  Kokkos::deep_copy (primary_vars.h2osno, S_->h2osno );
  Kokkos::deep_copy (primary_vars.h2osfc, S_->h2osfc );
  Kokkos::deep_copy (primary_vars.t_soisno, S_->t_soisno );
  Kokkos::deep_copy (primary_vars.t_grnd, S_->t_grnd );
  Kokkos::deep_copy (primary_vars.t_h2osfc, S_->t_h2osfc );
  Kokkos::deep_copy (primary_vars.t_h2osfc_bef, S_->t_h2osfc_bef );
  Kokkos::deep_copy (primary_vars.nrad, S_->nrad );
  Kokkos::deep_copy (primary_vars.dz, S_->dz );
  Kokkos::deep_copy (primary_vars.zsoi, S_->zsoi );
  Kokkos::deep_copy (primary_vars.zisoi, S_->zisoi );
}

std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> >
ELMInterface::getPrimaryVars()
{
  auto primary_vars =
    std::make_shared<PrimaryVars<ViewI1, ViewD1, ViewD2> > (ncols_);
  copyPrimaryVars(*primary_vars);
  return primary_vars;
}

} // namespace ELM
