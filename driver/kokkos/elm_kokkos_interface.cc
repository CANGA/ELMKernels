
#include "elm_kokkos_interface.hh"

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

  // initialize data containers and read time-invariant data from files
  // these calls only need to occur once @ beginning of simulation
  ELM::initialize_kokkos_elm(*S_.get(), fname_surfdata_,
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
  ELM::kokkos_init_timestep(S_.get(), dt_seconds, dt_start_date, fname_surfdata_);


  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // Main physics calls
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

  {
    // canhydro::fraction_wet
    ELM::kokkos_frac_wet(*S_.get());

    // call surface albedo and SNICAR kernels
    ELM::kokkos_albedo_snicar(*S_.get());

    // call canopy_hydrology kernels
    ELM::kokkos_canopy_hydrology(*S_.get(), dt_seconds);

    // call surface_radiation kernels
    ELM::kokkos_surface_radiation(*S_.get());

    // call canopy_temperature kernels
    ELM::kokkos_canopy_temperature(*S_.get());

    // call bareground_fluxes kernels
    ELM::kokkos_bareground_fluxes(*S_.get());

    // call canopy_fluxes kernels
    ELM::kokkos_canopy_fluxes(*S_.get(), dt_seconds);

    // call soil_temperature kernels
    ELM::kokkos_soil_temperature(*S_.get(), dt_seconds);

    // call snow_hydrology kernels
    ELM::kokkos_snow_hydrology(*S_.get(), dt_seconds, time_plus_half_dt);

    // call surface_fluxes kernels
    ELM::kokkos_surface_fluxes(*S_.get(), dt_seconds);

    ELM::kokkos_evaluate_conservation(*S_.get(), dt_seconds);
  }

  return false;
}

void
ELMInterface::copyPrimaryVars (
std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> > primary_vars)
{
  Kokkos::deep_copy (primary_vars->snl, S_->snl );
  Kokkos::deep_copy (primary_vars->snow_depth, S_->snow_depth );
  Kokkos::deep_copy (primary_vars->frac_sno, S_->frac_sno );
  Kokkos::deep_copy (primary_vars->int_snow, S_->int_snow );
  Kokkos::deep_copy (primary_vars->snw_rds, S_->snw_rds );
  Kokkos::deep_copy (primary_vars->h2osoi_liq, S_->h2osoi_liq );
  Kokkos::deep_copy (primary_vars->h2osoi_ice, S_->h2osoi_ice );
  Kokkos::deep_copy (primary_vars->h2osoi_vol, S_->h2osoi_vol );
  Kokkos::deep_copy (primary_vars->h2ocan, S_->h2ocan );
  Kokkos::deep_copy (primary_vars->h2osno, S_->h2osno );
  Kokkos::deep_copy (primary_vars->h2osfc, S_->h2osfc );
  Kokkos::deep_copy (primary_vars->t_soisno, S_->t_soisno );
  Kokkos::deep_copy (primary_vars->t_grnd, S_->t_grnd );
  Kokkos::deep_copy (primary_vars->t_h2osfc, S_->t_h2osfc );
  Kokkos::deep_copy (primary_vars->t_h2osfc_bef, S_->t_h2osfc_bef );
  Kokkos::deep_copy (primary_vars->nrad, S_->nrad );
  Kokkos::deep_copy (primary_vars->dz, S_->dz );
  Kokkos::deep_copy (primary_vars->zsoi, S_->zsoi );
  Kokkos::deep_copy (primary_vars->zisoi, S_->zisoi );
}

std::shared_ptr<PrimaryVars<ViewI1, ViewD1, ViewD2> >
ELMInterface::getPrimaryVars()
{
  auto primary_vars =
    std::make_shared<PrimaryVars<ViewI1, ViewD1, ViewD2> > (ncols_);
  copyPrimaryVars(primary_vars);
  return primary_vars;
}

} // namespace ELM
