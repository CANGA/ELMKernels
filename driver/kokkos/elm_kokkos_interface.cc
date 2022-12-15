
#include "elm_kokkos_interface.hh"

using namespace ELM::ELMdims;

ELM::ELMKernels::ELMKernels(size_t ncols) :
  ncols_(ncols), cell_per_col_(nlevgrnd)
{
  ELM::ELMKernels::init();
}

void ELM::ELMKernels::init()
{
  // hardwired
  const int n_procs = 1;
  auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
  // domain and processor topology info
  auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
        { 1, 1 },
        { 0, 0 });

  // number of forcing timesteps to store in device memory
  int atm_nsteps = 101;
  // starting time of forcing file
  const auto fstart = ELM::Utils::Date(1985, 1, 1);

  // the path/to/file for phenology and atm forcing inputs
  std::string input_dir = INPUT_DATA_DIR;
  fname_surfdata_ = input_dir +
    "E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc";
  fname_forc_ = input_dir +
    "pt-e3sm-inputdata/atm/datm7/1x1pt_US-Brw/cpl_bypass_full/all_hourly.nc";

  // ELM State
  S_ = std::make_shared<ELMStateType>(ncols_, dd, fname_forc_, fstart, atm_nsteps);

}

void ELM::ELMKernels::setup()
{
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


bool ELM::ELMKernels::advance(const Utils::Date& dt_start_date, const double& dt_seconds)
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
