
#include "init_timestep_kokkos.hh"
#include "atm_forcing_kokkos.hh"
#include "phenology_kokkos.hh"
#include "aerosol_kokkos.hh"
#include "helper_functions.hh"
#include "invoke_kernel.hh"

#include "day_length.h"
#include "incident_shortwave.h"
#include "init_timestep.h"

#include "conserved_quantity_evaluators.h"




void ELM::kokkos_init_timestep(ELMStateType *S,
                              const double dtime,
                              const Utils::Date& current,
                              const std::string& fname_surfdata)
{

double dtime_d = dtime / 86400.0;
Utils::Date t_centered(current);
t_centered.increment_seconds(static_cast<size_t>((dtime + 1.0)/2)); // round to nearest second

  // get coszen
  // only one value currently
  // will change when slope aspect modifier is completed
  auto decday = Utils::decimal_doy(current) + 1.0;
  double cosz = incident_shortwave::average_cosz(S->lat_r, S->lon_r, dtime, decday);
  Utils::assign(S->coszen, cosz);
  S->max_dayl = ELM::max_daylength(S->lat_r);
  S->dayl = ELM::daylength(
    S->lat_r, ELM::incident_shortwave::declination_angle_sin(current.doy + 1));

  // read phenology data if required
  // reader will read 3 months of data on first call
  // subsequent calls only read the newest months (when phen_data.need_data() == true)
  // and shift the index of the two remaining older months
  ELM::update_phenology(*S, current, fname_surfdata);

  // read new atm data if needed
  // get current time forcing values
  ELM::read_forcing(*S, current);
  ELM::get_forcing(*S, dtime_d, t_centered);

  // get aerosol mass (mss) and concentration in snowpack (cnc)
  ELM::invoke_aerosol_source(*S, dtime, t_centered);
  ELM::invoke_aerosol_concen_and_mass(*S, dtime);

  auto init_step_kernel = ELM_LAMBDA (const int& idx) {

    S->h2osno_old(idx) = S->h2osno(idx);

    S->dtbegin_column_h2o(idx) =
      ELM::conservation_eval::column_water_mass(S->h2ocan(idx),
        S->h2osno(idx), S->h2osfc(idx),
        Kokkos::subview(S->h2osoi_ice, idx, Kokkos::ALL),
        Kokkos::subview(S->h2osoi_liq, idx, Kokkos::ALL));

    ELM::init_timestep(S->Land.lakpoi, S->veg_active(idx),
                     S->frac_veg_nosno_alb(idx),
                     S->snl(idx), S->h2osno(idx),
                     Kokkos::subview(S->h2osoi_ice, idx, Kokkos::ALL),
                     Kokkos::subview(S->h2osoi_liq, idx, Kokkos::ALL),
                     S->do_capsnow(idx),
                     S->frac_veg_nosno(idx),
                     Kokkos::subview(S->frac_iceold, idx, Kokkos::ALL));
  }; // end init_step lambda
  apply_parallel_for(init_step_kernel, "kokkos_init_timestep", S->snl.extent(0));
}





