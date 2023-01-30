
#include "date_time.hh"
#include "aerosol_data.h"
#include "aerosol_physics.h"

#include "invoke_kernel.hh"
#include "aerosol_kokkos.hh"

void ELM::invoke_aerosol_source(ELMStateType& S,
                                const double& dtime,
                                const ELM::Utils::Date& model_time)
{
  auto aerosol_forc_flux = 
    S.aerosol_data->get_aerosol_source(model_time, dtime);
  ELM::aerosols::ComputeAerosolDeposition
    aerosol_source_object(aerosol_forc_flux,
                          S.snl, *S.aerosol_masses);

  apply_parallel_for(aerosol_source_object, "ComputeAerosolDeposition", S.num_columns);
}


void ELM::invoke_aerosol_concen_and_mass(ELMStateType& S, 
                                         const double& dtime)
{
  ELM::aerosols::ComputeAerosolConcenAndMass
    aerosol_c_mass_object(dtime, S.do_capsnow, S.snl, S.h2osoi_liq,
                          S.h2osoi_ice, S.snw_rds, S.qflx_snwcp_ice,
                          *S.aerosol_masses, *S.aerosol_concentrations);

  apply_parallel_for(aerosol_c_mass_object, "ComputeAerosolConcenAndMass", S.num_columns);
}
