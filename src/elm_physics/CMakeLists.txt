project(ELM_PHYSICS)

set(physics_src_files
day_length.cc
incident_shortwave.cc
)

set(physics_inc_files
aerosol_physics.h
aerosol_physics_impl.hh
atm_physics.h
atm_physics_impl.hh
bareground_fluxes.h
bareground_fluxes_impl.hh
canopy_fluxes.h
canopy_fluxes_impl.hh
canopy_hydrology.h
canopy_hydrology_impl.hh
canopy_temperature.h
canopy_temperature_impl.hh
day_length.h
elm_constants.h
friction_velocity.h
friction_velocity_impl.hh
incident_shortwave.h
init_snow_state.h
init_snow_state_impl.hh
init_soil_state.h
init_soil_state_impl.hh
init_timestep.h
init_timestep_impl.hh
init_topography.h
init_topography_impl.hh
pentadiagonal_solver.h
pentadiagonal_solver_impl.hh
phase_change.h
phase_change_impl.hh
phenology_physics.h
phenology_physics_impl.hh
photosynthesis.h
photosynthesis_impl.hh
qsat.h
qsat_impl.hh
snow_hydrology.h
snow_hydrology_impl.hh
snow_snicar.h
snow_snicar_impl.hh
soil_moist_stress.h
soil_moist_stress_impl.hh
soil_temp_lhs.h
soil_temp_lhs_impl.hh
soil_temp_rhs.h
soil_temp_rhs_impl.hh
soil_temperature.h
soil_temperature_impl.hh
soil_texture_hydraulic_model.h
soil_texture_hydraulic_model_impl.hh
soil_thermal_properties.h
soil_thermal_properties_impl.hh
surface_albedo.h
surface_albedo_impl.hh
surface_fluxes.h
surface_fluxes_impl.hh
surface_radiation.h
surface_radiation_impl.hh
surface_resistance.h
surface_resistance_impl.hh
transpiration.h
transpiration_impl.hh
)

add_library (elm_physics
${physics_src_files}
${physics_inc_files}
)

if (ENABLE_KOKKOS)
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_INPUT_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR} ${KOKKOS_INCLUDE_DIR})
target_link_libraries (elm_physics LINK_PUBLIC elm_input elm_utils Kokkos::kokkos)
else()
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_INPUT_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})
target_link_libraries (elm_physics LINK_PUBLIC elm_input elm_utils)
endif()

install(TARGETS elm_physics)
