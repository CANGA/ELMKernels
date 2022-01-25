project(ELM_PHYSICS)

add_library (elm_physics friction_velocity.cc qsat.cc bareground_fluxes.cc canopy_hydrology.cc canopy_temperature.cc InitTopography.cc photosynthesis.cc soil_moist_stress.cc 
#surface_albedo.cc CanopyHydrology_newphysics.cc CanHydro_newphys_internal_logic.cc CanHydro_newphys_external_logic.cc
surface_albedo.cc canopy_hydrology.cc
daylength.cc incident_shortwave.cc atm_forcing_physics.cc snicar_data.cc monthly_data.cc aerosol_data.cc)
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})

target_link_libraries (elm_physics LINK_PUBLIC elm_utils)
install(TARGETS elm_physics)
