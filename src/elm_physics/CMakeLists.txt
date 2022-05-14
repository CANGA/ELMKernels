project(ELM_PHYSICS)

add_library (elm_physics
#friction_velocity.cc
#qsat.cc
#bareground_fluxes.cc
#canopy_hydrology.cc
#canopy_temperature.cc
#init_topography.cc
#photosynthesis.cc
#soil_moist_stress.cc
#surface_albedo.cc
day_length.cc
incident_shortwave.cc
#atm_physics.cc
#snicar_data.cc
monthly_data.cc
aerosol_data.cc
#phenology_data.cc
#surface_fluxes.cc
#soil_texture_hydraulic_model.cc
)

#include_directories(${ELM_PHYSICS_SOURCE_DIR})
#include_directories(${ELM_UTILS_SOURCE_DIR})

if (ENABLE_KOKKOS)
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR} ${KOKKOS_INCLUDE_DIR})
target_link_libraries (elm_physics LINK_PUBLIC elm_utils Kokkos::kokkos)
else()
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})
target_link_libraries (elm_physics LINK_PUBLIC elm_utils)
endif()



install(TARGETS elm_physics)
