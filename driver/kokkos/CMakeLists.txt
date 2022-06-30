
set(kokkos_driver_src_files
atm_helpers.cc
kokkos_elm_initialize.cc
albedo_kokkos.cc
canopy_hydrology_kokkos.cc
surface_radiation_kokkos.cc
canopy_temperature_kokkos.cc
bareground_fluxes_kokkos.cc
canopy_fluxes_kokkos.cc
soil_temperature_kokkos.cc
new_kokkos_driver.cc
)

set(kokkos_driver_inc_files
atm_helpers.hh
kokkos_elm_initialize.hh
albedo_kokkos.hh
canopy_hydrology_kokkos.hh
surface_radiation_kokkos.hh
canopy_temperature_kokkos.hh
bareground_fluxes_kokkos.hh
canopy_fluxes_kokkos.hh
soil_temperature_kokkos.hh
)

add_executable (elm_kokkos
${kokkos_driver_src_files}
${kokkos_driver_inc_files}
)

target_include_directories (elm_kokkos PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_INPUT_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR} ${KOKKOS_INCLUDE_DIR} ${NetCDF_INCLUDE_DIR})
target_link_libraries (elm_kokkos LINK_PUBLIC elm_physics elm_input elm_utils Kokkos::kokkos ${NetCDF_LIBRARIES})


install(TARGETS elm_kokkos)
