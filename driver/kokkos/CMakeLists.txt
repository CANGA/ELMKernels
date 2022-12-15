project(ELM_KERNELS_KOKKOS)

set(kokkos_driver_src_files
aerosol_kokkos.cc
atm_forcing_kokkos.cc
phenology_kokkos.cc
initialize_elm_kokkos.cc
albedo_kokkos.cc
canopy_hydrology_kokkos.cc
conserved_quantity_kokkos.cc
surface_radiation_kokkos.cc
canopy_temperature_kokkos.cc
bareground_fluxes_kokkos.cc
canopy_fluxes_kokkos.cc
soil_temperature_kokkos.cc
snow_hydrology_kokkos.cc
surface_fluxes_kokkos.cc
init_timestep_kokkos.cc
elm_kokkos_interface.cc
#kokkos_driver.cc
)

set(kokkos_driver_inc_files
aerosol_kokkos.hh
atm_forcing_kokkos.hh
phenology_kokkos.hh
initialize_elm_kokkos.hh
albedo_kokkos.hh
canopy_hydrology_kokkos.hh
conserved_quantity_kokkos.hh
surface_radiation_kokkos.hh
canopy_temperature_kokkos.hh
bareground_fluxes_kokkos.hh
canopy_fluxes_kokkos.hh
soil_temperature_kokkos.hh
snow_hydrology_kokkos.hh
surface_fluxes_kokkos.hh
init_timestep_kokkos.hh
elm_kokkos_interface.hh
)

add_library (elm_kokkos
${kokkos_driver_src_files}
${kokkos_driver_inc_files}
)

add_compile_definitions(INPUT_DATA_DIR="${ELM_INPUT_DATA_PREFIX}/")

target_include_directories (elm_kokkos PUBLIC ${ELM_KERNELS_KOKKOS_SOURCE_DIR})
target_link_libraries (elm_kokkos LINK_PUBLIC elm_physics)

add_executable (elm_kokkos_exe
kokkos_driver.cc
)

target_include_directories (elm_kokkos_exe PUBLIC ${ELM_KERNELS_KOKKOS_SOURCE_DIR})
target_link_libraries (elm_kokkos_exe LINK_PUBLIC elm_kokkos)

install(TARGETS elm_kokkos elm_kokkos_exe)
install(FILES ${kokkos_driver_inc_files} DESTINATION include)
