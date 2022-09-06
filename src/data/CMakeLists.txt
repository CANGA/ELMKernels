project(ELM_DATA)

set(data_src_files
aerosol_data.cc
monthly_data.cc
)

set(data_inc_files
aerosol_data.h
aerosol_data_impl.hh
atm_data.h
atm_data_impl.hh
elm_state.h
elm_state_impl.hh
land_data.h
monthly_data.h
pft_data.h
pft_data_impl.hh
phenology_data.h
phenology_data_impl.hh
snicar_data.h
snicar_data_impl.hh
soil_data.h
soil_data_impl.hh
)

add_library (elm_data
${data_src_files}
${data_inc_files}
)

if (ENABLE_KOKKOS)
target_include_directories (elm_data PUBLIC ${ELM_DATA_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR} ${KOKKOS_INCLUDE_DIR})
target_link_libraries (elm_data LINK_PUBLIC elm_utils Kokkos::kokkos)
else()
target_include_directories (elm_data PUBLIC ${ELM_DATA_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})
target_link_libraries (elm_data LINK_PUBLIC elm_utils)
endif()

install(TARGETS elm_data)
