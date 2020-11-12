add_library (physics physics_dummy.cc fracwet.cc canopy_interception.cc canopy_ground_flux.cc canopy_irrigation.cc)
target_include_directories (physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
