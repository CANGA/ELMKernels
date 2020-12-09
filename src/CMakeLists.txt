add_library (physics physics_dummy.cc fracwet.cc canopy_interception.cc canopy_ground_flux.cc canopy_irrigation.cc coldstart_topography.cc coldstart_snow.cc snow_init.cc frac_h2osfc.cc canopy_shade_frac.cc SurfaceRadiationMod_functions.cc)
target_include_directories (physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
