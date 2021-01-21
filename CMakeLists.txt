cmake_minimum_required(VERSION 3.13.3)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

project(ELM_physics CXX)

# need netcdf
list(APPEND CMAKE_MODULE_PATH "${ELM_physics_SOURCE_DIR}/config/cmake")
include(AddImportedLibrary)

if (BUILD_SHARED_LIBS)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_SKIP_INSTALL_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) 
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_RPATH}")
else()
  SET(CMAKE_SKIP_INSTALL_RPATH ON CACHE BOOL "Turn off for static install." FORCE)
  SET(CMAKE_SKIP_RPATH ON CACHE BOOL "Turn off for static install." FORCE)
endif()


find_package(NetCDF REQUIRED)
message("Found NetCDF = ${NetCDF_C_LIBRARIES}")

add_subdirectory (src)
add_subdirectory (driver)

