cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "\n-------------------------------------------")
message(STATUS "\n-- CMake: Configuring and building ELMKernels.\n--")
message(STATUS "----------------------------------------")
project(ELMKernels)

enable_testing()

#project(ELM CXX)
#enable_language(C)

# need netcdf
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/config/cmake")
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

# options for drivers
# only build one
option(ENABLE_KOKKOS "Enable building with Kokkos driver" OFF)
if (ENABLE_KOKKOS)
  add_compile_definitions(ENABLE_KOKKOS)
  find_package(Kokkos REQUIRED)
  message(STATUS "Building with Kokkos driver")
  option(ENABLE_CC "Enable building with default C-style driver" OFF)
else()
  option(ENABLE_CC "Enable building with default C-style driver" ON)
  message(STATUS "Building with CC driver")
endif()

add_subdirectory (src)
add_subdirectory (driver)

#add_subdirectory (test)
#add_test (NAME CanopyHydrology COMMAND test_CanHydro)
