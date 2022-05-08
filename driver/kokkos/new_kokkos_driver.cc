
#include <iostream>
#include <string>

// utilities
#include "array.hh"
#include "read_input.hh"
#include "read_netcdf.hh"
#include "utils.hh"

// constants
#include "elm_constants.h"

// input data readers and structs
#include "pft_data.h"
#include "atm_data.h"
#include "soil_data.h"
#include "snicar_data.h"
//#include "satellite_phenology.h"

//#include "InitSoil.hh"
#include "init_soil_state.h"
#include "init_snow_state.h"


//#include "InitSnowLayers.hh"
//#include "InitTimestep.hh"
#include "init_timestep.h"
//#include "InitTopography.hh"
#include "init_topography.h"
#include "land_data.h"
//#include "read_atmosphere.h"
//#include "ReadTestData.hh"

#include "incident_shortwave.h"
#include "day_length.h"

#include "canopy_hydrology.h"
#include "surface_radiation.h"
#include "canopy_temperature.h"
#include "bareground_fluxes.h"
#include "canopy_fluxes.h"
#include "aerosol_data.h"
#include "aerosol_physics.h"
#include "phenology_data.h"
#include "surface_albedo.h"
#include "snow_snicar.h"
//#include "root_biophys.h"
#include "surface_fluxes.h"
#include "soil_texture_hydraulic_model.h"

#include "Kokkos_Core.hpp"

using ArrayB1 = Kokkos::View<bool *>;
using ArrayI1 = Kokkos::View<int *>;
using ArrayD1 = Kokkos::View<double *>;
using ArrayD2 = Kokkos::View<double **>;
using ArrayP1 = Kokkos::View<ELM::PSNVegData *>;

template <class Array_t> Array_t create(const std::string &name, int D0)
{ return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1)
{ return Array_t(name, D0, D1); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2)
{ return Array_t(name, D0, D1, D2); }
template <class Array_t, typename T> void assign(Array_t &arr, T val)
{ Kokkos::deep_copy(arr, val); }


using AtmForcType = ELM::AtmForcType;

template<AtmForcType ftype>
using atm_forc_util = ELM::AtmDataManager<ArrayD1, ArrayD2, ftype>;

template <AtmForcType ftype>
atm_forc_util<ftype> create_forc_util(const std::string& filename,
  const ELM::Utils::Date &file_start_time, const int ntimes,
  const int ncells)
{ return atm_forc_util<ftype>(filename, file_start_time, ntimes, ncells); }

int main(int argc, char **argv) {

{ // enclosing scope

  Kokkos::initialize(argc, argv);

  { // inner scope

    int MPI_COMM_WORLD;
    const int n_procs = 1;
    const int ncells = 1;
    const int ntimes = 1008;
    const int myrank = 0;
    const double dtime = 1800.0;
    const double dtime_d = 1800.0 / 86400.0;
    const auto start = ELM::Utils::Date(2014, 7, 15);
    ELM::LandType Land;
    Land.ltype = 1;
    Land.ctype = 1;
    Land.vtype = 12;

    std::string fname_surfdata(
    "/Users/80x/Software/kernel_test_E3SM/E3SM/components/elm/test_submodules/inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_forcanga_arcticgrass.nc");

    auto proc_decomp = ELM::Utils::square_numprocs(n_procs);
    auto dd = ELM::Utils::create_domain_decomposition_2D(proc_decomp,
          { 1, 1 },
          { 0, 0 });

    // forcing data
    auto forc_tbot = create<ArrayD1>("forc_tbot", ncells);
    auto forc_thbot = create<ArrayD1>("forc_thbot", ncells);
    auto forc_pbot = create<ArrayD1>("forc_pbot", ncells);
    auto forc_qbot = create<ArrayD1>("forc_qbot", ncells);
    auto forc_rh = create<ArrayD1>("forc_rh", ncells);
    auto forc_lwrad = create<ArrayD1>("forc_lwrad", ncells);
    auto forc_solai = create<ArrayD2>("forc_solai", ncells, 2);
    auto forc_solad = create<ArrayD2>("forc_solad", ncells, 2);
    auto forc_rain = create<ArrayD1>("forc_rain", ncells);
    auto forc_snow = create<ArrayD1>("forc_snow", ncells);
    auto forc_u = create<ArrayD1>("forc_u", ncells);
    auto forc_v = create<ArrayD1>("forc_v", ncells);
    auto forc_hgt = create<ArrayD1>("forc_hgt", ncells);
    auto forc_hgt_u = create<ArrayD1>("forc_hgt_u", ncells);
    auto forc_hgt_t = create<ArrayD1>("forc_hgt_t", ncells);
    auto forc_hgt_q = create<ArrayD1>("forc_hgt_q", ncells);
    auto forc_vp = create<ArrayD1>("forc_vp", ncells);
    auto forc_rho = create<ArrayD1>("forc_rho", ncells);
    auto forc_po2 = create<ArrayD1>("forc_po2", ncells);
    auto forc_pco2 = create<ArrayD1>("forc_pco2", ncells);    
    

    // phenology data
    // this is what sat phenology wants
    auto vtype = create<ArrayI1>("vtype", ncells); // pft type
    assign(vtype, 12);
    auto tlai = create<ArrayD1>("tlai", ncells);
    auto tsai = create<ArrayD1>("tsai", ncells);
    auto elai = create<ArrayD1>("elai", ncells);
    auto esai = create<ArrayD1>("esai", ncells);
    auto htop = create<ArrayD1>("htop", ncells);
    auto hbot = create<ArrayD1>("hbot", ncells);
    auto frac_veg_nosno_alb = create<ArrayI1>("frac_veg_nosno_alb", ncells);

    // snow variables
    auto snow_depth = create<ArrayD1>("snow_depth", ncells); // NEED VALUES! - probably always init at 0
    assign(snow_depth, 0.0);
    auto frac_sno = create<ArrayD1>("frac_sno", ncells);     // NEED VALUES!  \ if not glc, icemec, etc, always init these @ 0.0
    assign(frac_sno, 0.0);
    auto int_snow = create<ArrayD1>("int_snow", ncells);     // NEED VALUES!
    assign(int_snow, 0.0);


    ELM::Utils::Date current(start);
    ELM::PhenologyDataManager phen_data(dd, ncells, 17);

    std::map<std::string, ArrayD2::HostMirror> host_phen_views;
    host_phen_views["MONTHLY_LAI"] = Kokkos::create_mirror_view(phen_data.mlai_);
    host_phen_views["MONTHLY_SAI"] = Kokkos::create_mirror_view(phen_data.msai_);
    host_phen_views["MONTHLY_HEIGHT_TOP"] = Kokkos::create_mirror_view(phen_data.mhtop_);
    host_phen_views["MONTHLY_HEIGHT_BOT"] = Kokkos::create_mirror_view(phen_data.mhbot_);

    int idx = 0; // hardwire for ncells = 1
    for (int t = 0; t < ntimes; ++t) {

      // read phenology data if required
      // reader will read 3 months of data on first call
      // subsequent calls only read the newest months (when phen_data.need_data() == true)
      // and shift the index of the two remaining older months

      // copy device data to host
      // copying entire views is likely inefficient, but it's currently necessary
      // could be eliminated by shifting older month indices in parallel kernel
      // and reading new data into a mirror of a subview (or a subview of a mirror?)
      // then we would only need one copy from host view into the device view
      // instead of the two we currently have
      // will fix later - too infrequently run to cause concern  
      if (phen_data.need_data()) {
        Kokkos::deep_copy(host_phen_views["MONTHLY_LAI"], phen_data.mlai_);
        Kokkos::deep_copy(host_phen_views["MONTHLY_SAI"], phen_data.msai_);
        Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_TOP"], phen_data.mhtop_);
        Kokkos::deep_copy(host_phen_views["MONTHLY_HEIGHT_BOT"], phen_data.mhbot_);
      }

      // reads three months of data on first call
      // after first call, read new data if phen_data.need_new_data_ == true
      auto phen_updated = phen_data.read_data(host_phen_views, fname_surfdata, current, vtype); // if needed

      // copy host views to device
      // could be made more efficient, see above
      if (phen_updated) {
        Kokkos::deep_copy(phen_data.mlai_, host_phen_views["MONTHLY_LAI"]);
        Kokkos::deep_copy(phen_data.msai_, host_phen_views["MONTHLY_SAI"]);
        Kokkos::deep_copy(phen_data.mhtop_, host_phen_views["MONTHLY_HEIGHT_TOP"]);
        Kokkos::deep_copy(phen_data.mhbot_, host_phen_views["MONTHLY_HEIGHT_BOT"]);
      }

      // run parallel kernel to process phenology data
      phen_data.get_data(current, snow_depth, frac_sno, vtype, elai, esai,
        htop, hbot, tlai, tsai, frac_veg_nosno_alb);

      

      for (int i = 0; i < ncells; ++i) {
      
        std::cout << "elai: " << elai(i) << std::endl;
        std::cout << "esai: " << esai(i) << std::endl;
        std::cout << "tlai: " << tlai(i) << std::endl;
        std::cout << "tsai: " << tsai(i) << std::endl;
        std::cout << "htop: " << htop(i) << std::endl;
        std::cout << "hbot: " << hbot(i) << std::endl;
        std::cout << "frac_veg_nosno_alb: " << frac_veg_nosno_alb(i) << std::endl;
      }

      current.increment_seconds(1800);
    }

  } // inner scope

  Kokkos::finalize();
} // enclosing scope

return 0;
}
