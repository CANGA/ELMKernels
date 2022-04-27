
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

template <class Array_t> Array_t create(const std::string &name, int D0) { return Array_t(name, D0); }
template <class Array_t, class T> Array_t create(const std::string &name, int D0, T val) { return Array_t(name, D0, val); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1) { return Array_t(name, D0, D1); }
template <class Array_t, class T> Array_t create(const std::string &name, int D0, int D1, T val) { return Array_t(name, D0, D1, val); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1, int D2) { return Array_t(name, D0, D1, D2); }
template <class Array_t, class T> Array_t create(const std::string &name, int D0, int D1, int D2, T val) { return Array_t(name, D0, D1, D2, val); }
template <class Array_t, typename T> void assign(Array_t &arr, T val) { Kokkos::deep_copy(arr, val); }

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
    auto snow_depth = create<ArrayD1>("snow_depth", ncells); // NEED VALUES! - probably always init at 0
    auto frac_sno = create<ArrayD1>("frac_sno", ncells);     // NEED VALUES!  \ if not glc, icemec, etc, always init these @ 0.0


    ELM::Utils::Date current(start);
    ELM::PhenologyDataManager phen_data(ncells, 17);


    int idx = 0; // hardwire for ncells = 1
    for (int t = 0; t < ntimes; ++t) {
  
      phen_data.read_data(dd, fname_surfdata, current, vtype); // if needed
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
