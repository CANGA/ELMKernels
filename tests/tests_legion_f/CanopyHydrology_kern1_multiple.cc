//
// This example plays with geometric regions in Legion by figuring out one way
// to do geometric regions that are grid_cell x PFT.
//
// The first strategy is a 2D Rect IndexSpace
//

#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>

#include "legion.h"

#include "data.hh"
#include "tasks.hh"
#include "CanopyHydrology.hh"

using namespace Legion;

void top_level_task(const Task *task,
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, Runtime *runtime)
{
  std::cout << "LOG: Executing Top Level Task" << std::endl;

  const int n_pfts = 17;
  const int n_times_max = 31 * 24 * 2; // max days per month times hours per
                                       // day * half hour timestep
  const int n_grid_cells = 24;
  const int n_parts = 4;

  // -----------------------------------------------------------------------------
  // SETUP Phase
  // -----------------------------------------------------------------------------
  //
  // Create data
  //
  // grid cell x pft data for phenology
  auto phenology_fs_names = std::vector<std::string>{ "elai", "esai" };
  Data<2> phenology("phenology", ctx, runtime,
                    Point<2>(n_grid_cells, n_pfts), Point<2>(n_parts,1),
                    phenology_fs_names);
  
  // n_times_max x n_grid_cells forcing data
  auto forcing_fs_names = std::vector<std::string>{
    "forc_rain", "forc_snow", "forc_air_temp", "forc_irrig"};
  Data<2> forcing("forcing", ctx, runtime,
                  Point<2>(n_times_max, n_grid_cells), Point<2>(1,n_parts),
                  forcing_fs_names);

  // grid cell x pft water state and flux outputs
  auto flux_fs_names = std::vector<std::string>{
    "qflx_prec_intr", "qflx_irrig", "qflx_prec_grnd", "qflx_snwcp_liq",
    "qflx_snwcp_ice", "qflx_snow_grnd_patch", "qflx_rain_grnd", "h2ocan"};
  Data<2> flux("flux", ctx, runtime,
               Point<2>(n_grid_cells, n_pfts), Point<2>(n_parts,1),
               flux_fs_names);

  // -----------------------------------------------------------------------------
  // Initialization Phase
  // -----------------------------------------------------------------------------
  // launch task to read phenology
  std::cout << "LOG: Launching Init Phenology" << std::endl;
  InitPhenology().launch(ctx, runtime, phenology);

  // launch task to read forcing
  std::cout << "LOG: Launching Init Forcing" << std::endl;
  auto forc_future = InitForcing().launch(ctx, runtime, forcing);
  int n_times = forc_future.get_result<int>();
  
  // -----------------------------------------------------------------------------
  // Run Phase
  // -----------------------------------------------------------------------------
  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_multiple.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  soln_file << std::setprecision(16) << 0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << std::endl;

  // create a color space for indexed launching.  This is what a Data1D
  // color_space would look like.
  auto color_space = Rect<1>(Point<1>(0), Point<1>(n_parts-1));
  
  std::vector<Future> futures;
  for (int i=0; i!=n_times; ++i) {
    // launch interception
    CanopyHydrology_Interception().launch(ctx, runtime, color_space,
            phenology, forcing, flux, i);
    
    // launch accumulator for h2ocan
    futures.push_back(SumMinMaxReduction().launch(ctx, runtime, flux, "h2ocan"));
  }

  int i = 0;
  for (auto future : futures) {
    i++;
    //
    // write out to file
    //  
    auto sum_min_max = future.get_result<std::array<double,3>>();
    soln_file << std::setprecision(16) << i << "\t" << sum_min_max[0]
              << "\t" << sum_min_max[1]
              << "\t" << sum_min_max[2] << std::endl;
  }
}


// Main just calls top level task
int main(int argc, char **argv)
{
  Runtime::set_top_level_task_id(TaskIDs::TOP_LEVEL_TASK);

  {
    TaskVariantRegistrar registrar(TaskIDs::TOP_LEVEL_TASK, "top_level");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<top_level_task>(registrar, "top_level");
  }

  CanopyHydrology_Interception::preregister();
  InitForcing::preregister();
  InitPhenology::preregister();
  SumMinMaxReduction::preregister();

  return Runtime::start(argc, argv);
}
  



  
  
