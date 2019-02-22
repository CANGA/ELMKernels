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
  const int n_levels_soil_col = 5;  // NOTE: this will change once we get soil,
                                    // currently this is just n_lev_snow
  

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
  //
  // NOTE: Combine this with phenology I think?  Or look into ELM and how they
  // split these things.
  auto flux_fs_names = std::vector<std::string>{
    "qflx_prec_intr", "qflx_irrig", "qflx_prec_grnd", "qflx_snwcp_liq",
    "qflx_snwcp_ice", "qflx_snow_grnd_patch", "qflx_rain_grnd",
    "h2ocan", "fwet", "fdry"};
  Data<2> flux("flux", ctx, runtime,
                    Point<2>(n_grid_cells, n_pfts), Point<2>(n_parts,1),
                    flux_fs_names);

  // grid cell x soil/snow column
  auto soil_fs_names = std::vector<std::string>{
    "swe_old", "h2osoi_liq", "h2osoi_ice", "t_soisno", "frac_iceold"};
  Data<2> soil("soil", ctx, runtime,
               Point<2>(n_grid_cells, n_levels_soil_col), Point<2>(n_parts,1),
               soil_fs_names);

  // grid cell data
  // NOTE this has an int in it, not just doubles!
  auto grid_cell_fs_names_double = std::vector<std::string>{
    "t_grnd", "h2osno", "snow_depth", "integrated_snow",
    "h2osfc", "frac_h2osfc", "qflx_snow_grnd_col", "qflx_snow_h2osfc",
    "qflx_h2osfc2topsoi", "qflx_floodc", "frac_snow_eff", "frac_sno"};
  Data<1> surface("surface", ctx, runtime,
                  Point<1>(n_grid_cells), Point<1>(n_parts));
  for (auto fname : grid_cell_fs_names_double)
    surface.addField<double>(fname);
  surface.addField<int>("snow_level");
  surface.finalize();
  
  
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

  // create a color space for indexed launching.  Can just launch of the
  // existing 1D data structure's color_space
  auto color_space = surface.color_domain;
  
  std::vector<Future> futures;
  for (int i=0; i!=n_times; ++i) {
    // launch interception

    // NOTE: This is where it would be interesting to have launches over PFTs
    // as well as over grid cells.  This would allow to explore whether it
    // makes sense to keep all PFTs of the same type together or the other
    // direction.  This would require a second color_space, and a Data
    // structure that allowed multiple partitionings of the same
    // logical_region.  We would have one partitioning that included
    // partitioning over PFTs and grid cells, and one that only partitioned
    // over grid cells (for, e.g. SumOverPFTs launch below).  For now,
    // however, we'll just launch this decomposed over the same grid-cell only
    // partitioning as the other tasks.
    CanopyHydrology_Interception().launch(ctx, runtime, color_space,
            phenology, forcing, flux, i);

    // NOTE: WRITE ME!
    CanopyHydrology_FracWet().launch(ctx, runtime, color_space,
            phenology, flux);

    // launch integration kernel/task to, for each grid cell, sum over PFTs.
    // NOTE: WRITE ME!  This task should be generic!
    SumOverPFTs().launch(ctx, runtime, color_space,
                         flux, "qflx_snow_grnd_patch",
                         surface, "qflx_snow_grnd_col");

    // launch water balance kernel
    // NOTE: WRITE ME!
    CanopyHydrology_SnowWater().launch(ctx, runtime, color_space,
            forcing, surface, i);

    // launch fraction of water to surface
    // NOTE: WRITE ME!
    CanopyHydrology_FracH2OSfc().launch(ctx, runtime, color_space, surface);

    // NOTE: Figure out how to evaluate the success of this test!  launch
    // accumulators?  Print something to file?  Can we make
    // SumMinMaxReduction() both an actual reduction and dimension
    // independent?
    futures.push_back(SumMinMaxReduction().launch(ctx, runtime, flux, "h2ocan"));
    futures.push_back(SumMinMaxReduction().launch(ctx, runtime, surface, "frac_h2osfc"));
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

  InitForcing::preregister();
  InitPhenology::preregister();
  SumMinMaxReduction::preregister();
  SumOverPFTs::preregister();
  CanopyHydrology_Interception::preregister();
  CanopyHydrology_FracWet::preregister();
  CanopyHydrology_SnowWater::preregister();
  CanopyHydrology_FracH2OSfc::preregister();

  return Runtime::start(argc, argv);
}
  



  
  
