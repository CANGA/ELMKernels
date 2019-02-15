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

#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"
#include "legion.h"

using namespace Legion;

namespace TaskIDs {
enum TaskIDs {
  TOP_LEVEL_TASK,
  INIT_PHENOLOGY,
  INIT_FORCING
};
} // namespace

namespace FieldIDs {
enum FieldIDs {
  ELAI,
  ESAI,
  FORC_AIR_TEMP,
  FORC_RAIN,
  FORC_SNOW,
  FORC_IRRIG,
  QFLX_PREC_INTR,
  QFLX_IRRIG,
  QFLX_PREC_GRND,
  QFLX_SNWCP_LIQ,
  QFLX_SNWCP_ICE,
  QFLX_SNOW_GRND_PATCH,
  QFLX_RAIN_GRND
};
} // namespace

static const int n_months = 12;
static const int n_pfts = 17;
static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                            // day * half hour timestep
static const int n_grid_cells = 24;


void InitPhenology(const Task *task,
                   const std::vector<PhysicalRegion> &regions,
                   Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
  assert(task->regions[0].privilege_fields.size() == 2); // LAI, SAI

  std::cout << "Executing InitPhenology task" << std::endl;
  const FieldAccessor<WRITE_DISCARD,double,2> elai(regions[0], FieldIDs::ELAI);
  const FieldAccessor<WRITE_DISCARD,double,2> esai(regions[0], FieldIDs::ESAI);
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", n_months, n_pfts, 0, elai, esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", n_months, n_pfts, n_months, elai, esai);
}  

int InitForcing(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
  assert(task->regions[0].privilege_fields.size() == 4); // rain, snow, temp, irrig

  std::cout << "Executing InitForcing task" << std::endl;

  // init rain, snow, and air temp through reader
  const FieldAccessor<WRITE_DISCARD,double,2> rain(regions[0], FieldIDs::FORC_RAIN);
  const FieldAccessor<WRITE_DISCARD,double,2> snow(regions[0], FieldIDs::FORC_SNOW);
  const FieldAccessor<WRITE_DISCARD,double,2> air_temp(regions[0], FieldIDs::FORC_AIR_TEMP);
  int n_times = ELM::Utils::read_forcing("../links/forcing", n_max_times, 0, n_grid_cells,
          rain, snow, air_temp);

  // init irrig to zero
  const FieldAccessor<WRITE_DISCARD,double,2> irrig(regions[0], FieldIDs::FORC_IRRIG);
  for (size_t t=0; t!=n_max_times; ++t) {
    for (size_t g=0; g!=n_grid_cells; ++g) {
      irrig[t][g] = 0.;
    }
  }

  return n_times;
}  


void top_level_task(const Task *task,
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, Runtime *runtime)
{
  std::cout << "Executing Top Level Task" << std::endl;

  // -----------------------------------------------------------------------------
  // SETUP Phase
  // -----------------------------------------------------------------------------
  //
  // Create index spaces
  //
  // create a domain and index space for PFT-state
  const Rect<2> pft_state_domain(Point<2>(0,0), Point<2>(n_grid_cells-1, n_pfts-1));
  IndexSpace pft_state_is = runtime->create_index_space(ctx, pft_state_domain); 
  printf("Created index space for PFTs: %x\n", pft_state_is.get_id());

  // create a domain and index space for forcing data
  const Rect<2> forcing_domain(Point<2>(0,0), Point<2>(n_max_times-1,n_grid_cells-1));
  IndexSpace forcing_is = runtime->create_index_space(ctx, forcing_domain); 
  printf("Created index space for forcing data: %x\n", forcing_is.get_id());

  // create a domain and index space for grid cell-state
  const Rect<1> gc_state_domain(Point<1>(0), Point<1>(n_grid_cells-1));
  IndexSpace gc_state_is = runtime->create_index_space(ctx, gc_state_domain); 
  printf("Created index space for grid-cells: %x\n", gc_state_is.get_id());
  

  //
  // Create field spaces
  //
  // phenology field space
  FieldSpace phenology_fs = runtime->create_field_space(ctx);
  auto phenology_fs_ids = std::vector<FieldIDs::FieldIDs>{ FieldIDs::ELAI, FieldIDs::ESAI };
  printf("Created field space for phenology: %x\n", phenology_fs.get_id());
  {
    FieldAllocator allocator = runtime->create_field_allocator(ctx, phenology_fs);
    for (auto id : phenology_fs_ids) allocator.allocate_field(sizeof(double), id);
  }

  // forcing field space
  FieldSpace forcing_fs = runtime->create_field_space(ctx);
  auto forcing_fs_ids = std::vector<FieldIDs::FieldIDs>{
    FieldIDs::FORC_RAIN, FieldIDs::FORC_SNOW, FieldIDs::FORC_AIR_TEMP, FieldIDs::FORC_IRRIG};
  printf("Created field space for forcing: %x\n", forcing_fs.get_id());
  {
    FieldAllocator allocator = runtime->create_field_allocator(ctx, forcing_fs);
    for (auto id : forcing_fs_ids) allocator.allocate_field(sizeof(double), id);
  }

  // grid-cell flux data field space
  FieldSpace flux_fs = runtime->create_field_space(ctx);
  auto flux_fs_ids = std::vector<FieldIDs::FieldIDs>{
    FieldIDs::QFLX_PREC_GRND, FieldIDs::QFLX_IRRIG, FieldIDs::QFLX_SNWCP_LIQ, FieldIDs::QFLX_SNWCP_ICE,
    FieldIDs::QFLX_SNOW_GRND_PATCH, FieldIDs::QFLX_RAIN_GRND};
  printf("Created field space for PFT-level fluxes: %x\n", flux_fs.get_id());
  {
    FieldAllocator allocator = runtime->create_field_allocator(ctx, flux_fs);
    for (auto id : flux_fs_ids) allocator.allocate_field(sizeof(double), id);
  }
  
  
  //
  // Physical Regions
  //  
  // Logical region is the cross product of IndexSpace and FieldSpace --
  // create logical regions, physical regions, and the block til memory is
  // available.
  LogicalRegion phenology_lr = runtime->create_logical_region(ctx, pft_state_is, phenology_fs);
  RegionRequirement phenology_req(phenology_lr, READ_WRITE, EXCLUSIVE, phenology_lr);
  for (auto id : phenology_fs_ids) phenology_req.add_field(id);
  InlineLauncher phenology_region_launcher(phenology_req);
  PhysicalRegion phenology_region = runtime->map_region(ctx, phenology_region_launcher);

  LogicalRegion forcing_lr = runtime->create_logical_region(ctx, forcing_is, forcing_fs);
  RegionRequirement forcing_req(forcing_lr, READ_WRITE, EXCLUSIVE, forcing_lr);
  for (auto id : forcing_fs_ids) forcing_req.add_field(id);
  InlineLauncher forcing_region_launcher(forcing_req);
  PhysicalRegion forcing_region = runtime->map_region(ctx, forcing_region_launcher);

  LogicalRegion flux_lr = runtime->create_logical_region(ctx, pft_state_is, flux_fs);
  RegionRequirement flux_req(flux_lr, READ_WRITE, EXCLUSIVE, flux_lr);
  for (auto id : flux_fs_ids) flux_req.add_field(id);
  InlineLauncher flux_region_launcher(flux_req);
  PhysicalRegion flux_region = runtime->map_region(ctx, flux_region_launcher);

  
  // -----------------------------------------------------------------------------
  // END SETUP
  // -----------------------------------------------------------------------------
  
  // -----------------------------------------------------------------------------
  // Initialization Phase
  // -----------------------------------------------------------------------------
  // launch task to read phenology
  std::cout << "Launching Init Phenology" << std::endl;
  phenology_region.wait_until_valid();
  TaskLauncher phenology_launcher(TaskIDs::INIT_PHENOLOGY, TaskArgument(NULL, 0));
  phenology_launcher.add_region_requirement(
      RegionRequirement(phenology_lr, WRITE_DISCARD, EXCLUSIVE, phenology_lr));
  for (auto id : phenology_fs_ids) phenology_launcher.add_field(0,id);
  runtime->execute_task(ctx, phenology_launcher);

  // launch task to read forcing
  std::cout << "Launching Init Forcing" << std::endl;
  forcing_region.wait_until_valid();
  TaskLauncher forcing_launcher(TaskIDs::INIT_FORCING, TaskArgument(NULL, 0));
  forcing_launcher.add_region_requirement(
      RegionRequirement(forcing_lr, WRITE_DISCARD, EXCLUSIVE, forcing_lr));
  for (auto id : forcing_fs_ids) forcing_launcher.add_field(0,id);
  runtime->execute_task(ctx, forcing_launcher);

  // MUST DO BEFORE FLUX IS READY
  flux_region.wait_until_valid();

  
  // get an accessor and check phenology
  const FieldAccessor<READ_ONLY,double,2> elai(phenology_region, FieldIDs::ELAI);
  std::cout << "elai = " << elai[5][7] << "," << elai[6][9] << std::endl;

  // get an accessor and check forcing
  const FieldAccessor<READ_ONLY,double,2> air_temp(forcing_region, FieldIDs::FORC_AIR_TEMP);
  std::cout << "air temp = " << air_temp[5][7] << "," << air_temp[33][14] << std::endl;

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

  {
    TaskVariantRegistrar registrar(TaskIDs::INIT_PHENOLOGY, "initialize_phenology");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<InitPhenology>(registrar, "initialize_phenology");
  }

  {
    TaskVariantRegistrar registrar(TaskIDs::INIT_FORCING, "initialize_forcing");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<int,InitForcing>(registrar, "initialize_forcing");
  }

  return Runtime::start(argc, argv);
}
  



  
  
