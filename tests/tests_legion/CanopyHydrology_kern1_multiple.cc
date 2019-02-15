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

#include "utils.hh"
#include "readers.hh"
#include "CanopyHydrology.hh"
#include "legion.h"

using namespace Legion;

namespace TaskIDs {
enum TaskIDs {
  TOP_LEVEL_TASK,
  INIT_PHENOLOGY,
  INIT_FORCING,
  CANOPY_HYDROLOGY_INTERCEPTION,
  UTIL_SUM_MIN_MAX_REDUCTION
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
  QFLX_RAIN_GRND,
  H2O_CAN  
};
} // namespace

static const int n_months = 12;
static const int n_pfts = 17;
static const int n_max_times = 31 * 24 * 2; // max days per month times hours per
                                            // day * half hour timestep
static const int n_grid_cells = 24;


std::array<double,3> SumMinMaxReduction(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
  assert(task->regions[0].privilege_fields.size() == 1);
  FieldID fid = *(task->regions[0].privilege_fields.begin());

  FieldAccessor<READ_ONLY,double,2,coord_t,
                Realm::AffineAccessor<double,2,coord_t> > field(regions[0], fid);
  Rect<2> rect = runtime->get_index_space_domain(ctx,
          task->regions[0].region.get_index_space());

  std::array<double,3> sum_min_max = {0., 0., 0.};
  for (PointInRectIterator<2> pir(rect); pir(); pir++) {  
    auto val = field[*pir];
    sum_min_max[0] += val;
    sum_min_max[1] = std::min(sum_min_max[1], val);
    sum_min_max[2] = std::max(sum_min_max[2], val);
  }
  return sum_min_max;
}
  


void InitPhenology(const Task *task,
                   const std::vector<PhysicalRegion> &regions,
                   Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
  assert(task->regions[0].privilege_fields.size() == 2); // LAI, SAI

  std::cout << "LOG: Executing InitPhenology task" << std::endl;
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

  std::cout << "LOG: Executing InitForcing task" << std::endl;

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


void CanopyHydrology_Interception_task(const Task *task,
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, Runtime *runtime)
{
  assert(regions.size() == 3);
  assert(task->regions.size() == 3);

  // process args / parameters
  int lcv_time;
  double dtime, dewmx;
  int ltype, ctype, frac_veg_nosno;
  bool urbpoi, do_capsnow;
  using args_t = std::tuple<int, double, int, int, bool, bool, double, int>;
  std::tie(lcv_time, dtime, ltype, ctype, urbpoi, do_capsnow, dewmx, frac_veg_nosno) =
      *((args_t*) task->args);

  // get accessors
  using AffineAccessorRO = FieldAccessor<READ_ONLY,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  using AffineAccessorRW = FieldAccessor<READ_WRITE,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  
  // -- forcing
  const AffineAccessorRO forc_rain(regions[0], FieldIDs::FORC_RAIN);
  const AffineAccessorRO forc_snow(regions[0], FieldIDs::FORC_SNOW);
  const AffineAccessorRO forc_irrig(regions[0], FieldIDs::FORC_IRRIG);

  // -- phenology
  const AffineAccessorRO elai(regions[1], FieldIDs::ELAI);
  const AffineAccessorRO esai(regions[1], FieldIDs::ESAI);

  // -- output
  const AffineAccessorRW qflx_prec_intr(regions[2], FieldIDs::QFLX_PREC_INTR);
  const AffineAccessorRW qflx_irrig(regions[2], FieldIDs::QFLX_IRRIG);
  const AffineAccessorRW qflx_prec_grnd(regions[2], FieldIDs::QFLX_PREC_GRND);
  const AffineAccessorRW qflx_snwcp_liq(regions[2], FieldIDs::QFLX_SNWCP_LIQ);
  const AffineAccessorRW qflx_snwcp_ice(regions[2], FieldIDs::QFLX_SNWCP_ICE);
  const AffineAccessorRW qflx_snow_grnd_patch(regions[2], FieldIDs::QFLX_SNOW_GRND_PATCH);
  const AffineAccessorRW qflx_rain_grnd(regions[2], FieldIDs::QFLX_RAIN_GRND);
  const AffineAccessorRW h2ocan(regions[2], FieldIDs::H2O_CAN);


  int n_irrig_steps_left = 0.;  // NOTE: still not physical quite sure what to do with this one.
  
  for (size_t g = 0; g != n_grid_cells; ++g) {
    for (size_t p = 0; p != n_pfts; ++p) {
      ELM::CanopyHydrology_Interception(dtime,
              forc_rain[lcv_time][g], forc_snow[lcv_time][g], forc_irrig[lcv_time][g],
              ltype, ctype, urbpoi, do_capsnow,
              elai[g][p], esai[g][p], dewmx, frac_veg_nosno,
              h2ocan[g][p], n_irrig_steps_left,
              qflx_prec_intr[g][p], qflx_irrig[g][p], qflx_prec_grnd[g][p],
              qflx_snwcp_liq[g][p], qflx_snwcp_ice[g][p],
              qflx_snow_grnd_patch[g][p], qflx_rain_grnd[g][p]);
    }
  }
}
  
  

void top_level_task(const Task *task,
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, Runtime *runtime)
{
  std::cout << "LOG: Executing Top Level Task" << std::endl;

  // -----------------------------------------------------------------------------
  // SETUP Phase
  // -----------------------------------------------------------------------------
  //
  // Create index spaces
  //
  // create a domain and index space for PFT-state
  const Rect<2> pft_state_domain(Point<2>(0,0), Point<2>(n_grid_cells-1, n_pfts-1));
  IndexSpace pft_state_is = runtime->create_index_space(ctx, pft_state_domain); 
  std::cout << "LOG: Created index space for PFTs: " <<  pft_state_is.get_id() << std::endl;

  // create a domain and index space for forcing data
  const Rect<2> forcing_domain(Point<2>(0,0), Point<2>(n_max_times-1,n_grid_cells-1));
  IndexSpace forcing_is = runtime->create_index_space(ctx, forcing_domain); 
  std::cout << "LOG: Created index space for forcing data: " <<  forcing_is.get_id() << std::endl;

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
    FieldIDs::QFLX_PREC_INTR, FieldIDs::QFLX_IRRIG, FieldIDs::QFLX_PREC_GRND,
    FieldIDs::QFLX_SNWCP_LIQ, FieldIDs::QFLX_SNWCP_ICE,
    FieldIDs::QFLX_SNOW_GRND_PATCH, FieldIDs::QFLX_RAIN_GRND,
    FieldIDs::H2O_CAN};
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
  std::cout << "LOG: Launching Init Phenology" << std::endl;
  //  phenology_region.wait_until_valid(); // actually likely not needed, implicit wait on accessor
  TaskLauncher phenology_launcher(TaskIDs::INIT_PHENOLOGY, TaskArgument(NULL, 0));
  phenology_launcher.add_region_requirement(
      RegionRequirement(phenology_lr, WRITE_DISCARD, EXCLUSIVE, phenology_lr));
  for (auto id : phenology_fs_ids) phenology_launcher.add_field(0,id);
  runtime->execute_task(ctx, phenology_launcher);

  // launch task to read forcing
  std::cout << "LOG: Launching Init Forcing" << std::endl;
  //  forcing_region.wait_until_valid(); // actually likely not needed, implicit wait on accessor
  TaskLauncher forcing_launcher(TaskIDs::INIT_FORCING, TaskArgument(NULL, 0));
  forcing_launcher.add_region_requirement(
      RegionRequirement(forcing_lr, WRITE_DISCARD, EXCLUSIVE, forcing_lr));
  for (auto id : forcing_fs_ids) forcing_launcher.add_field(0,id);
  auto forcing_future = runtime->execute_task(ctx, forcing_launcher);
  int n_times = forcing_future.get_result<int>();

  // launch task to call interception
  std::cout << "LOG: Launching Init Forcing" << std::endl;

  // -- fixed magic parameters as arguments
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  const double dewmx = 0.1;
  const double dtime = 1800.0;

  std::ofstream soln_file;
  soln_file.open("test_CanopyHydrology_kern1_multiple.soln");
  soln_file << "Time\t Total Canopy Water\t Min Water\t Max Water" << std::endl;
  soln_file << std::setprecision(16) << 0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << std::endl;

  std::vector<Future> futures;
  for (int i=0; i!=n_times; ++i) {
    auto args = std::make_tuple(i, dtime, ltype, ctype, urbpoi, do_capsnow, dewmx, frac_veg_nosno);
    TaskLauncher interception_launcher(TaskIDs::CANOPY_HYDROLOGY_INTERCEPTION,
            TaskArgument(&args, sizeof(args)));

    // -- permissions on forcing
    interception_launcher.add_region_requirement(
        RegionRequirement(forcing_lr, READ_ONLY, EXCLUSIVE, forcing_lr));
    interception_launcher.add_field(0, FieldIDs::FORC_RAIN);
    interception_launcher.add_field(0, FieldIDs::FORC_SNOW);
    interception_launcher.add_field(0, FieldIDs::FORC_IRRIG);

    // -- permissions on phenology
    interception_launcher.add_region_requirement(
        RegionRequirement(phenology_lr, READ_ONLY, EXCLUSIVE, phenology_lr));
    interception_launcher.add_field(1, FieldIDs::ELAI);
    interception_launcher.add_field(1, FieldIDs::ESAI);

    // -- permissions on output
    interception_launcher.add_region_requirement(
        RegionRequirement(flux_lr, READ_WRITE, EXCLUSIVE, flux_lr));
    for (auto id : flux_fs_ids) interception_launcher.add_field(2, id);

    // -- launch the interception
    runtime->execute_task(ctx, interception_launcher);

    // launch accumulator for h2ocan
    TaskLauncher accumlate_launcher(TaskIDs::UTIL_SUM_MIN_MAX_REDUCTION, TaskArgument());
    accumlate_launcher.add_region_requirement(
        RegionRequirement(flux_lr, READ_ONLY, EXCLUSIVE, flux_lr));
    accumlate_launcher.add_field(0, FieldIDs::H2O_CAN);
    futures.push_back(runtime->execute_task(ctx, accumlate_launcher));
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

  // clean up resources
  soln_file.close();
  runtime->destroy_logical_region(ctx, phenology_lr);
  runtime->destroy_logical_region(ctx, forcing_lr);
  runtime->destroy_logical_region(ctx, flux_lr);
  runtime->destroy_field_space(ctx, phenology_fs);
  runtime->destroy_field_space(ctx, forcing_fs);
  runtime->destroy_field_space(ctx, flux_fs);
  runtime->destroy_index_space(ctx, pft_state_is);
  runtime->destroy_index_space(ctx, forcing_is);
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


  {
    TaskVariantRegistrar registrar(TaskIDs::CANOPY_HYDROLOGY_INTERCEPTION, "CanopyHydrology_interception");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<CanopyHydrology_Interception_task>(registrar, "CanopyHydrology_interception");
  }

  {
    TaskVariantRegistrar registrar(TaskIDs::UTIL_SUM_MIN_MAX_REDUCTION, "sum_min_max_reduction");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<std::array<double,3>,SumMinMaxReduction>(registrar, "sum_min_max_reduction");
  }
  
  return Runtime::start(argc, argv);
}
  



  
  
