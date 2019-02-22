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
#include "domains.hh"
#include "tasks.hh"

using namespace Legion;


void InitPhenology(const Task *task,
                   const std::vector<PhysicalRegion> &regions,
                   Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(regions.size() == 1);
  assert(task->regions[0].instance_fields.size() == 2); // LAI, SAI

  std::cout << "LOG: Executing InitPhenology task" << std::endl;
  const FieldAccessor<WRITE_DISCARD,double,2> elai(regions[0],
          task->regions[0].instance_fields[0]);
  const FieldAccessor<WRITE_DISCARD,double,2> esai(regions[0],
          task->regions[0].instance_fields[1]);

  Rect<2> my_bounds = Domain(runtime->get_index_space_domain(
      regions[0].get_logical_region().get_index_space()));
  coord_t n_grid_cells = my_bounds.hi[0] - my_bounds.lo[0] + 1;
  coord_t n_pfts = my_bounds.hi[1] - my_bounds.lo[1] + 1;
  
  assert(n_grid_cells == 24); // hard coded as two reads of 2x 12 increments
  ELM::Utils::read_phenology("../links/surfacedataWBW.nc", 12, n_pfts, 0, elai, esai);
  ELM::Utils::read_phenology("../links/surfacedataBRW.nc", 12, n_pfts, 12, elai, esai);
}

int InitForcing(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(regions.size() == 1);
  assert(task->regions[0].instance_fields.size() == 4); // rain, snow, temp, irrig

  std::cout << "LOG: Executing InitForcing task" << std::endl;
  Rect<2> my_bounds = Domain(runtime->get_index_space_domain(regions[0].get_logical_region().get_index_space()));
  coord_t n_times_max = my_bounds.hi[0] - my_bounds.lo[0] + 1;
  coord_t n_grid_cells = my_bounds.hi[1] - my_bounds.lo[1] + 1;
  
  // init rain, snow, and air temp through reader
  const FieldAccessor<WRITE_DISCARD,double,2> rain(regions[0], task->regions[0].instance_fields[0]);
  const FieldAccessor<WRITE_DISCARD,double,2> snow(regions[0], task->regions[0].instance_fields[1]);
  const FieldAccessor<WRITE_DISCARD,double,2> air_temp(regions[0], task->regions[0].instance_fields[2]);
  const FieldAccessor<WRITE_DISCARD,double,2> irrig(regions[0], task->regions[0].instance_fields[3]);
  int n_times = ELM::Utils::read_forcing("../links/forcing", n_times_max, 0, n_grid_cells,
          rain, snow, air_temp);

  // init irrig to zero
  for (size_t t=0; t!=n_times_max; ++t) {
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
  assert(regions.size() == 3);
  assert(task->regions[0].instance_fields.size() == 3);
  std::cout << "LOG: Executing Interception task" << std::endl;

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
  std::cout << "rain, snow, irrig = "
            << task->regions[0].instance_fields[0] << ","
            << task->regions[0].instance_fields[1] << ","
            << task->regions[0].instance_fields[2] << std::endl;
  const AffineAccessorRO forc_rain(regions[0], task->regions[0].instance_fields[0]);
  const AffineAccessorRO forc_snow(regions[0], task->regions[0].instance_fields[1]);
  const AffineAccessorRO forc_irrig(regions[0], task->regions[0].instance_fields[2]);

  // -- phenology
  const AffineAccessorRO elai(regions[1], task->regions[1].instance_fields[0]);
  const AffineAccessorRO esai(regions[1], task->regions[1].instance_fields[1]);

  // -- output
  const AffineAccessorRW qflx_prec_intr(regions[2], task->regions[2].instance_fields[0]);
  const AffineAccessorRW qflx_irrig(regions[2], task->regions[2].instance_fields[1]);
  const AffineAccessorRW qflx_prec_grnd(regions[2], task->regions[2].instance_fields[2]);
  const AffineAccessorRW qflx_snwcp_liq(regions[2], task->regions[2].instance_fields[3]);
 
  const AffineAccessorRW qflx_snwcp_ice(regions[2], task->regions[2].instance_fields[4]);
  const AffineAccessorRW qflx_snow_grnd_patch(regions[2], task->regions[2].instance_fields[5]); 
  const AffineAccessorRW qflx_rain_grnd(regions[2], task->regions[2].instance_fields[6]);
  const AffineAccessorRW h2ocan(regions[2], task->regions[2].instance_fields[7]);

  LogicalRegion lr = regions[2].get_logical_region();
  IndexSpaceT<2> is(lr.get_index_space());
  Rect<2> bounds = Domain(runtime->get_index_space_domain(is));

  std::cout << "LOG: With bounds: " << bounds.lo << "," << bounds.hi << std::endl;
  
  int n_irrig_steps_left = 0.;  // NOTE: still not physical quite sure what to do with this one.
  
  for (size_t g = bounds.lo[0]; g != bounds.hi[0]+1; ++g) {
    for (size_t p = bounds.lo[1]; p != bounds.hi[1]+1; ++p) {
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
  Data2D phenology(n_grid_cells, n_pfts, n_parts, "phenology", phenology_fs_names,
                   ctx, runtime);
  std::cout << "LOG: Launching Init Phenology" << std::endl;
  TaskLauncher phenology_launcher(TaskIDs::INIT_PHENOLOGY, TaskArgument(NULL, 0));
  phenology_launcher.add_region_requirement(
      RegionRequirement(phenology.logical_region, WRITE_DISCARD, EXCLUSIVE,
                        phenology.logical_region));
  for (auto name : phenology.field_names) phenology_launcher.add_field(0,phenology.field_ids[name]);
  runtime->execute_task(ctx, phenology_launcher);
  
  // n_times_max x n_grid_cells forcing data
  auto forcing_fs_ids = std::vector<std::string>{
    "forc_rain", "forc_snow", "forc_air_temp", "forc_irrig"};
  Data2D_Transposed forcing(n_grid_cells, n_times_max, n_parts, "forcing",
                            forcing_fs_ids, ctx, runtime);
  
  // grid cell x pft water state and flux outputs
  auto flux_fs_ids = std::vector<std::string>{
    "qflx_prec_intr", "qflx_irrig", "qflx_prec_grnd", "qflx_snwcp_liq",
    "qflx_snwcp_ice", "qflx_snow_grnd_patch", "qflx_rain_grnd", "h2ocan"};
  Data2D flux(n_grid_cells, n_pfts, n_parts, "flux", flux_fs_ids, ctx, runtime);

  // create a color space for indexed launching.  This is what a Data1D
  // color_space would look like.
  auto color_space = Rect<1>(Point<1>(0), Point<1>(n_parts-1));
  
  // -----------------------------------------------------------------------------
  // Initialization Phase
  // -----------------------------------------------------------------------------
  // launch task to read phenology

  

  // launch task to read forcing
  std::cout << "LOG: Launching Init Forcing" << std::endl;
  TaskLauncher forcing_launcher(TaskIDs::INIT_FORCING, TaskArgument(NULL, 0));
  forcing_launcher.add_region_requirement(
      RegionRequirement(forcing.logical_region, WRITE_DISCARD, EXCLUSIVE,
                        forcing.logical_region));
  for (auto name : forcing.field_names) forcing_launcher.add_field(0,forcing.field_ids[name]);
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

  // DEBUG HACK --ETC
  //n_times = 12;
  // END DEBUG HACK
  std::vector<Future> futures;
  for (int i=0; i!=n_times; ++i) {
    auto args = std::make_tuple(i, dtime, ltype, ctype, urbpoi, do_capsnow, dewmx, frac_veg_nosno);

    ArgumentMap arg_map;
    IndexLauncher interception_launcher(TaskIDs::CANOPY_HYDROLOGY_INTERCEPTION,
            color_space, TaskArgument(&args, sizeof(args)), arg_map);

    // -- permissions on forcing
    interception_launcher.add_region_requirement(
        RegionRequirement(forcing.logical_partition, forcing.projection_id,
                          READ_ONLY, EXCLUSIVE, forcing.logical_region));
    interception_launcher.add_field(0, forcing.field_ids["forc_rain"]);
    interception_launcher.add_field(0, forcing.field_ids["forc_snow"]);
    interception_launcher.add_field(0, forcing.field_ids["forc_irrig"]);

    // -- permissions on phenology
    interception_launcher.add_region_requirement(
        RegionRequirement(phenology.logical_partition, phenology.projection_id,
                          READ_ONLY, EXCLUSIVE, phenology.logical_region));
    interception_launcher.add_field(1, phenology.field_ids["elai"]);
    interception_launcher.add_field(1, phenology.field_ids["esai"]);

    // -- permissions on output
    interception_launcher.add_region_requirement(
        RegionRequirement(flux.logical_partition, flux.projection_id,
                          READ_WRITE, EXCLUSIVE, flux.logical_region));
    for (auto name : flux.field_names) interception_launcher.add_field(2,flux.field_ids[name]);

    // -- launch the interception
    runtime->execute_index_space(ctx, interception_launcher);

    // launch accumulator for h2ocan
    SumMinMaxReduction launcher;
    futures.push_back(launcher.launch(ctx, runtime, flux, "h2ocan"));
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

  SumMinMaxReduction::preregister();

  Runtime::preregister_projection_functor(Data2D::projection_id,
          new Data2D::LocalProjectionFunction());
  Runtime::preregister_projection_functor(Data2D_Transposed::projection_id,
          new Data2D_Transposed::LocalProjectionFunction());
  
  
  return Runtime::start(argc, argv);
}
  



  
  
