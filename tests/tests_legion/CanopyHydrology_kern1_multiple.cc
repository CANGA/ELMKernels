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


std::array<double,3> SumMinMaxReduction(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
  assert(task->regions[0].privilege_fields.size() == 1);
  std::cout << "LOG: Executing SumMinMax Task" << std::endl;
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

  Rect<2> my_bounds = Domain(runtime->get_index_space_domain(regions[0].get_logical_region().get_index_space()));
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
  assert(task->regions.size() == 1);
  assert(task->regions[0].privilege_fields.size() == 4); // rain, snow, temp, irrig

  std::cout << "LOG: Executing InitForcing task" << std::endl;
  Rect<2> my_bounds = Domain(runtime->get_index_space_domain(regions[0].get_logical_region().get_index_space()));
  coord_t n_times_max = my_bounds.hi[0] - my_bounds.lo[0] + 1;
  coord_t n_grid_cells = my_bounds.hi[1] - my_bounds.lo[1] + 1;
  
  // init rain, snow, and air temp through reader
  const FieldAccessor<WRITE_DISCARD,double,2> rain(regions[0], FieldIDs::FORC_RAIN);
  const FieldAccessor<WRITE_DISCARD,double,2> snow(regions[0], FieldIDs::FORC_SNOW);
  const FieldAccessor<WRITE_DISCARD,double,2> air_temp(regions[0], FieldIDs::FORC_AIR_TEMP);
  int n_times = ELM::Utils::read_forcing("../links/forcing", n_times_max, 0, n_grid_cells,
          rain, snow, air_temp);

  // init irrig to zero
  const FieldAccessor<WRITE_DISCARD,double,2> irrig(regions[0], FieldIDs::FORC_IRRIG);
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
  assert(task->regions.size() == 3);
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
  auto phenology_fs_ids = std::vector<unsigned>{ FieldIDs::ELAI, FieldIDs::ESAI };
  Data2D phenology(n_grid_cells, n_pfts, n_parts, "phenology", phenology_fs_ids,
                   ctx, runtime);

  // ntimes x grid cells forcing data
  auto forcing_fs_ids = std::vector<unsigned>{
    FieldIDs::FORC_RAIN, FieldIDs::FORC_SNOW, FieldIDs::FORC_AIR_TEMP, FieldIDs::FORC_IRRIG};
  Data2D_Transposed forcing(n_grid_cells, n_times_max, n_parts, "forcing", forcing_fs_ids,
                            ctx, runtime);
  
  // grid cell x pft water state and flux outputs
  auto flux_fs_ids = std::vector<unsigned>{
    FieldIDs::QFLX_PREC_INTR, FieldIDs::QFLX_IRRIG, FieldIDs::QFLX_PREC_GRND,
    FieldIDs::QFLX_SNWCP_LIQ, FieldIDs::QFLX_SNWCP_ICE,
    FieldIDs::QFLX_SNOW_GRND_PATCH, FieldIDs::QFLX_RAIN_GRND,
    FieldIDs::H2O_CAN};
  Data2D flux(n_grid_cells, n_pfts, n_parts, "flux", flux_fs_ids, ctx, runtime);

  // create a color space for indexed launching.  This is what a Data1D
  // color_space would look like.
  auto color_space = Rect<1>(Point<1>(0), Point<1>(n_parts-1));
  
  // -----------------------------------------------------------------------------
  // Initialization Phase
  // -----------------------------------------------------------------------------
  // launch task to read phenology
  std::cout << "LOG: Launching Init Phenology" << std::endl;
  TaskLauncher phenology_launcher(TaskIDs::INIT_PHENOLOGY, TaskArgument(NULL, 0));
  phenology_launcher.add_region_requirement(
      RegionRequirement(phenology.logical_region, WRITE_DISCARD, EXCLUSIVE,
                        phenology.logical_region));
  for (auto id : phenology_fs_ids) phenology_launcher.add_field(0,id);
  runtime->execute_task(ctx, phenology_launcher);

  // launch task to read forcing
  std::cout << "LOG: Launching Init Forcing" << std::endl;
  TaskLauncher forcing_launcher(TaskIDs::INIT_FORCING, TaskArgument(NULL, 0));
  forcing_launcher.add_region_requirement(
      RegionRequirement(forcing.logical_region, WRITE_DISCARD, EXCLUSIVE,
                        forcing.logical_region));
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
    interception_launcher.add_field(0, FieldIDs::FORC_RAIN);
    interception_launcher.add_field(0, FieldIDs::FORC_SNOW);
    interception_launcher.add_field(0, FieldIDs::FORC_IRRIG);

    // -- permissions on phenology
    interception_launcher.add_region_requirement(
        RegionRequirement(phenology.logical_partition, phenology.projection_id,
                          READ_ONLY, EXCLUSIVE, phenology.logical_region));
    interception_launcher.add_field(1, FieldIDs::ELAI);
    interception_launcher.add_field(1, FieldIDs::ESAI);

    // -- permissions on output
    interception_launcher.add_region_requirement(
        RegionRequirement(flux.logical_partition, flux.projection_id,
                          READ_WRITE, EXCLUSIVE, flux.logical_region));
    for (auto id : flux_fs_ids) interception_launcher.add_field(2, id);

    // -- launch the interception
    runtime->execute_index_space(ctx, interception_launcher);

    // launch accumulator for h2ocan
    // NOTE: need to somehow make this a reduction or something?  Shouldn't be on the full region!
    TaskLauncher accumlate_launcher(TaskIDs::UTIL_SUM_MIN_MAX_REDUCTION, TaskArgument());
    accumlate_launcher.add_region_requirement(
        RegionRequirement(flux.logical_region, READ_ONLY, EXCLUSIVE, flux.logical_region));
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

  Runtime::preregister_projection_functor(Data2D::projection_id,
          new Data2D::LocalProjectionFunction());
  Runtime::preregister_projection_functor(Data2D_Transposed::projection_id,
          new Data2D_Transposed::LocalProjectionFunction());
  
  
  return Runtime::start(argc, argv);
}
  



  
  
