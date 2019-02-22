#include "readers.hh"
#include "CanopyHydrology.hh"
#include "legion.h"
#include "tasks.hh"

using namespace Legion;

//
// SumMinMaxReduction task
//
// =============================================================================

Future
SumMinMaxReduction::launch(Context ctx, Runtime *runtime,
                           Data<2>& domain, const std::string& fname)
{
  TaskLauncher accumlate_launcher(taskid, TaskArgument());
  accumlate_launcher.add_region_requirement(
      RegionRequirement(domain.logical_region, READ_ONLY, EXCLUSIVE,
                        domain.logical_region));
  accumlate_launcher.add_field(0, domain.field_ids[fname]);
  return runtime->execute_task(ctx, accumlate_launcher);
}

std::array<double,3>
SumMinMaxReduction::cpu_execute_task(const Task *task,
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

void
SumMinMaxReduction::preregister(TaskID new_taskid) {
  // taskid = (taskid == AUTO_GENERATE_ID ?
  //           Legion::Runtime::generate_static_task_id() :
  //             new_taskid);
            
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<std::array<double,3>,cpu_execute_task>(registrar, name.c_str());
}    

TaskID SumMinMaxReduction::taskid = TaskIDs::UTIL_SUM_MIN_MAX_REDUCTION;
std::string SumMinMaxReduction::name = "sum_min_max_reduction";



//
// InitPhenology task
//
// =============================================================================
Future
InitPhenology::launch(Context ctx, Runtime *runtime, Data<2>& data)
{
  TaskLauncher phenology_launcher(taskid, TaskArgument(NULL, 0));
  phenology_launcher.add_region_requirement(
      RegionRequirement(data.logical_region, WRITE_DISCARD, EXCLUSIVE,
                        data.logical_region));

  phenology_launcher.add_field(0, data.field_ids["elai"]);
  phenology_launcher.add_field(0, data.field_ids["esai"]);
  return runtime->execute_task(ctx, phenology_launcher);
}

void
InitPhenology::cpu_execute_task(const Task *task,
                   const std::vector<PhysicalRegion> &regions,
                   Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
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

void
InitPhenology::preregister(TaskID new_taskid)
{
  // taskid = (taskid == AUTO_GENERATE_ID ?
  //           Legion::Runtime::generate_static_task_id() :
  //             new_taskid);
            
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<cpu_execute_task>(registrar, name.c_str());
}    

TaskID InitPhenology::taskid = TaskIDs::INIT_PHENOLOGY;
std::string InitPhenology::name = "init_phenology";



//
// InitForcing task
//
// =============================================================================
Future
InitForcing::launch(Context ctx, Runtime *runtime, Data<2>& data)
{

  std::cout << "LOG: Launching Init Forcing" << std::endl;
  
  TaskLauncher forcing_launcher(taskid, TaskArgument(NULL, 0));
  forcing_launcher.add_region_requirement(
      RegionRequirement(data.logical_region, WRITE_DISCARD, EXCLUSIVE,
                        data.logical_region));

  forcing_launcher.add_field(0, data.field_ids["forc_rain"]);
  forcing_launcher.add_field(0, data.field_ids["forc_snow"]);
  forcing_launcher.add_field(0, data.field_ids["forc_air_temp"]);
  forcing_launcher.add_field(0, data.field_ids["forc_irrig"]);
  return runtime->execute_task(ctx, forcing_launcher);
}

int
InitForcing::cpu_execute_task(const Task *task,
                   const std::vector<PhysicalRegion> &regions,
                   Context ctx, Runtime *runtime)
{
  assert(regions.size() == 1);
  assert(task->regions.size() == 1);
  assert(task->regions[0].instance_fields.size() == 4);

  std::cout << "LOG: Executing InitForcing task" << std::endl;
  Rect<2> my_bounds = Domain(runtime->get_index_space_domain(
      regions[0].get_logical_region().get_index_space()));
  coord_t n_times_max = my_bounds.hi[0] - my_bounds.lo[0] + 1;
  coord_t n_grid_cells = my_bounds.hi[1] - my_bounds.lo[1] + 1;
  
  // init rain, snow, and air temp through reader
  const FieldAccessor<WRITE_DISCARD,double,2> rain(regions[0],
          task->regions[0].instance_fields[0]);
  const FieldAccessor<WRITE_DISCARD,double,2> snow(regions[0],
          task->regions[0].instance_fields[1]);
  const FieldAccessor<WRITE_DISCARD,double,2> air_temp(regions[0],
          task->regions[0].instance_fields[2]);
  const FieldAccessor<WRITE_DISCARD,double,2> irrig(regions[0],
          task->regions[0].instance_fields[3]);
  int n_times = ELM::Utils::read_forcing("../links/forcing",
          n_times_max, 0, n_grid_cells,
          rain, snow, air_temp);

  // init irrig to zero
  for (size_t t=0; t!=n_times_max; ++t) {
    for (size_t g=0; g!=n_grid_cells; ++g) {
      irrig[t][g] = 0.;
    }
  }
  return n_times;
}  

void
InitForcing::preregister(TaskID new_taskid)
{
  // taskid = (taskid == AUTO_GENERATE_ID ?
  //           Legion::Runtime::generate_static_task_id() :
  //             new_taskid);
            
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<int,cpu_execute_task>(registrar, name.c_str());
}    

TaskID InitForcing::taskid = TaskIDs::INIT_FORCING;
std::string InitForcing::name = "init_forcing";


//
// CanopyHydrology Interception task manager
//
// =============================================================================

FutureMap
CanopyHydrology_Interception::launch(Context ctx, Runtime *runtime,
        Rect<1>& color_space,
        Data<2>& phenology,
        Data<2>& forcing,
        Data<2>& flux,
        int itime)
{
  // launch task to call interception
  // -- fixed magic parameters as arguments
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int frac_veg_nosno = 1;
  const double dewmx = 0.1;
  const double dtime = 1800.0;

  auto args = std::make_tuple(itime, dtime, ltype, ctype, urbpoi,
          do_capsnow, dewmx, frac_veg_nosno);
  ArgumentMap arg_map;
  IndexLauncher interception_launcher(taskid,
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
  for (auto fname : flux.field_names)
    interception_launcher.add_field(2,flux.field_ids[fname]);

  // -- launch the interception
  return runtime->execute_index_space(ctx, interception_launcher);
}

void
CanopyHydrology_Interception::cpu_execute_task(const Task *task,
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


void
CanopyHydrology_Interception::preregister(TaskID new_taskid)
{
  // taskid = (taskid == AUTO_GENERATE_ID ?
  //           Legion::Runtime::generate_static_task_id() :
  //             new_taskid);
            
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<cpu_execute_task>(registrar, name.c_str());
}    

TaskID CanopyHydrology_Interception::taskid = TaskIDs::CANOPY_HYDROLOGY_INTERCEPTION;
std::string CanopyHydrology_Interception::name = "canopy_hydrology_interception";
