#include <iostream>
#include <numeric>
#include <string>
#include <functional>
#include "readers.hh"
#include "CanopyHydrology_cc.hh"
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
  std::vector<std::string> output{"qflx_prec_intr", "qflx_irrig",
                                  "qflx_prec_grnd", "qflx_snwcp_liq",
                                  "qflx_snwcp_ice", "qflx_snow_grnd_patch",
                                  "qflx_rain_grnd","h2ocan"};
  for (auto fname : output)
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





//
// CanopyHydrology FracWet task manager
//
// =============================================================================

FutureMap
CanopyHydrology_FracWet::launch(Context ctx, Runtime *runtime,
        Rect<1>& color_space,
        Data<2>& phenology,
        Data<2>& flux)
{
  const int frac_veg_nosno = 1;
  const double dewmx = 0.1;
  double fwet = 0., fdry = 0.;
  auto args = std::make_tuple(dewmx, frac_veg_nosno, fwet, fdry);
  ArgumentMap arg_map;
  IndexLauncher interception_launcher(taskid,
          color_space, TaskArgument(&args, sizeof(args)), arg_map);

  // -- permissions on phenology
  interception_launcher.add_region_requirement(
      RegionRequirement(phenology.logical_partition, phenology.projection_id,
                        READ_ONLY, EXCLUSIVE, phenology.logical_region));
  interception_launcher.add_field(0, phenology.field_ids["elai"]);
  interception_launcher.add_field(0, phenology.field_ids["esai"]);
  
  // -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(flux.logical_partition, flux.projection_id,
                        READ_WRITE, EXCLUSIVE, flux.logical_region));
  std::vector<std::string> output{"h2ocan"};
  for (auto fname : output)
    interception_launcher.add_field(2,flux.field_ids[fname]);

  // -- launch the interception
  return runtime->execute_index_space(ctx, interception_launcher);
}

void
CanopyHydrology_FracWet::cpu_execute_task(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 2);
  assert(regions.size() == 2);
  assert(task->regions[0].instance_fields.size() == 2);
  std::cout << "LOG: Executing FracWet task" << std::endl;

  // process args / parameters
  double dewmx, fwet, fdry;
  int frac_veg_nosno;
  using args_t = std::tuple<double, int, double, double>;
  std::tie(dewmx, frac_veg_nosno, fwet, fdry) =
      *((args_t*) task->args);

  // get accessors
  using AffineAccessorRO = FieldAccessor<READ_ONLY,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  using AffineAccessorRW = FieldAccessor<READ_WRITE,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  
  // -- phenology
  std::cout << "elai,esai = "
            << task->regions[0].instance_fields[0] << ","
            << task->regions[0].instance_fields[1] <<  std::endl;
  
  const AffineAccessorRO elai(regions[0], task->regions[0].instance_fields[0]);
  const AffineAccessorRO esai(regions[1], task->regions[1].instance_fields[1]);

  // -- output
  const AffineAccessorRW h2ocan(regions[1], task->regions[1].instance_fields[0]);

  LogicalRegion lr = regions[2].get_logical_region();
  IndexSpaceT<2> is(lr.get_index_space());
  Rect<2> bounds = Domain(runtime->get_index_space_domain(is));

  std::cout << "LOG: With bounds: " << bounds.lo << "," << bounds.hi << std::endl;
    
  for (size_t g = bounds.lo[0]; g != bounds.hi[0]+1; ++g) {
    for (size_t p = bounds.lo[1]; p != bounds.hi[1]+1; ++p) {
      ELM::CanopyHydrology_FracWet(frac_veg_nosno,
              h2ocan[g][p], elai[g][p], esai[g][p], dewmx, fwet, fdry);
    }
  }
}


void
CanopyHydrology_FracWet::preregister(TaskID new_taskid)
{
              
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<cpu_execute_task>(registrar, name.c_str());
}    

TaskID CanopyHydrology_FracWet::taskid = TaskIDs::CANOPY_HYDROLOGY_FRACWET;
std::string CanopyHydrology_FracWet::name = "canopy_hydrology_fracwet";





//
// Sum Over PFTS task manager
//
// =============================================================================

FutureMap
SumOverPFTs::launch(Context ctx, Runtime *runtime,
        Rect<1>& color_space,
        Data<2>& flux,
        Data<1>& surface)
{
  auto args = std::make_tuple(NULL);
  ArgumentMap arg_map;
  IndexLauncher interception_launcher(taskid,
          color_space, TaskArgument(&args, sizeof(args)), arg_map);

   
  // -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(flux.logical_partition, flux.projection_id,
                        READ_WRITE, EXCLUSIVE, flux.logical_region));
  std::vector<std::string> output{"qflx_snow_grnd_patch"};
  for (auto fname : output)
    interception_launcher.add_field(2,flux.field_ids[fname]);

  // -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(surface.logical_partition, surface.projection_id,
                        READ_WRITE, EXCLUSIVE, surface.logical_region));
  std::vector<std::string> output1{"qflx_snow_grnd_col"};
  for (auto fname1 : output1)
    interception_launcher.add_field(1,surface.field_ids[fname1]);

  // -- launch the interception
  return runtime->execute_index_space(ctx, interception_launcher);
}

void
SumOverPFTs::cpu_execute_task(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  
  assert(regions.size() == 2);
  assert(regions.size() == 2);
  assert(task->regions[0].instance_fields.size() == 2);
  std::cout << "LOG: Executing SumOverPFTs task" << std::endl;

 // get accessors
  using AffineAccessorRO = FieldAccessor<READ_ONLY,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  using AffineAccessorWO = FieldAccessor<WRITE_DISCARD,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  // -- output
  const AffineAccessorWO qflx_snow_grnd_patch(regions[0], task->regions[0].instance_fields[0]);
  
  // -- output
  const AffineAccessorRO qflx_snow_grnd_col(regions[1], task->regions[1].instance_fields[0]);
  
  LogicalRegion lr = regions[0].get_logical_region();
  IndexSpaceT<2> is(lr.get_index_space());
  Rect<2> bounds = Domain(runtime->get_index_space_domain(is));

  std::cout << "LOG: With bounds: " << bounds.lo << "," << bounds.hi << std::endl;
  

   //int n = sizeof(qflx_snow_grnd_col) / sizeof(qflx_snow_grnd_patch[0]); 
   
  for (size_t g = bounds.lo[0]; g != bounds.hi[0]+1; ++g) {
  double sum = 0 ;    
  for (size_t p = bounds.lo[1]; p != bounds.hi[1]+1; ++p) {
      sum += qflx_snow_grnd_patch[g][p];
     }
     qflx_snow_grnd_col[g] = sum ;
  }
  
}


void
SumOverPFTs::preregister(TaskID new_taskid)
{
             
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<cpu_execute_task>(registrar, name.c_str());
}    

TaskID SumOverPFTs::taskid = TaskIDs::SUM_OVER_PFTS;
std::string SumOverPFTs::name = "sum_over_pfts";





//
// CanopyHydrology SnowWater task manager
//
// =============================================================================

FutureMap
CanopyHydrology_SnowWater::launch(Context ctx, Runtime *runtime,
        Rect<1>& color_space,
        Data<2>& forcing,
        Data<2>& soil,
        Data<1>& surface,
        int itime)
{
  // launch task to call interception
  // -- fixed magic parameters as arguments
  const int ctype = 1;
  const int ltype = 1;
  const bool urbpoi = false;
  const bool do_capsnow = false;
  const int oldfflag = 0;
  const double dtime = 1800.0;
  const double qflx_snow_melt = 0.;
  double qflx_floodg = 0.0;
  const double n_melt = 0.7;

  auto args = std::make_tuple(itime,dtime, qflx_floodg,ltype, ctype, urbpoi, do_capsnow, oldfflag, qflx_snow_melt, n_melt); 
  ArgumentMap arg_map;
  IndexLauncher interception_launcher(taskid,
          color_space, TaskArgument(&args, sizeof(args)), arg_map);

    // -- permissions on forcing
  interception_launcher.add_region_requirement(
      RegionRequirement(forcing.logical_partition, forcing.projection_id,
                        READ_ONLY, EXCLUSIVE, forcing.logical_region));
  interception_launcher.add_field(0, forcing.field_ids["forc_air_temp"]);
  
  // -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(soil.logical_partition, soil.projection_id,
                        READ_WRITE, EXCLUSIVE, soil.logical_region));
  std::vector<std::string> output{"swe_old", "h2osoi_liq", "h2osoi_ice", "t_soisno", "frac_iceold",
                                  "dz", "z", "zi"};
  for (auto fname : output)
    interception_launcher.add_field(2,soil.field_ids[fname]);

  // -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(surface.logical_partition, surface.projection_id,
                        READ_WRITE, EXCLUSIVE, surface.logical_region));
  std::vector<std::string> output1{"t_grnd", "h2osno", "snow_depth", "integrated_snow",
     "qflx_snow_grnd_col", "qflx_snow_h2osfc",
     "qflx_floodc", "frac_snow_eff", "frac_sno", "frac_h2osfc", "snow_level"};
  for (auto fname1 : output1)
    interception_launcher.add_field(1,surface.field_ids[fname1]);

  // -- launch the interception
  return runtime->execute_index_space(ctx, interception_launcher);
}

void
CanopyHydrology_SnowWater::cpu_execute_task(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 3);
  assert(regions.size() == 3);
  assert(task->regions[0].instance_fields.size() == 3);
  std::cout << "LOG: Executing SnowWater task" << std::endl;
  int lcv_time;
  double dtime, qflx_snow_melt , qflx_floodg , n_melt ;
  int ltype, ctype, oldfflag;
  bool urbpoi, do_capsnow;
  using args_t = std::tuple<int,double,double,int,int,bool,bool,int,double,double>;
  std::tie(lcv_time,dtime, qflx_floodg, ltype, ctype, urbpoi, do_capsnow, oldfflag, qflx_snow_melt, n_melt) =
      *((args_t*) task->args);

  // get accessors
  using AffineAccessorRO = FieldAccessor<READ_ONLY,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  using AffineAccessorRW = FieldAccessor<READ_WRITE,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  using AffineAccessorRO1 = FieldAccessor<READ_ONLY,double,1,coord_t,
                                         Realm::AffineAccessor<double,1,coord_t> >;
  using AffineAccessorRW1 = FieldAccessor<READ_WRITE,double,1,coord_t,
                                         Realm::AffineAccessor<double,1,coord_t> >;
  using AffineAccessorRW1_int = FieldAccessor<READ_WRITE,int,1,coord_t,
                                         Realm::AffineAccessor<int,1,coord_t> >;

  
  // -- forcing
  std::cout << "air_temp = "
            << task->regions[0].instance_fields[0] << std::endl;
  const AffineAccessorRO forc_air_temp(regions[0], task->regions[0].instance_fields[0]);
  
  // -- output
  const AffineAccessorRW swe_old(regions[1], task->regions[1].instance_fields[0]);
  const AffineAccessorRW h2osoi_liq(regions[1], task->regions[1].instance_fields[1]);
  const AffineAccessorRW h2osoi_ice(regions[1], task->regions[1].instance_fields[2]);
  const AffineAccessorRW t_soisno(regions[1], task->regions[1].instance_fields[3]);
  const AffineAccessorRW frac_iceold(regions[1], task->regions[1].instance_fields[4]);
  //const AffineAccessorRW snl(regions[1], task->regions[1].instance_fields[5]); // surface memeber 1D
  const AffineAccessorRW dz(regions[1], task->regions[1].instance_fields[5]);
  const AffineAccessorRW z(regions[1], task->regions[1].instance_fields[6]);
  const AffineAccessorRW zi(regions[1], task->regions[1].instance_fields[7]);


  // -- output
  const AffineAccessorRW1 t_grnd(regions[2], task->regions[2].instance_fields[0]);
  const AffineAccessorRW1 h2osno(regions[2], task->regions[2].instance_fields[1]);
  const AffineAccessorRW1 snow_depth(regions[2], task->regions[2].instance_fields[2]);
  const AffineAccessorRW1 integrated_snow(regions[2], task->regions[2].instance_fields[3]);
  const AffineAccessorRW1 qflx_snow_grnd_col(regions[2], task->regions[2].instance_fields[4]);
  const AffineAccessorRW1 qflx_snow_h2osfc(regions[2], task->regions[2].instance_fields[5]); 
  const AffineAccessorRW1 qflx_floodc(regions[2], task->regions[2].instance_fields[6]);
  const AffineAccessorRW1 frac_snow_eff(regions[2], task->regions[2].instance_fields[7]);
  const AffineAccessorRW1 frac_sno(regions[2], task->regions[2].instance_fields[8]);
  const AffineAccessorRW1 frac_h2osfc(regions[2], task->regions[2].instance_fields[9]);
  const AffineAccessorRW1_int snow_level(regions[2], task->regions[2].instance_fields[10]);

  LogicalRegion lr = regions[2].get_logical_region();
  IndexSpaceT<2> is(lr.get_index_space());
  Rect<2> bounds = Domain(runtime->get_index_space_domain(is));

  std::cout << "LOG: With bounds: " << bounds.lo << "," << bounds.hi << std::endl;
  
  int newnode = 0.;  
  
  for (size_t g = bounds.lo[0]; g != bounds.hi[0]+1; ++g) {
    
      ELM::CanopyHydrology_SnowWater(dtime, qflx_floodg,
              ltype, ctype, urbpoi, do_capsnow, oldfflag,
              forc_air_temp[lcv_time][g], t_grnd[g],
              qflx_snow_grnd_col[g], qflx_snow_melt, n_melt, frac_h2osfc[g],
              snow_depth[g], h2osno[g], integrated_snow[g], swe_old[g],
              h2osoi_liq[g], h2osoi_ice[g], t_soisno[g], frac_iceold[g],
              snow_level[g], dz[g], z[g], zi[g], newnode,
              qflx_floodc[g], qflx_snow_h2osfc[g], frac_snow_eff[g], frac_sno[g]);
    
  }
}


void
CanopyHydrology_SnowWater::preregister(TaskID new_taskid)
{
  // taskid = (taskid == AUTO_GENERATE_ID ?
  //           Legion::Runtime::generate_static_task_id() :
  //             new_taskid);
            
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<cpu_execute_task>(registrar, name.c_str());
}    

TaskID CanopyHydrology_SnowWater::taskid = TaskIDs::CANOPY_HYDROLOGY_SNOWWATER;
std::string CanopyHydrology_SnowWater::name = "canopy_hydrology_snowwater";






//
// CanopyHydrology FracH2OSfc task manager
//
// =============================================================================

FutureMap
CanopyHydrology_FracH2OSfc::launch(Context ctx, Runtime *runtime,
        Rect<1>& color_space,
        Data<2>& soil,
        Data<1>& surface)
{
  
  const int ltype = 1;
  const double dtime = 1800.0;
  const double micro_sigma = 0.1;
  const double min_h2osfc = 1.0e-8;

  auto args = std::make_tuple(dtime, min_h2osfc, micro_sigma ,ltype);
  ArgumentMap arg_map;
  IndexLauncher interception_launcher(taskid,
          color_space, TaskArgument(&args, sizeof(args)), arg_map);

// -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(soil.logical_partition, soil.projection_id,
                        READ_WRITE, EXCLUSIVE, soil.logical_region));
  std::vector<std::string> output{"h2osoi_liq"};
  for (auto fname : output)
    interception_launcher.add_field(2,soil.field_ids[fname]);
  
  // -- permissions on output
  interception_launcher.add_region_requirement(
      RegionRequirement(surface.logical_partition, surface.projection_id,
                        READ_WRITE, EXCLUSIVE, surface.logical_region));
  std::vector<std::string> output1{"h2osno", "h2osfc", "frac_h2osfc", "qflx_h2osfc2topsoi", "frac_snow_eff", "frac_sno"}; //h2osnow rw thing
  for (auto fname1 : output1)
    interception_launcher.add_field(1,surface.field_ids[fname1]);

  // -- launch the interception
  return runtime->execute_index_space(ctx, interception_launcher);
}

void
CanopyHydrology_FracH2OSfc::cpu_execute_task(const Task *task,
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  assert(regions.size() == 2);
  assert(regions.size() == 2);
  assert(task->regions[0].instance_fields.size() == 2);
  std::cout << "LOG: Executing FracH2OSfc task" << std::endl;

  double dtime, micro_sigma , min_h2osfc ;
  int ltype;
  using args_t = std::tuple<double,double,double,int>;
  std::tie(dtime, min_h2osfc, micro_sigma ,ltype) =
      *((args_t*) task->args);

  // get accessors
  using AffineAccessorRO = FieldAccessor<READ_ONLY,double,2,coord_t,
                                         Realm::AffineAccessor<double,2,coord_t> >;
  using AffineAccessorRO1 = FieldAccessor<READ_ONLY,double,1,coord_t,
                                         Realm::AffineAccessor<double,1,coord_t> >;                                         
  using AffineAccessorRW = FieldAccessor<READ_WRITE,double,1,coord_t,
                                         Realm::AffineAccessor<double,1,coord_t> >;
  using AffineAccessorWO = FieldAccessor<WRITE_DISCARD,double,1,coord_t,
                                         Realm::AffineAccessor<double,1,coord_t> >;
  // -- output
  const AffineAccessorRO h2osoi_liq(regions[0], task->regions[0].instance_fields[0]);
  
  // -- output
  const AffineAccessorRW h2osno(regions[1], task->regions[1].instance_fields[0]);
  const AffineAccessorRW h2osfc(regions[1], task->regions[1].instance_fields[1]);
  const AffineAccessorWO frac_h2osfc(regions[1], task->regions[1].instance_fields[2]);
  const AffineAccessorWO qflx_h2osfc2topsoi(regions[1], task->regions[1].instance_fields[3]);
  const AffineAccessorRW frac_snow_eff(regions[1], task->regions[1].instance_fields[4]);
  const AffineAccessorRW frac_sno(regions[1], task->regions[1].instance_fields[5]);
  
  LogicalRegion lr = regions[1].get_logical_region();
  IndexSpaceT<1> is(lr.get_index_space());
  Rect<1> bounds = Domain(runtime->get_index_space_domain(is));

  std::cout << "LOG: With bounds: " << bounds.lo << "," << bounds.hi << std::endl;
  
   
  for (size_t g = bounds.lo[0]; g != bounds.hi[0]+1; ++g) {
    //for (size_t p = bounds.lo[1]; p != bounds.hi[1]+1; ++p) {
         ELM::CanopyHydrology_FracH2OSfc(dtime, min_h2osfc, ltype, micro_sigma,
              h2osno[g], h2osfc[g], h2osoi_liq[g][0], frac_sno[g], frac_snow_eff[g],qflx_h2osfc2topsoi[g], frac_h2osfc[g] );
   //} 
  }
}


void
CanopyHydrology_FracH2OSfc::preregister(TaskID new_taskid)
{
             
  TaskVariantRegistrar registrar(taskid, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<cpu_execute_task>(registrar, name.c_str());
}    

TaskID CanopyHydrology_FracH2OSfc::taskid = TaskIDs::CANOPY_HYDROLOGY_FRACH2OS;
std::string CanopyHydrology_FracH2OSfc::name = "canopy_hydrology_frach2os";