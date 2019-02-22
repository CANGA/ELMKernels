//! --------------------------------------------------------------------------------
//
// Arcos -- Legion
//
// Author: Ethan Coon (coonet@ornl.gov)
// License: BSD
//
// A TaskManager is a struct for holding several static methods used
// in pushing tasks onto the Legion runtime.  All TaskManagers must supply:
//
// struct TaskManager {
//   static Legion::TaskID taskid;
//   static void preregister(Legion::TaskID new_taskid = AUTO_GENERATE_ID);
//   static Legion::Future launch(Legion::Context ctx, Legion::Runtime *runtime, ...);
//   static return_t cpu_execute_task(const Legion::Task *task,
// 		        const std::vector<Legion::PhysicalRegion> &regions,
// 		        Legion::Context ctx, Legion::Runtime *runtime);
// };
//
//
// preregister_task() should call some variant of
//	LegionRuntime::HighLevel::register_X() or
//	LegionRuntime::HighLevel::preregister_Y()
// This function is called exactly once on each task.
//
// launch() should take whatever arguments (parameters, future lists,
// etc) needed to bundle and spawn the task.  This is called each time
// a task is passed off to the runtime queue.
//
// cpu_execute_task() should be the actual task implementation, which unpacks
// the futures/args and does the work.  This interface is fixed by Legion,
// with the exception of the return type.
//
// ---------------------------------------------------------------------------------

#ifndef TASKS_HH_
#define TASKS_HH_


#include <array>
#include "legion.h"
#include "data.hh"

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

//
// Reduction on a 2D domain, calculates sum, min, and max value of a field.
//
// NOTE: need to somehow make this a reduction or something?  Shouldn't be
// global on the full region!
struct SumMinMaxReduction {
  Future launch(Context ctx, Runtime *runtime,
                Data<2>& domain, const std::string& fname);  
  static std::array<double,3> cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
};  


//
// Reads phenology data from netcdf files.
//
// If this grows, will need to make it index-launched.
struct InitPhenology {
  Future launch(Context ctx, Runtime *runtime, Data<2>& data);
  static void cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
};  


//
// Reads forcing data from netcdf files.
//
// If this grows, will need to make it index-launched.
struct InitForcing {
  Future launch(Context ctx, Runtime *runtime, Data<2>& data);
  static int cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
};  


//
// Task manager for CanopyHydrology interception
//
struct CanopyHydrology_Interception {
  FutureMap launch(Context ctx, Runtime *runtime,
                   Rect<1>& color_space,
                   Data<2>& phenology_data,
                   Data<2>& forcing_data,
                   Data<2>& flux,
                   int itime);
  static void cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
};  




#endif
