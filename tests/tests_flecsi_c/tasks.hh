/*
    @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
   /@@/////  /@@          @@////@@ @@////// /@@
   /@@       /@@  @@@@@  @@    // /@@       /@@
   /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
   /@@////   /@@/@@@@@@@/@@       ////////@@/@@
   /@@       /@@/@@//// //@@    @@       /@@/@@
   /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
   //       ///  //////   //////  ////////  //
   Copyright (c) 2016, Los Alamos National Security, LLC
   All rights reserved.
                                                                              */

#ifndef TASKS_HH_
#define TASKS_HH_


#include <array>
#include "data.hh"

#include <flecsi-tutorial/specialization/mesh/mesh.h>
#include <flecsi/data/data.h>
#include <flecsi/execution/execution.h>

using namespace flecsi;



namespace TaskIDs {
enum TaskIDs {
  TOP_LEVEL_TASK,
  INIT_PHENOLOGY,
  INIT_FORCING,
  CANOPY_HYDROLOGY_INTERCEPTION,
  CANOPY_HYDROLOGY_FRACWET,
  SUM_OVER_PFTS,
  CANOPY_HYDROLOGY_SNOWWATER,
  CANOPY_HYDROLOGY_FRACH2OS,
  UTIL_SUM_MIN_MAX_REDUCTION,
  UTIL_SUM_MIN_MAX_REDUCTION1D
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
// Reduction on a 1D domain, calculates sum, min, and max value of a field.
//
// NOTE: need to somehow make this a reduction or something?  Shouldn't be
// global on the full region!
struct SumMinMaxReduction1D {
  Future launch(Context ctx, Runtime *runtime,
                Data<1>& domain, const std::string& fname);  
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





//
// Task manager for CanopyHydrology FracWet
//
struct CanopyHydrology_FracWet {
  FutureMap launch(Context ctx, Runtime *runtime,
                   Rect<1>& color_space,
                   Data<2>& phenology_data,
                   Data<2>& flux);
  static void cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
}; 



//
// Task manager for SumOverPFTs
//
struct SumOverPFTs {
  FutureMap launch(Context ctx, Runtime *runtime,
                  Rect<1>& color_space,
                  Data<2>& flux,
                  Data<1>& surface);
  static void cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
};


//
// Task manager for CanopyHydrology SnowWater
//
struct CanopyHydrology_SnowWater {
  FutureMap launch(Context ctx, Runtime *runtime,
        Rect<1>& color_space,
        Data<2>& forcing_data,
        Data<2>& soil,
        Data<1>& surface,
        int itime);
  static void cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
}; 



//
// Task manager for CanopyHydrology FracH2OSfc
//
struct CanopyHydrology_FracH2OSfc {
  FutureMap launch(Context ctx, Runtime *runtime,
                   Rect<1>& color_space,
                   Data<2>& soil,
                   Data<1>& surface);
  static void cpu_execute_task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister(TaskID taskid=AUTO_GENERATE_ID);
  static TaskID taskid;
  static std::string name;
}; 

#endif
