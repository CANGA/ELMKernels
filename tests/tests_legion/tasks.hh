#ifndef TASKS_HH_
#define TASKS_HH_


#include <array>
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

//
// Reduction on a 2D domain, calculates sum, min, and max value of a field.
//
// NOTE: need to somehow make this a reduction or something?  Shouldn't be
// global on the full region!
struct SumMinMaxReduction {
  Future launch(Context ctx, Runtime *runtime,
                Data2D& domain, const std::string& fname);  
  static std::array<double,3> task(const Task *task,
          const std::vector<PhysicalRegion> &regions,
          Context ctx, Runtime *runtime);
  static void preregister();
  static TaskID id;
  static std::string name;
};  

#endif
