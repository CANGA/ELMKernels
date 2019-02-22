#include "legion.h"
#include "tasks.hh"

using namespace Legion;

Future
SumMinMaxReduction::launch(Context ctx, Runtime *runtime,
                           Data2D& domain, const std::string& fname)
{
  TaskLauncher accumlate_launcher(id, TaskArgument());
  accumlate_launcher.add_region_requirement(
      RegionRequirement(domain.logical_region, READ_ONLY, EXCLUSIVE,
                        domain.logical_region));
  accumlate_launcher.add_field(0, domain.field_ids[fname]);
  return runtime->execute_task(ctx, accumlate_launcher);
}

std::array<double,3>
SumMinMaxReduction::task(const Task *task,
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
SumMinMaxReduction::preregister() {
  TaskVariantRegistrar registrar(id, name.c_str());
  registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
  registrar.set_leaf();
  Runtime::preregister_task_variant<std::array<double,3>,task>(registrar, name.c_str());
}    

TaskID SumMinMaxReduction::id = TaskIDs::UTIL_SUM_MIN_MAX_REDUCTION;
std::string SumMinMaxReduction::name = "sum_min_max_reduction";


