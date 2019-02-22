#ifndef DOMAINS_HH_
#define DOMAINS_HH_

#include "legion.h"
using namespace Legion;

namespace ProjectionIDs {
enum ProjectionIDs {
  IDENTITY,
  TO_1D_FROM_2D,
  TO_1D_FROM_2D_TRANSPOSED  
};
} // namespace



struct Data1D {
  Data1D(coord_t n_grid_cells_, coord_t n_parts_,
         const std::vector<std::string>& field_names_,
         const std::string& name_,
         Context ctx_, Runtime* runtime_) :
      name(name_),
      field_names(field_names_),
      n_grid_cells(n_grid_cells_),
      domain(Point<1>(0), Point<1>(n_grid_cells_-1)),
      n_parts(n_parts_),
      color_domain(Point<1>(0), Point<1>(n_parts_-1)),
      ctx(ctx_),
      runtime(runtime_)
  {
    index_space = runtime->create_index_space<1, coord_t>(ctx, domain);

    auto block = Point<1>(n_grid_cells / n_parts);
    index_partition = runtime->create_partition_by_blockify<1,coord_t>(ctx, index_space, block);
    runtime->attach_name(index_partition, (std::string("index_partition_")+name).c_str());

    field_space = runtime->create_field_space(ctx);
    {
      FieldAllocator allocator = runtime->create_field_allocator(ctx, field_space);
      for (auto fname : field_names)
        field_ids[fname] = allocator.allocate_field(sizeof(double), AUTO_GENERATE_ID);
    }
    runtime->attach_name(field_space, (std::string("field_space_")+name).c_str());

    logical_region = runtime->create_logical_region(ctx, index_space, field_space);
    runtime->attach_name(logical_region, (std::string("logical_region_")+name).c_str());
    logical_partition = runtime->get_logical_partition(ctx, logical_region, index_partition);
    runtime->attach_name(logical_partition, (std::string("logical_partition_")+name).c_str());
  }

  ~Data1D() {
    runtime->destroy_logical_partition(ctx, logical_partition);
    runtime->destroy_logical_region(ctx, logical_region);
    runtime->destroy_field_space(ctx, field_space);
    runtime->destroy_index_partition(ctx, index_partition);
    runtime->destroy_index_space(ctx, index_space);
  }

  // name the container
  std::string name;
  std::vector<std::string> field_names;

  // index space and partition
  coord_t n_grid_cells;
  const Rect<1> domain;
  IndexSpaceT<1,coord_t> index_space;

  coord_t n_parts;
  const Rect<1> color_domain;
  IndexPartitionT<1, coord_t> index_partition;

  // field space
  std::map<std::string, unsigned> field_ids;
  FieldSpace field_space;

  // cross product of index and field
  LogicalRegion logical_region;
  LogicalPartition logical_partition;

  // projection from this color_space into the 1D default color_space
  const static unsigned projection_id = ProjectionIDs::IDENTITY;
  Context ctx;
  Runtime* runtime;
};


struct Data2D {
  Data2D(coord_t n_grid_cells_, coord_t n_second_, coord_t n_parts_,
         const std::string& name_,
         const std::vector<std::string>& field_names_,
         Context ctx_, Runtime* runtime_) :
      name(name_),
      field_names(field_names_),
      n_grid_cells(n_grid_cells_),
      n_second(n_second_),
      domain(Point<2>(0,0), Point<2>(n_grid_cells_-1, n_second_-1)),
      n_parts(n_parts_),
      color_domain(Point<2>(0,0), Point<2>(n_parts_-1,0)),
      ctx(ctx_),
      runtime(runtime_)
  {
    index_space = runtime->create_index_space<2,coord_t>(ctx, domain);

    auto block = Point<2>(n_grid_cells / n_parts, n_second);
    index_partition = runtime->create_partition_by_blockify<2,coord_t>(ctx, index_space, block);
    runtime->attach_name(index_partition, (std::string("index_partition_")+name).c_str());

    field_space = runtime->create_field_space(ctx);
    {
      FieldAllocator allocator = runtime->create_field_allocator(ctx, field_space);
      for (auto fname : field_names_)
        field_ids[fname] = allocator.allocate_field(sizeof(double), AUTO_GENERATE_ID);
    }
    runtime->attach_name(field_space, (std::string("field_space_")+name).c_str());

    logical_region = runtime->create_logical_region(ctx, index_space, field_space);
    runtime->attach_name(logical_region, (std::string("logical_region_")+name).c_str());
    logical_partition = runtime->get_logical_partition(ctx, logical_region, index_partition);
    runtime->attach_name(logical_partition, (std::string("logical_partition_")+name).c_str());
  }

  ~Data2D() {
    runtime->destroy_logical_partition(ctx, logical_partition);
    runtime->destroy_logical_region(ctx, logical_region);
    runtime->destroy_field_space(ctx, field_space);
    runtime->destroy_index_partition(ctx, index_partition);
    runtime->destroy_index_space(ctx, index_space);
  }
  
  std::string name;
  std::vector<std::string> field_names;

  coord_t n_grid_cells, n_second;
  const Rect<2> domain;
  IndexSpaceT<2,coord_t> index_space;
  IndexPartitionT<2, coord_t> index_partition;

  std::map<std::string, unsigned> field_ids;
  FieldSpace field_space;

  LogicalRegion logical_region;
  LogicalPartition logical_partition;

  coord_t n_parts;
  const Rect<2> color_domain;

  Context ctx;
  Runtime* runtime;

  class LocalProjectionFunction : public ProjectionFunctor{
   public:
    LocalProjectionFunction() {}
    LocalProjectionFunction(Runtime *rt) {}
    virtual ~LocalProjectionFunction() {}
   public:
    virtual LogicalRegion project(const Mappable *mappable, unsigned index,
            LogicalRegion upper_bound,
            const DomainPoint &point) override {
      assert(false);
    }
    virtual LogicalRegion project(const Mappable *mappable, unsigned index,
            LogicalPartition upper_bound,
            const DomainPoint &point) override {
      Point<2,coord_t> color; color[0] = point[0]; color[1] = 0;
      return runtime->get_logical_subregion_by_color(upper_bound, color);
    }
    virtual bool is_exclusive(void) const override { return true; }
    virtual unsigned get_depth(void) const override { return 0; }
  };

  const static unsigned projection_id = ProjectionIDs::TO_1D_FROM_2D;

};
  

struct Data2D_Transposed {
  Data2D_Transposed(coord_t n_grid_cells_, coord_t n_second_, coord_t n_parts_,
                    const std::string& name_,
                    const std::vector<std::string>& field_names_,
                    Context ctx_, Runtime* runtime_) :
      name(name_),
      field_names(field_names_),
      n_grid_cells(n_grid_cells_),
      n_second(n_second_),
      domain(Point<2>(0,0), Point<2>(n_second_-1, n_grid_cells_-1)),
      n_parts(n_parts_),
      color_domain(Point<2>(0,0), Point<2>(0,n_parts_-1)),
      ctx(ctx_),
      runtime(runtime_)
  {
    index_space = runtime->create_index_space<2,coord_t>(ctx, domain);

    auto block = Point<2>(n_second, n_grid_cells / n_parts);
    index_partition = runtime->create_partition_by_blockify<2,coord_t>(ctx, index_space, block);
    runtime->attach_name(index_partition, (std::string("index_partition_")+name).c_str());

    field_space = runtime->create_field_space(ctx);
    {
      FieldAllocator allocator = runtime->create_field_allocator(ctx, field_space);
      for (auto fname : field_names_)
        field_ids[fname] = allocator.allocate_field(sizeof(double), AUTO_GENERATE_ID);

      std::cout << "Field IDs (" << name_ << "):";
      for (auto fid : field_ids) {
        std::cout << fid.first << "," << fid.second << ";";
      }
      std::cout << std::endl;
        
    }
    runtime->attach_name(field_space, (std::string("field_space_")+name).c_str());

    logical_region = runtime->create_logical_region(ctx, index_space, field_space);
    runtime->attach_name(logical_region, (std::string("logical_region_")+name).c_str());
    logical_partition = runtime->get_logical_partition(ctx, logical_region, index_partition);
    runtime->attach_name(logical_partition, (std::string("logical_partition_")+name).c_str());
  }

  ~Data2D_Transposed() {
    runtime->destroy_logical_partition(ctx, logical_partition);
    runtime->destroy_logical_region(ctx, logical_region);
    runtime->destroy_field_space(ctx, field_space);
    runtime->destroy_index_partition(ctx, index_partition);
    runtime->destroy_index_space(ctx, index_space);
  }

  
  std::string name;
  std::vector<std::string> field_names;

  coord_t n_grid_cells, n_second;
  const Rect<2> domain;
  IndexSpaceT<2,coord_t> index_space;
  IndexPartitionT<2, coord_t> index_partition;

  std::map<std::string,unsigned> field_ids;
  FieldSpace field_space;

  LogicalRegion logical_region;
  LogicalPartition logical_partition;

  coord_t n_parts;
  const Rect<2> color_domain;

  static const unsigned projection_id = ProjectionIDs::TO_1D_FROM_2D_TRANSPOSED;

  Context ctx;
  Runtime* runtime;
  
  class LocalProjectionFunction : public ProjectionFunctor{
   public:
    LocalProjectionFunction() {}
    LocalProjectionFunction(Runtime *rt) {}
    virtual ~LocalProjectionFunction() {}
   public:
    virtual LogicalRegion project(const Mappable *mappable, unsigned index,
            LogicalRegion upper_bound,
            const DomainPoint &point) override {
      assert(false);
    }
    virtual LogicalRegion project(const Mappable *mappable, unsigned index,
            LogicalPartition upper_bound,
            const DomainPoint &point) override {
      Point<2,coord_t> color; color[0] = 0; color[1] = point[0];
      return runtime->get_logical_subregion_by_color(upper_bound, color);
    }
    virtual bool is_exclusive(void) const override { return true; }
    virtual unsigned get_depth(void) const override { return 0; }
  };
  
};
  

#endif
