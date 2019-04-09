#ifndef DOMAINS_HH_
#define DOMAINS_HH_

#include "legion.h"
using namespace Legion;

// projection
template<size_t DIM>
class LocalProjectionFunction : public ProjectionFunctor {
 public:
  LocalProjectionFunction(const Point<DIM>& partition, Runtime* rt) :
      ProjectionFunctor(rt),
      partition_(partition) {}
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
    auto color = Point<DIM,coord_t>::ZEROES();
    size_t j = 0;
    for (size_t i = 0; i!=DIM; ++i) {
      if (partition_[i] > 1) {
        assert(j <= point.get_dim());
        color[i] = point[j];
        j++;
      }
    }
    return runtime->get_logical_subregion_by_color(upper_bound, (DomainPoint) color);
  }
  virtual bool is_exclusive(void) const override { return true; }
  virtual unsigned get_depth(void) const override { return 0; }

 private:
  Point<DIM> partition_;
};



template<size_t DIM>
struct Data {
  Data(const std::string& name_,
         Context ctx_, Runtime* runtime_) :
      name(name_),
      ctx(ctx_),
      runtime(runtime_)
  {
    field_space = runtime->create_field_space(ctx);
    allocator_ = runtime->create_field_allocator(ctx, field_space);
  }

  Data(const std::string& name_,
       Context ctx_, Runtime* runtime_,
       const Point<DIM> shape_,
       const Point<DIM> blocked_shape_) :
      Data(name_, ctx_, runtime_)
  {
    initShape(shape_, blocked_shape_);
  }

  Data(const std::string& name_,
       Context ctx_, Runtime* runtime_,
       const Point<DIM> shape_,
       const Point<DIM> blocked_shape_,
       const std::vector<std::string>& field_names_) :
      Data(name_, ctx_, runtime_, shape_, blocked_shape_)
  {
    addFields<double>(field_names_);
    finalize();
  }

  // member functions
  void initShape(const Point<DIM> shape_, const Point<DIM> blocked_shape_)
  {
    shape = shape_;
    blocked_shape = blocked_shape_;

    domain = Rect<DIM,coord_t>(Point<DIM>::ZEROES(), shape - Point<DIM>::ONES());
    index_space = runtime->create_index_space<DIM>(ctx, domain);

    color_domain = Rect<DIM>(Point<DIM>::ZEROES(), blocked_shape - Point<DIM>::ONES());
    Point<DIM> block;
    for (int i = 0; i!=DIM; ++i) block[i] = shape[i] / blocked_shape[i];
    index_partition = runtime->create_partition_by_blockify<DIM>(ctx, index_space, block);
    runtime->attach_name(index_partition, (std::string("index_partition_")+name).c_str());
  }

  template<typename T>
  FieldID allocateField_() {
    return allocator_.allocate_field(sizeof(T), AUTO_GENERATE_ID);
  }
  
  template<typename T> 
  void addField(const std::string& fieldname) {
    FieldID fid = allocateField_<T>();
    field_ids[fieldname] = fid;
    field_names.push_back(fieldname);
  }

  template<typename T> 
  void addFields(const std::vector<std::string>& fieldnames) {
    for (auto fname : fieldnames) {
      addField<T>(fname);
    }
  }

  template<typename... Ts>
  void addVariadicFields(const std::vector<std::string>& fieldnames) {
    // unpack the template parameter, adding the right types
    FieldID fids[]{allocateField_<Ts>()...};
    size_t i = 0;
    for (auto fname : fieldnames) {
      field_ids[fname] = fids[i];
      field_names.push_back(fname);
      i++;
    }
  }

  void finalize() {
    allocator_ = FieldAllocator();
    runtime->attach_name(field_space, (std::string("field_space_")+name).c_str());
    
    // // DEBUG CRUFT
    // std::vector<FieldID> my_fields;
    // runtime->get_field_space_fields(ctx, field_space, my_fields);
    // std::cout << "Allocated field space of size: " << my_fields.size() << std::endl;
    // std::cout << "num ids: " << field_ids.size() << std::endl;
    // std::cout << "num names: " << field_names.size() << std::endl;
    // // END DEBUG CRUFT
    
    // create logical regions, partitions
    logical_region = runtime->create_logical_region(ctx, index_space, field_space);
    runtime->attach_name(logical_region, (std::string("logical_region_")+name).c_str());
    logical_partition = runtime->get_logical_partition(ctx, logical_region, index_partition);
    runtime->attach_name(logical_partition, (std::string("logical_partition_")+name).c_str());

    // register the projection
    projection_id = runtime->generate_dynamic_projection_id();
    auto* pf = new LocalProjectionFunction<DIM>(blocked_shape, runtime);
    runtime->register_projection_functor(projection_id,pf);
            
  }

  ~Data() {
    runtime->destroy_logical_partition(ctx, logical_partition);
    runtime->destroy_logical_region(ctx, logical_region);
    runtime->destroy_field_space(ctx, field_space);
    runtime->destroy_index_partition(ctx, index_partition);
    runtime->destroy_index_space(ctx, index_space);
  }

  // name the container
  std::string name;
  Context ctx;
  Runtime* runtime;

  // index space
  Point<DIM> shape;
  Point<DIM> blocked_shape;
  Rect<DIM> domain;
  IndexSpaceT<DIM> index_space;

  // index partition
  coord_t n_parts;
  Rect<DIM> color_domain;
  IndexPartitionT<DIM, coord_t> index_partition;

  // field space
  std::vector<std::string> field_names;
  std::map<std::string, FieldID> field_ids;
  FieldSpace field_space;
  FieldAllocator allocator_;

  // cross product of index and field
  LogicalRegion logical_region;
  LogicalPartition logical_partition;

  ProjectionID projection_id;
};

#endif
