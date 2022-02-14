
#pragma once

#include <algorithm>

namespace ELM::phenology {

// functor to calculate phenology parameters for time = model_time
template <typename ArrayI1, typename ArrayD1, typename ArrayD2> struct ComputePhenology {
  ComputePhenology(const ArrayD2& mlai, const ArrayD2& msai, const ArrayD2& mhtop, const ArrayD2& mhbot,
                   const ArrayD1& snow_depth, const ArrayD1& frac_sno, const ArrayI1& vtype, const double& wt1,
                   const double wt2, const int start_idx, ArrayD1& elai, ArrayD1& esai, ArrayD1& htop, ArrayD1& hbot,
                   ArrayD1& tlai, ArrayD1& tsai, ArrayI1& frac_veg_nosno_alb);

  void operator()(const int i) const;

private:
  ArrayD2 mlai_, msai_, mhtop_, mhbot_;
  ArrayD1 snow_depth_, frac_sno_;
  ArrayI1 vtype_;
  double wt1_, wt2_;
  int start_idx_;
  ArrayD1 elai_, esai_, htop_, hbot_, tlai_, tsai_;
  ArrayI1 frac_veg_nosno_alb_;
};

} // namespace ELM::phenology

#include "phenology_physics_impl.hh"
