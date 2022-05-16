
#pragma once

// this is derived from SatellitePhenologyMod.F90
namespace ELM::phenology {

// functor to calculate phenology parameters for time = model_time
template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
ComputePhenology<ArrayI1, ArrayD1, ArrayD2>::ComputePhenology(const ArrayD2 mlai, const ArrayD2 msai,
                                                              const ArrayD2 mhtop, const ArrayD2 mhbot,
                                                              const ArrayD1 snow_depth, const ArrayD1 frac_sno,
                                                              const ArrayI1 vtype, const double wt1, const double wt2,
                                                              const int start_idx, ArrayD1 elai, ArrayD1 esai,
                                                              ArrayD1 htop, ArrayD1 hbot, ArrayD1 tlai,
                                                              ArrayD1 tsai, ArrayI1 frac_veg_nosno_alb)
    : mlai_(mlai), msai_(msai), mhtop_(mhtop), mhbot_(mhbot), snow_depth_(snow_depth), frac_sno_(frac_sno),
      vtype_(vtype), wt1_(wt1), wt2_(wt2), start_idx_(start_idx), elai_(elai), esai_(esai), htop_(htop), hbot_(hbot),
      tlai_(tlai), tsai_(tsai), frac_veg_nosno_alb_(frac_veg_nosno_alb) {}

template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
ACCELERATE
void ComputePhenology<ArrayI1, ArrayD1, ArrayD2>::operator()(const int i) const {
  // leaf phenology
  // Set leaf and stem areas based on day of year
  // Interpolate leaf area index, stem area index, and vegetation heights
  // between two monthly values using weights, (wt1, wt2)
  if (vtype_(i) != PFT::noveg) {
    tlai_(i) = wt1_ * mlai_(start_idx_, i) + wt2_ * mlai_(start_idx_ + 1, i);
    tsai_(i) = wt1_ * msai_(start_idx_, i) + wt2_ * msai_(start_idx_ + 1, i);
    htop_(i) = wt1_ * mhtop_(start_idx_, i) + wt2_ * mhtop_(start_idx_ + 1, i);
    hbot_(i) = wt1_ * mhbot_(start_idx_, i) + wt2_ * mhbot_(start_idx_ + 1, i);
  } else {
    tlai_(i) = 0.0;
    tsai_(i) = 0.0;
    htop_(i) = 0.0;
    hbot_(i) = 0.0;
  }

  // adjust lai and sai for burying by snow. if exposed lai and sai
  // are less than 0.05, set equal to zero to prevent numerical
  // problems associated with very small lai and sai.
  // snow burial fraction for short vegetation (e.g. grasses) as in
  // Wang and Zeng, 2007.
  double fb;
  if (vtype_(i) > PFT::noveg && vtype_(i) <= PFT::nbrdlf_dcd_brl_shrub) {
    double ol = std::min(std::max(snow_depth_(i) - hbot_(i), 0.0), htop_(i) - hbot_(i));
    fb = 1.0 - ol / std::max(1.e-06, htop_(i) - hbot_(i));
  } else {
    // 0.2m is assumed depth of snow required for complete burial of grasses
    fb = 1.0 - std::max(std::min(snow_depth_(i), 0.2), 0.0) / 0.2;
  }

  // area weight by snow covered fraction
  elai_(i) = std::max(tlai_(i) * (1.0 - frac_sno_(i)) + tlai_(i) * fb * frac_sno_(i), 0.0);
  esai_(i) = std::max(tsai_(i) * (1.0 - frac_sno_(i)) + tsai_(i) * fb * frac_sno_(i), 0.0);
  if (elai_(i) < 0.05) {
    elai_(i) = 0.0;
  }
  if (esai_(i) < 0.05) {
    esai_(i) = 0.0;
  }
  // Fraction of vegetation free of snow
  if ((elai_(i) + esai_(i)) >= 0.05) {
    frac_veg_nosno_alb_(i) = 1;
  } else {
    frac_veg_nosno_alb_(i) = 0;
  }
}

} // namespace ELM::phenology
