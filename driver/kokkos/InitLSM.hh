
#pragma once

#include <InitTopography.hh>
#include <InitSnowLayers.hh>

using ArrayI1 = Kokkos::View<int *>;
using ArrayD1 = Kokkos::View<double *>;
using ArrayD2 = Kokkos::View<double **>;

struct CallInitSnowLayers {

  CallInitSnowLayers (ArrayD1 &snow_depth_, ELM::LandType &Land_, ArrayI1 &snl_, 
                      ArrayD2 &dz_, ArrayD2 &z_,ArrayD2 &zi_) 
  : snow_depth(snow_depth_), Land(Land_), snl(snl_), dz(dz_), z(z_), zi(zi_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {

    ELM::InitSnowLayers(snow_depth[i], Land.lakpoi, snl[i], Kokkos::subview(dz, i, Kokkos::ALL), 
      Kokkos::subview(z, i, Kokkos::ALL), Kokkos::subview(zi, i, Kokkos::ALL));
  }

private:
  ELM::LandType Land;
  ArrayI1 snl;
  ArrayD1 snow_depth; 
  ArrayD2 dz, z, zi;
}; 


struct CallInitTopography {

  CallInitTopography (ELM::LandType &Land_, ArrayD1 &topo_slope_, ArrayD1 &topo_std_, ArrayD1 &n_melt_, 
                      ArrayD1 &micro_sigma_)
  : Land(Land_), topo_slope(topo_slope_), topo_std(topo_std_), n_melt(n_melt_), micro_sigma(micro_sigma_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {

    ELM::InitTopoSlope(topo_slope);
    ELN::InitMicroTopo(Land.ltype, topo_slope, topo_std,
                        n_melt, micro_sigma);
  }

private:
  ELM::LandType Land;
  ArrayD1 topo_slope, topo_std, n_melt, micro_sigma;
}; 
