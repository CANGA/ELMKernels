
#include "aerosol_data.h"
#include "aerosol_physics.h"
#include "elm_constants.h"
#include "monthly_data.h"
#include "read_input.hh"

#include <array>
#include <cmath>
#include <functional>
#include <tuple>

#include "invoke_kernel.hh"

ELM::AerosolMasses::AerosolMasses(const int ncells)
    : mss_bcphi("mss_bcphi", ncells, ELM::nlevsno, 0.0), mss_bcpho("mss_bcpho", ncells, ELM::nlevsno, 0.0),
      mss_dst1("mss_dst1", ncells, ELM::nlevsno, 0.0), mss_dst2("mss_dst2", ncells, ELM::nlevsno, 0.0),
      mss_dst3("mss_dst3", ncells, ELM::nlevsno, 0.0), mss_dst4("mss_dst4", ncells, ELM::nlevsno, 0.0) {}

ELM::AerosolConcentrations::AerosolConcentrations(const int ncells)
    : mss_cnc_bcphi("mss_cnc_bcphi", ncells, ELM::nlevsno, 0.0),
      mss_cnc_bcpho("mss_cnc_bcpho", ncells, ELM::nlevsno, 0.0),
      mss_cnc_dst1("mss_cnc_dst1", ncells, ELM::nlevsno, 0.0), mss_cnc_dst2("mss_cnc_dst2", ncells, ELM::nlevsno, 0.0),
      mss_cnc_dst3("mss_cnc_dst3", ncells, ELM::nlevsno, 0.0), mss_cnc_dst4("mss_cnc_dst4", ncells, ELM::nlevsno, 0.0) {
}

ELM::AerosolDataManager::AerosolDataManager()
    : bcdep_("bcdep", 12), bcpho_("bcpho", 12), bcphi_("bcphi", 12), dst1_1_("dst1_1", 12), dst1_2_("dst1_2", 12),
      dst2_1_("dst2_1", 12), dst2_2_("dst2_2", 12), dst3_1_("dst3_1", 12), dst3_2_("dst3_2", 12),
      dst4_1_("dst4_1", 12), dst4_2_("dst4_2", 12) {}

void ELM::AerosolDataManager::read_variable_slice(const Comm_type& comm, const std::string& filename,
                                                  const std::string& varname, const size_t& lon_idx,
                                                  const size_t& lat_idx, h_ArrayD1& arr) {
  Array<double, 3> file_data(12, 1, 1); // one year of monthly data
  std::array<size_t, 3> start{0, lat_idx, lon_idx};
  std::array<size_t, 3> count{12, 1, 1};
  IO::read_netcdf(comm, filename, varname, start, count, file_data.data());

  for (int i = 0; i < 12; ++i) {
    arr(i) = file_data(i, 0, 0);
  }
}

void ELM::AerosolDataManager::read_data(std::map<std::string, ArrayD1::HostMirror>& aerosol_views, 
  const Comm_type& comm, const std::string& filename, const double& lon_d, const double& lat_d) {

  auto [lon_idx, lat_idx] = get_nearest_indices(comm, filename, lon_d, lat_d);

  for (auto& [varname, arr] : aerosol_views) {
    read_variable_slice(comm, filename, varname, lon_idx, lat_idx, arr);
  }
}

std::pair<size_t, size_t> ELM::AerosolDataManager::get_nearest_indices(const Comm_type& comm,
                                                                       const std::string& filename, const double& lon_d,
                                                                       const double& lat_d) {
  Array<double, 1> file_lon(144);
  Array<double, 1> file_lat(96);
  std::array<size_t, 1> start{0};
  std::array<size_t, 1> count_lon{144};
  std::array<size_t, 1> count_lat{96};
  IO::read_netcdf(comm, filename, "lon", start, count_lon, file_lon.data());
  IO::read_netcdf(comm, filename, "lat", start, count_lat, file_lat.data());
  double mindist = 99999.0;
  size_t lon_idx;
  size_t lat_idx;
  for (size_t thisx = 0; thisx < 144; ++thisx) {
    for (size_t thisy = 0; thisy < 96; ++thisy) {

      if (lon_d < 0.0) {
        if (file_lon(thisx) >= 180.0) {
          file_lon(thisx) -= 360.0;
        }
      } else if (lon_d >= 180.0) {
        if (file_lon(thisx) < 0.0)
          file_lon(thisx) += 360.0;
      }

      double dlon2 = pow(file_lon(thisx) - lon_d, 2.0);
      double dlat2 = pow(file_lat(thisy) - lat_d, 2.0);
      double thisdist = 100.0 * pow(dlon2 + dlat2, 0.5);

      if (thisdist < mindist) {
        mindist = thisdist;
        lon_idx = thisx;
        lat_idx = thisy;
      }
    }
  }
  return std::make_pair(lon_idx, lat_idx);
}

auto ELM::AerosolDataManager::get_aerosol_source(const Utils::Date& model_time, const double& dtime) {
  auto [wt1, wt2] = monthly_data::monthly_data_weights(model_time);
  auto [m1, m2] = monthly_data::month_indices(model_time);
  double forc_bcphi = (wt1 * bcphi_(m1) + wt2 * bcphi_(m2)) * dtime;
  double forc_bcpho = (wt1 * (bcdep_(m1) + bcpho_(m1)) + wt2 * (bcdep_(m2) + bcpho_(m2))) * dtime;
  double forc_dst1 = (wt1 * (dst1_1_(m1) + dst1_2_(m1)) + wt2 * (dst1_1_(m2) + dst1_2_(m2))) * dtime;
  double forc_dst2 = (wt1 * (dst2_1_(m1) + dst2_2_(m1)) + wt2 * (dst2_1_(m2) + dst2_2_(m2))) * dtime;
  double forc_dst3 = (wt1 * (dst3_1_(m1) + dst3_2_(m1)) + wt2 * (dst3_1_(m2) + dst3_2_(m2))) * dtime;
  double forc_dst4 = (wt1 * (dst4_1_(m1) + dst4_2_(m1)) + wt2 * (dst4_1_(m2) + dst4_2_(m2))) * dtime;
  return std::make_tuple(forc_bcphi, forc_bcpho, forc_dst1, forc_dst2, forc_dst3, forc_dst4);
}

void ELM::AerosolDataManager::invoke_aerosol_source(const Utils::Date& model_time, const double& dtime,
                                                    const ArrayI1& snl, AerosolMasses& aerosol_masses) {
  auto aerosol_forc_flux = get_aerosol_source(model_time, dtime);
  aerosols::ComputeAerosolDeposition aerosol_source_object(aerosol_forc_flux, snl, aerosol_masses);
  
  const std::string name("ComputeAerosolDeposition");
  invoke_kernel(aerosol_source_object, std::make_tuple(snl.extent(0)), name);

}
