
#include "aerosol_data.h"

std::pair<size_t, size_t>
ELM::aerosol_utils::get_nearest_indices(const Comm_type& comm,
                                        const std::string& filename,
                                        const double& lon_d, const double& lat_d)
{
  Array<double, 1> file_lon(144);
  Array<double, 1> file_lat(96);
  const std::array<size_t, 1> start{0};
  const std::array<size_t, 1> count_lon{144};
  const std::array<size_t, 1> count_lat{96};
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
