// Load data from NetCDF initialization file.





namespace ELM {

template<typename ArrayD1>
void ReadLandData(const std::string &dir, const std::string &basename, const Utils::DomainDecomposition<2> &dd, ArrayD1 &topo_slope,
  ArrayD1 &topo_std) {

  IO::read_distributed_scalar<double>(dir, basename, "SLOPE", dd, topo_slope);
  IO::read_distributed_scalar<double>(dir, basename, "STD_ELEV", dd, topo_std);
}

template <typename ArrayD1, typename ArrayI1, typename ArrayD2>
void ReadTestInitData(const std::string &dir, const std::string &basename, const Utils::DomainDecomposition<2> &dd,
                      ArrayD2 &t_soisno, ArrayI1 &snl, ArrayD1 &snow_depth, ArrayD1 &h2ocan, ArrayD2 &h2osoi_liq,
                      ArrayD2 &h2osoi_ice, ArrayD1 &h2osno, ArrayD2 &snw_rds, ArrayD1 &int_snow) {

  IO::read_distributed_array<double>(dir, basename, "T_SOISNO", dd, t_soisno);
  IO::read_distributed_scalar<int>(dir, basename, "SNLSNO", snl);
  for (auto &sl : snl)
    sl *= -1; // --  need to multiply by -1

  IO::read_distributed_scalar<double>(dir, basename, "SNOW_DEPTH", snow_depth);
  IO::read_distributed_scalar<double>(dir, basename, "H2OCAN", 2, h2ocan); // hardwired for pft == 2

  IO::read_distributed_array<double>(dir, basename, "H2OSOI_LIQ", dd, h2osoi_liq);
  IO::read_distributed_array<double>(dir, basename, "H2OSOI_ICE", dd, h2osoi_ice);

  IO::read_distributed_scalar<double>(dir, basename, "H2OSNO", h2osno);
  IO::read_distributed_array<double>(dir, basename, "snw_rds", dd, snw_rds);
  IO::read_distributed_scalar<double>(dir, basename, "INT_SNOW", int_snow);
}

} // namespace ELM