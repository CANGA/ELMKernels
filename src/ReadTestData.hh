// Load data from NetCDF initialization file.





namespace ELM {

template<typename ArrayD1>
void ReadLandData(const std::string &dir, const std::string &basename, const Utils::DomainDecomposition<2> &dd, ArrayD1 &topo_slope,
  ArrayD1 &topo_std) {

  IO::read_distributed_scalar<double>(dir, basename, "SLOPE", dd, topo_slope);
  IO::read_distributed_scalar<double>(dir, basename, "STD_ELEV", dd, topo_std);
}

template<typename ArrayD1, typename ArrayD2>
void ReadTestInitData(const std::string &dir, const std::string &basename, const Utils::DomainDecomposition<2> &dd, 
  ArrayD2 &t_soisno, ArrayD1 &snl) {


  IO::read_distributed_array<double>(dir, basename, "T_SOISNO", dd, t_soisno);
  IO::read_distributed_scalar<int>(dir, basename, "SNLSNO", snl);// --  need to multiply by -1
}





} // namespace ELM