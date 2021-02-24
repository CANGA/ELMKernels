#pragma once

#include "clm_constants.h"
#include "read_input.hh"
#include "utils.hh"
#include <assert.h>
namespace ELM {

static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

// Read and initialize vegetation (PFT) constants
template <typename ArrayS1, typename ArrayD1>
void ReadPFTConstants(const std::string &data_dir, const std::string &basename_pfts, const int maxpfts,
                      ArrayS1 &pftnames, ArrayD1 &z0mr, ArrayD1 &displar, ArrayD1 &dleaf, ArrayD1 &c3psn) {

  std::array<std::string, 25> expected_pftnames = {"not_vegetated",
                                                   "needleleaf_evergreen_temperate_tree",
                                                   "needleleaf_evergreen_boreal_tree",
                                                   "needleleaf_deciduous_boreal_tree",
                                                   "broadleaf_evergreen_tropical_tree",
                                                   "broadleaf_evergreen_temperate_tree",
                                                   "broadleaf_deciduous_tropical_tree",
                                                   "broadleaf_deciduous_temperate_tree",
                                                   "broadleaf_deciduous_boreal_tree",
                                                   "broadleaf_evergreen_shrub",
                                                   "broadleaf_deciduous_temperate_shrub",
                                                   "broadleaf_deciduous_boreal_shrub",
                                                   "c3_arctic_grass",
                                                   "c3_non-arctic_grass",
                                                   "c4_grass",
                                                   "c3_crop",
                                                   "c3_irrigated",
                                                   "corn",
                                                   "irrigated_corn",
                                                   "spring_temperate_cereal",
                                                   "irrigated_spring_temperate_cereal",
                                                   "winter_temperate_cereal",
                                                   "irrigated_winter_temperate_cereal",
                                                   "soybean",
                                                   "irrigated_soybean"};

  int MPI_COMM_WORLD;

  // rhol, rhos, taul, taus are xxx[pfts][numrad] - need to figure out read into one ArrayD2, or split into 1D

  // const auto maxpfts = ELM::IO::get_maxpfts(MPI_COMM_WORLD, data_dir, basename_pfts, "z0mr");
  // ELM::IO::read_pft(data_dir, basename_pfts, "pftname", pftnames);
  const int strlen = 40;
  ELM::IO::read_names(data_dir, basename_pfts, "pftname", strlen, pftnames);

  for (int i = 0; i != pftnames.extent(0); i++)
    assert(pftnames[i] == expected_pftnames[i] && "pftname does not match expected");

  ELM::IO::read_pft_var(data_dir, basename_pfts, "z0mr", z0mr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "displar", displar);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "dleaf", dleaf);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "c3psn", c3psn);
  // ELM::IO::read_pft_var(data_dir, basename_pfts, "rholvis", rhol[0]);
  // ELM::IO::read_pft_var(data_dir, basename_pfts, "rholnir", rhol[1]);
  // ELM::IO::read_pft_var(data_dir, basename_pfts, "z0mr", z0mr);
}

} // namespace ELM