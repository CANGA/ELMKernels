#pragma once

#include "ELMConstants.h"
#include "read_input.hh"
#include "utils.hh"
#include <assert.h>
namespace ELM {

static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

// Read and initialize vegetation (PFT) constants
template <typename ArrayS1, typename ArrayD1>
void ReadPFTConstants(const std::string &data_dir, const std::string &basename_pfts, ArrayS1 &pftnames, ArrayD1 &z0mr,
                      ArrayD1 &displar, ArrayD1 &dleaf, ArrayD1 &c3psn, ArrayD1 &xl, ArrayD1 &roota_par,
                      ArrayD1 &rootb_par, ArrayD1 &slatop, ArrayD1 &leafcn, ArrayD1 &flnr, ArrayD1 &smpso,
                      ArrayD1 &smpsc, ArrayD1 &fnitr, ArrayD1 &fnr, ArrayD1 &act25, ArrayD1 &kcha, ArrayD1 &koha,
                      ArrayD1 &cpha, ArrayD1 &vcmaxha, ArrayD1 &jmaxha, ArrayD1 &tpuha, ArrayD1 &lmrha,
                      ArrayD1 &vcmaxhd, ArrayD1 &jmaxhd, ArrayD1 &tpuhd, ArrayD1 &lmrhd, ArrayD1 &lmrse, ArrayD1 &qe,
                      ArrayD1 &theta_cj, ArrayD1 &bbbopt, ArrayD1 &mbbopt, ArrayD1 &nstor, ArrayD1 &br_xr,
                      ArrayD1 &tc_stress, ArrayD1 rholvis, ArrayD1 rholnir, ArrayD1 rhosvis, ArrayD1 rhosnir,
                      ArrayD1 taulvis, ArrayD1 taulnir, ArrayD1 tausvis, ArrayD1 tausnir) {

  const std::array<std::string, 25> expected_pftnames = {"not_vegetated",
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

  // read pftnames
  const int strlen = 40;
  ELM::IO::read_names(data_dir, basename_pfts, "pftname", strlen, pftnames);

  // check to make sure order is as expected
  for (int i = 0; i != pftnames.extent(0); i++)
    assert(pftnames[i] == expected_pftnames[i] && "pftname does not match expected pftname");

  // read pft constants
  ELM::IO::read_pft_var(data_dir, basename_pfts, "z0mr", z0mr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "displar", displar);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "dleaf", dleaf);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "c3psn", c3psn);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "xl", xl);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "roota_par", roota_par);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rootb_par", rootb_par);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "slatop", slatop);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "leafcn", leafcn);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "flnr", flnr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "smpso", smpso);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "smpsc", smpsc);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "fnitr", fnitr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "fnr", fnr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "act25", act25);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "kcha", kcha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "koha", koha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "cpha", cpha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "vcmaxha", vcmaxha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "jmaxha", jmaxha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tpuha", tpuha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "lmrha", lmrha);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "vcmaxhd", vcmaxhd);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "jmaxhd", jmaxhd);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tpuhd", tpuhd);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "lmrhd", lmrhd);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "lmrse", lmrse);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "qe", qe);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "theta_cj", theta_cj);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "bbbopt", bbbopt);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "mbbopt", mbbopt);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "nstor", nstor);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "br_xr", br_xr);

  ELM::IO::read_pft_var(data_dir, basename_pfts, "rholvis", rholvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rholnir", rholnir);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rhosvis", rhosvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rhosnir", rhosnir);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "taulvis", taulvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "taulnir", taulnir);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tausvis", tausvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tausnir", tausnir);

  ELM::IO::read_pft_var(data_dir, basename_pfts, "tc_stress", tc_stress);
}

} // namespace ELM
