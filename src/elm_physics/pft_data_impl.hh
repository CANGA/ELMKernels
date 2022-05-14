
#pragma once

namespace ELM {

// default constructor - construct array objects in initializer list
// in case ArrayType can't be default constructed
template <typename ArrayD1, typename ArrayD2>
PFTData<ArrayD1, ArrayD2>::PFTData()
    : fnr("fnr", ELM::numpft), act25("act25", ELM::numpft), kcha("kcha", ELM::numpft), koha("koha", ELM::numpft),
      cpha("cpha", ELM::numpft), vcmaxha("vcmaxha", ELM::numpft), jmaxha("jmaxha", ELM::numpft),
      tpuha("tpuha", ELM::numpft), lmrha("lmrha", ELM::numpft), vcmaxhd("vcmaxhd", ELM::numpft),
      jmaxhd("jmaxhd", ELM::numpft), tpuhd("tpuhd", ELM::numpft), lmrhd("lmrhd", ELM::numpft),
      lmrse("lmrse", ELM::numpft), qe("qe", ELM::numpft), theta_cj("theta_cj", ELM::numpft),
      bbbopt("bbbopt", ELM::numpft), mbbopt("mbbopt", ELM::numpft), c3psn("c3psn", ELM::numpft),
      slatop("slatop", ELM::numpft), leafcn("leafcn", ELM::numpft), flnr("flnr", ELM::numpft),
      fnitr("fnitr", ELM::numpft), dleaf("dleaf", ELM::numpft), smpso("smpso", ELM::numpft),
      smpsc("smpsc", ELM::numpft), tc_stress("tc_stress", 1), z0mr("z0mr", ELM::numpft),
      displar("displar", ELM::numpft), xl("xl", ELM::numpft), roota_par("roota_par", ELM::numpft),
      rootb_par("rootb_par", ELM::numpft), rholvis("rholvis", ELM::numpft), rholnir("rholnir", ELM::numpft),
      rhosvis("rhosvis", ELM::numpft), rhosnir("rhosnir", ELM::numpft), taulvis("taulvis", ELM::numpft),
      taulnir("taulnir", ELM::numpft), tausvis("tausvis", ELM::numpft), tausnir("tausnir", ELM::numpft)
    {}

template <typename ArrayD1, typename ArrayD2>
ACCELERATED
PFTDataPSN PFTData<ArrayD1, ArrayD2>::get_pft_psn(const int pft) const {
  PFTDataPSN psn_pft_data;
  psn_pft_data.fnr = fnr(pft);
  psn_pft_data.act25 = act25(pft);
  psn_pft_data.kcha = kcha(pft);
  psn_pft_data.koha = koha(pft);
  psn_pft_data.cpha = cpha(pft);
  psn_pft_data.vcmaxha = vcmaxha(pft);
  psn_pft_data.jmaxha = jmaxha(pft);
  psn_pft_data.tpuha = tpuha(pft);
  psn_pft_data.lmrha = lmrha(pft);
  psn_pft_data.vcmaxhd = vcmaxhd(pft);
  psn_pft_data.jmaxhd = jmaxhd(pft);
  psn_pft_data.tpuhd = tpuhd(pft);
  psn_pft_data.lmrhd = lmrhd(pft);
  psn_pft_data.lmrse = lmrse(pft);
  psn_pft_data.qe = qe(pft);
  psn_pft_data.theta_cj = theta_cj(pft);
  psn_pft_data.bbbopt = bbbopt(pft);
  psn_pft_data.mbbopt = mbbopt(pft);
  psn_pft_data.c3psn = c3psn(pft);
  psn_pft_data.slatop = slatop(pft);
  psn_pft_data.leafcn = leafcn(pft);
  psn_pft_data.flnr = flnr(pft);
  psn_pft_data.fnitr = fnitr(pft);
  psn_pft_data.dleaf = dleaf(pft);
  psn_pft_data.smpso = smpso(pft);
  psn_pft_data.smpsc = smpsc(pft);
  psn_pft_data.tc_stress = tc_stress(0);
  return psn_pft_data;
}

template <typename ArrayD1, typename ArrayD2>
ACCELERATED
PFTDataAlb PFTData<ArrayD1, ArrayD2>::get_pft_alb(const int pft) const {
  PFTDataAlb alb_pft_data;
  alb_pft_data.rhol[0] = rholvis(pft);
  alb_pft_data.rhol[1] = rholnir(pft);
  alb_pft_data.rhos[0] = rhosvis(pft);
  alb_pft_data.rhos[1] = rhosnir(pft);
  alb_pft_data.taul[0] = taulvis(pft);
  alb_pft_data.taul[1] = taulnir(pft);
  alb_pft_data.taus[0] = tausvis(pft);
  alb_pft_data.taus[1] = tausnir(pft);
  alb_pft_data.xl = xl(pft);
  return alb_pft_data;
}


template <typename h_ArrayD1>
void read_pft_data(std::map<std::string, h_ArrayD1>& pft_views,
                                              const Comm_type& comm, const std::string& fname_pft)
{

  ELM::Array<std::string, 1> pftnames("pftnames", ELM::mxpft);
  const std::array<std::string, ELM::mxpft> expected_pftnames = 
    { "not_vegetated",
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
      "irrigated_soybean"
    };

  // read pftnames
  const int strlen = 40;
  ELM::IO::read_names(comm, fname_pft, "pftname", strlen, pftnames);

  // check to make sure order is as expected
  for (int i = 0; i != pftnames.extent(0); i++)
    assert(pftnames[i] == expected_pftnames[i] && "pftname does not match expected pftname");

  // read pft constants
  for (auto& [varname, arr] : pft_views) {
    ELM::IO::read_pft_var(comm, fname_pft, varname, arr);
  }
}

} // namespace ELM
