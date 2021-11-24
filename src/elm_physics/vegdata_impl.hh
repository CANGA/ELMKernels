
#pragma once
#include "vegdata.h"

namespace ELM {

// default constructor - construct array objects in initializer list
// in case ArrayType can't be default constructed
template <typename ArrayD1, typename ArrayD2>
VegData<ArrayD1, ArrayD2>::VegData() :
    fnr("fnr", ELM::numpft), act25("act25", ELM::numpft), kcha("kcha", ELM::numpft),
    koha("koha", ELM::numpft), cpha("cpha", ELM::numpft), vcmaxha("vcmaxha", ELM::numpft),
    jmaxha("jmaxha", ELM::numpft), tpuha("tpuha", ELM::numpft), lmrha("lmrha", ELM::numpft),
    vcmaxhd("vcmaxhd", ELM::numpft), jmaxhd("jmaxhd", ELM::numpft), tpuhd("tpuhd", ELM::numpft),
    lmrhd("lmrhd", ELM::numpft), lmrse("lmrse", ELM::numpft), qe("qe", ELM::numpft),
    theta_cj("theta_cj", ELM::numpft), bbbopt("bbbopt", ELM::numpft), mbbopt("mbbopt", ELM::numpft),
    c3psn("c3psn", ELM::numpft), slatop("slatop", ELM::numpft), leafcn("leafcn", ELM::numpft),
    flnr("flnr", ELM::numpft), fnitr("fnitr", ELM::numpft), dleaf("dleaf", ELM::numpft),
    smpso("smpso", ELM::numpft), smpsc("smpsc", ELM::numpft), tc_stress("tc_stress", 1),
    z0mr("z0mr", ELM::numpft), displar("displar", ELM::numpft), xl("xl", ELM::numpft),
    roota_par("roota_par", ELM::numpft), rootb_par("rootb_par", ELM::numpft), 
    rholvis("rholvis", ELM::numpft), rholnir("rholnir", ELM::numpft),
    rhosvis("rhosvis", ELM::numpft), rhosnir("rhosnir", ELM::numpft),
    taulvis("taulvis", ELM::numpft), taulnir("taulnir", ELM::numpft),
    tausvis("tausvis", ELM::numpft), tausnir("tausnir", ELM::numpft)   { }


template <typename ArrayD1, typename ArrayD2>
void VegData<ArrayD1, ArrayD2>::read_veg_data(const std::string &data_dir, const std::string &basename_pfts) {

  ELM::Array<std::string, 1> pftnames("pftnames", ELM::mxpft);
  const std::array<std::string, ELM::mxpft> expected_pftnames = {"not_vegetated",
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
  ELM::IO::read_pft_var(data_dir, basename_pfts, "c3psn", c3psn);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "slatop", slatop);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "leafcn", leafcn);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "flnr", flnr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "fnitr", fnitr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "dleaf", dleaf);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "smpso", smpso);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "smpsc", smpsc);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tc_stress", tc_stress);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "z0mr", z0mr);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "displar", displar);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "xl", xl);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "roota_par", roota_par);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rootb_par", rootb_par);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rholvis", rholvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rholnir", rholnir);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rhosvis", rhosvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "rhosnir", rhosnir);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "taulvis", taulvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "taulnir", taulnir);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tausvis", tausvis);
  ELM::IO::read_pft_var(data_dir, basename_pfts, "tausnir", tausnir);
}


template <typename ArrayD1, typename ArrayD2>
PSNVegData VegData<ArrayD1, ArrayD2>::get_pft_psnveg (int vegtype) {
  PSNVegData psnvegdata;
  psnvegdata.fnr = fnr(vegtype);
  psnvegdata.act25 = act25(vegtype);
  psnvegdata.kcha = kcha(vegtype);
  psnvegdata.koha = koha(vegtype);
  psnvegdata.cpha = cpha(vegtype);
  psnvegdata.vcmaxha = vcmaxha(vegtype);
  psnvegdata.jmaxha = jmaxha(vegtype);
  psnvegdata.tpuha = tpuha(vegtype);
  psnvegdata.lmrha = lmrha(vegtype);
  psnvegdata.vcmaxhd = vcmaxhd(vegtype);
  psnvegdata.jmaxhd = jmaxhd(vegtype);
  psnvegdata.tpuhd = tpuhd(vegtype);
  psnvegdata.lmrhd = lmrhd(vegtype);
  psnvegdata.lmrse = lmrse(vegtype);
  psnvegdata.qe = qe(vegtype);
  psnvegdata.theta_cj = theta_cj(vegtype);
  psnvegdata.bbbopt = bbbopt(vegtype);
  psnvegdata.mbbopt = mbbopt(vegtype);
  psnvegdata.c3psn = c3psn(vegtype);
  psnvegdata.slatop = slatop(vegtype);
  psnvegdata.leafcn = leafcn(vegtype);
  psnvegdata.flnr = flnr(vegtype);
  psnvegdata.fnitr = fnitr(vegtype);
  psnvegdata.dleaf = dleaf(vegtype);
  psnvegdata.smpso = smpso(vegtype);
  psnvegdata.smpsc = smpsc(vegtype);
  psnvegdata.tc_stress = tc_stress(0);
  return psnvegdata;
}

template <typename ArrayD1, typename ArrayD2>
AlbedoVegData VegData<ArrayD1, ArrayD2>::get_pft_albveg(int vegtype) {
 AlbedoVegData albedovegdata;
 albedovegdata.rhol[0] = rholvis(vegtype);
 albedovegdata.rhol[1] = rholnir(vegtype);
 albedovegdata.rhos[0] = rhosvis(vegtype);
 albedovegdata.rhos[1] = rhosnir(vegtype);
 albedovegdata.taul[0] = taulvis(vegtype);
 albedovegdata.taul[1] = taulnir(vegtype);
 albedovegdata.taus[0] = tausvis(vegtype);
 albedovegdata.taus[1] = tausnir(vegtype);
 albedovegdata.xl = xl(vegtype);
 return albedovegdata;
}




} // namespace ELM
