
#pragma once

template <typename ArrayD2>
ELM::AerosolMasses<ArrayD2>::AerosolMasses(const size_t& ncells)
    : mss_bcphi("mss_bcphi", ncells, nlevsno_),
      mss_bcpho("mss_bcpho", ncells, nlevsno_),
      mss_dst1("mss_dst1", ncells, nlevsno_),
      mss_dst2("mss_dst2", ncells, nlevsno_),
      mss_dst3("mss_dst3", ncells, nlevsno_),
      mss_dst4("mss_dst4", ncells, nlevsno_)
    {}

template <typename ArrayD2>
ELM::AerosolConcentrations<ArrayD2>::AerosolConcentrations(const size_t& ncells)
    : mss_cnc_bcphi("mss_cnc_bcphi", ncells, nlevsno_),
      mss_cnc_bcpho("mss_cnc_bcpho", ncells, nlevsno_),
      mss_cnc_dst1("mss_cnc_dst1", ncells, nlevsno_),
      mss_cnc_dst2("mss_cnc_dst2", ncells, nlevsno_),
      mss_cnc_dst3("mss_cnc_dst3", ncells, nlevsno_),
      mss_cnc_dst4("mss_cnc_dst4", ncells, nlevsno_)
    {}

namespace ELM::aerosols {

template <typename T, typename ArrayI1, typename ArrayD2>
ComputeAerosolDeposition<T, ArrayI1, ArrayD2>::
ComputeAerosolDeposition(const T& aerosol_forc,
                         const ArrayI1 snl,
                         AerosolMasses<ArrayD2>& aerosol_masses)
    : snl_{snl}, aerosol_masses_{aerosol_masses}
    {
      const auto& [forc_bcphi, forc_bcpho, forc_dst1, forc_dst2, forc_dst3, forc_dst4] = aerosol_forc;
      forc_bcphi_ = forc_bcphi;
      forc_bcpho_ = forc_bcpho;
      forc_dst1_ = forc_dst1;
      forc_dst2_ = forc_dst2;
      forc_dst3_ = forc_dst3;
      forc_dst4_ = forc_dst4;
    }

template <typename T, typename ArrayI1, typename ArrayD2>
ACCELERATE
void ComputeAerosolDeposition<T, ArrayI1, ArrayD2>::
operator()(const int i) const {
  if (snl_(i) > 0) {
    const int j = ELMdims::nlevsno() - snl_(i);
    aerosol_masses_.mss_bcphi(i, j) += forc_bcphi_;
    aerosol_masses_.mss_bcpho(i, j) += forc_bcpho_;
    aerosol_masses_.mss_dst1(i, j) += forc_dst1_;
    aerosol_masses_.mss_dst2(i, j) += forc_dst2_;
    aerosol_masses_.mss_dst3(i, j) += forc_dst3_;
    aerosol_masses_.mss_dst4(i, j) += forc_dst4_;
  }
}

template <typename ArrayB1, typename ArrayI1, typename ArrayD1, typename ArrayD2>
ComputeAerosolConcenAndMass<ArrayB1, ArrayI1, ArrayD1, ArrayD2>::
ComputeAerosolConcenAndMass(const double& dtime, const ArrayB1 do_capsnow,
                            const ArrayI1 snl, const ArrayD2 h2osoi_liq,
                            const ArrayD2 h2osoi_ice, const ArrayD2 snw_rds,
                            const ArrayD1 qflx_snwcp_ice, AerosolMasses<ArrayD2>& aerosol_masses,
                            AerosolConcentrations<ArrayD2>& aerosol_concentrations)
    : dtime_{dtime}, do_capsnow_{do_capsnow}, snl_{snl}, h2osoi_liq_{h2osoi_liq}, h2osoi_ice_{h2osoi_ice},
      snw_rds_{snw_rds}, qflx_snwcp_ice_{qflx_snwcp_ice}, aerosol_masses_{aerosol_masses},
      aerosol_concentrations_{aerosol_concentrations} {}

template <typename ArrayB1, typename ArrayI1, typename ArrayD1, typename ArrayD2>
ACCELERATE
void ComputeAerosolConcenAndMass<ArrayB1, ArrayI1, ArrayD1, ArrayD2>::
operator()(const int i) const {
  for (int sl = 0; sl < ELMdims::nlevsno(); ++sl) {
    double snowmass = h2osoi_ice_(i, sl) + h2osoi_liq_(i, sl);
    if (sl == ELMdims::nlevsno() - snl_(i) && do_capsnow_(i)) {
      const double snowcap_scl_fct = snowmass / (snowmass + qflx_snwcp_ice_(i) * dtime_);
      aerosol_masses_.mss_bcpho(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_bcphi(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst1(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst2(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst3(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst4(i, sl) *= snowcap_scl_fct;
    }

    if (sl < ELMdims::nlevsno() - snl_(i)) {
      snw_rds_(i, sl) = 0.0;
      aerosol_masses_.mss_bcpho(i, sl) = 0.0;
      aerosol_masses_.mss_bcphi(i, sl) = 0.0;
      aerosol_masses_.mss_dst1(i, sl) = 0.0;
      aerosol_masses_.mss_dst2(i, sl) = 0.0;
      aerosol_masses_.mss_dst3(i, sl) = 0.0;
      aerosol_masses_.mss_dst4(i, sl) = 0.0;
      snowmass = 1.e-12; // protect against divide by 0
    }

    aerosol_concentrations_.mss_cnc_bcphi(i, sl) = aerosol_masses_.mss_bcphi(i, sl) / snowmass;
    aerosol_concentrations_.mss_cnc_bcpho(i, sl) = aerosol_masses_.mss_bcpho(i, sl) / snowmass;
    aerosol_concentrations_.mss_cnc_dst1(i, sl) = aerosol_masses_.mss_dst1(i, sl) / snowmass;
    aerosol_concentrations_.mss_cnc_dst2(i, sl) = aerosol_masses_.mss_dst2(i, sl) / snowmass;
    aerosol_concentrations_.mss_cnc_dst3(i, sl) = aerosol_masses_.mss_dst3(i, sl) / snowmass;
    aerosol_concentrations_.mss_cnc_dst4(i, sl) = aerosol_masses_.mss_dst4(i, sl) / snowmass;
  }
}

} // namespace ELM::aerosols
