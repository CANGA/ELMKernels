
#pragma once

template <typename ArrayD2>
ELM::AerosolMasses<ArrayD2>::AerosolMasses(const size_t& ncells)
    : mss_bcphi("mss_bcphi", ncells, ELM::nlevsno),
      mss_bcpho("mss_bcpho", ncells, ELM::nlevsno),
      mss_dst1("mss_dst1", ncells, ELM::nlevsno),
      mss_dst2("mss_dst2", ncells, ELM::nlevsno),
      mss_dst3("mss_dst3", ncells, ELM::nlevsno),
      mss_dst4("mss_dst4", ncells, ELM::nlevsno)
    {}

template <typename ArrayD2>
ELM::AerosolConcentrations<ArrayD2>::AerosolConcentrations(const size_t& ncells)
    : mss_cnc_bcphi("mss_cnc_bcphi", ncells, ELM::nlevsno),
      mss_cnc_bcpho("mss_cnc_bcpho", ncells, ELM::nlevsno),
      mss_cnc_dst1("mss_cnc_dst1", ncells, ELM::nlevsno),
      mss_cnc_dst2("mss_cnc_dst2", ncells, ELM::nlevsno),
      mss_cnc_dst3("mss_cnc_dst3", ncells, ELM::nlevsno),
      mss_cnc_dst4("mss_cnc_dst4", ncells, ELM::nlevsno)
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
    const int j = ELM::nlevsno - snl_(i);
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
  for (int sl = 0; sl < ELM::nlevsno; ++sl) {
    double snowmass = h2osoi_ice_(i, sl) + h2osoi_liq_(i, sl);
    if (sl == ELM::nlevsno - snl_(i) && do_capsnow_(i)) {
      const double snowcap_scl_fct = snowmass / (snowmass + qflx_snwcp_ice_(i) * dtime_);
      aerosol_masses_.mss_bcpho(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_bcphi(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst1(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst2(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst3(i, sl) *= snowcap_scl_fct;
      aerosol_masses_.mss_dst4(i, sl) *= snowcap_scl_fct;
    }

    if (sl < ELM::nlevsno - snl_(i)) {
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

template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
void invoke_aerosol_source(const Utils::Date& model_time, const double& dtime, const ArrayI1 snl,
                           const AerosolDataManager<ArrayD1>& aerosol_data,
                           AerosolMasses<ArrayD2>& aerosol_masses)
{
  auto aerosol_forc_flux = aerosol_data.get_aerosol_source(model_time, dtime);
  ComputeAerosolDeposition aerosol_source_object(aerosol_forc_flux, snl, aerosol_masses);
  
  invoke_kernel(aerosol_source_object, std::make_tuple(snl.extent(0)), "ComputeAerosolDeposition");

}

template <typename ArrayB1, typename ArrayI1, typename ArrayD1, typename ArrayD2>
void invoke_aerosol_concen_and_mass(const double& dtime, const ArrayB1 do_capsnow, const ArrayI1 snl,
                                    const ArrayD2 h2osoi_liq, const ArrayD2 h2osoi_ice, const ArrayD2 snw_rds,
                                    const ArrayD1 qflx_snwcp_ice, AerosolMasses<ArrayD2>& aerosol_masses,
                                    AerosolConcentrations<ArrayD2>& aerosol_concentrations)
{
  ComputeAerosolConcenAndMass aerosol_c_mass_object(dtime, do_capsnow, snl, h2osoi_liq, h2osoi_ice, snw_rds,
                                                    qflx_snwcp_ice, aerosol_masses, aerosol_concentrations);

  invoke_kernel(aerosol_c_mass_object, std::make_tuple(snl.extent(0)), "ComputeAerosolConcenAndMass");
}

} // namespace ELM::aerosols
