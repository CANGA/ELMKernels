
#pragma once

#include "aerosol_data_new.h"


template<typename ArrayD1>
ELM::aero_data::AerosolFileInput<ArrayD1>::AerosolFileInput(int ncols)
:
  bcphi("bcphi", ncols), bcpho("bcpho", ncols),
  bcdep("bcdep", ncols), dst1_1("dst1_1", ncols),
  dst1_2("dst1_2", ncols), dst2_1("dst2_1", ncols),
  dst2_2("dst2_2", ncols), dst3_1("dst3_1", ncols),
  dst3_2("dst3_2", ncols), dst4_1("dst4_1", ncols),
  dst4_2("dst4_2", ncols)
{}

template<typename ArrayD1>
ELM::aero_data::AerosolFileInput<ArrayD1>::AerosolFileInput(
  const ArrayD1& bcphi_in, const ArrayD1& bcpho_in, const ArrayD1& bcdep_in,
  const ArrayD1& dst1_1_in, const ArrayD1& dst1_2_in, const ArrayD1& dst2_1_in,
  const ArrayD1& dst2_2_in, const ArrayD1& dst3_1_in, const ArrayD1& dst3_2_in,
  const ArrayD1& dst4_1_in, const ArrayD1& dst4_2_in)
:
  bcphi(bcphi_in), bcpho(bcpho_in),
  bcdep(bcdep_in), dst1_1(dst1_1_in),
  dst1_2(dst1_2_in), dst2_1(dst2_1_in),
  dst2_2(dst2_2_in), dst3_1(dst3_1_in),
  dst3_2(dst3_2_in), dst4_1(dst4_1_in),
  dst4_2(dst4_2_in)
{}

template <typename ArrayD2>
ELM::aero_data::AerosolMassess<ArrayD2>::AerosolMassess(int ncols)
:
  mss_bcphi("mss_bcphi", ncols, nlev_sno),
  mss_bcpho("mss_bcpho", ncols, nlev_sno),
  mss_dst1("mss_dst1", ncols, nlev_sno),
  mss_dst2("mss_dst2", ncols, nlev_sno),
  mss_dst3("mss_dst3", ncols, nlev_sno),
  mss_dst4("mss_dst4", ncols, nlev_sno)
{}

template <typename ArrayD2>
ELM::aero_data::AerosolConcentrationss<ArrayD2>::AerosolConcentrationss(int ncols)
:
  cnc_bcphi("cnc_bcphi", ncols, nlev_sno),
  cnc_bcpho("cnc_bcpho", ncols, nlev_sno),
  cnc_dst1("cnc_dst1", ncols, nlev_sno),
  cnc_dst2("cnc_dst2", ncols, nlev_sno),
  cnc_dst3("cnc_dst3", ncols, nlev_sno),
  cnc_dst4("cnc_dst4", ncols, nlev_sno)
{}
