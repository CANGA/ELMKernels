
#pragma once

#include "elm_constants.h"

namespace ELM::aero_data {

  // convenience struct
  // holds aerosol forcing file data streams interpolated in (time, lat, lon) 
  template<typename ArrayD1>
  struct AerosolFileInput {
    ArrayD1 bcphi;
    ArrayD1 bcpho;
    ArrayD1 bcdep;
    ArrayD1 dst1_1;
    ArrayD1 dst1_2;
    ArrayD1 dst2_1;
    ArrayD1 dst2_2;
    ArrayD1 dst3_1;
    ArrayD1 dst3_2;
    ArrayD1 dst4_1;
    ArrayD1 dst4_2;
    AerosolFileInput(int ncols);
    AerosolFileInput(const ArrayD1& bcphi_in, const ArrayD1& bcpho_in, const ArrayD1& bcdep_in,
                     const ArrayD1& dst1_1_in, const ArrayD1& dst1_2_in, const ArrayD1& dst2_1_in,
                     const ArrayD1& dst2_2_in, const ArrayD1& dst3_1_in, const ArrayD1& dst3_2_in,
                     const ArrayD1& dst4_1_in, const ArrayD1& dst4_2_in);
  };


  template <typename ArrayD2>
  struct AerosolMasses {
    ArrayD2 mss_bcphi;
    ArrayD2 mss_bcpho;
    ArrayD2 mss_dst1;
    ArrayD2 mss_dst2;
    ArrayD2 mss_dst3;
    ArrayD2 mss_dst4;
    AerosolMasses(int ncols);
  };

  template <typename ArrayD2>
  struct AerosolConcentrations {
    ArrayD2 cnc_bcphi;
    ArrayD2 cnc_bcpho;
    ArrayD2 cnc_dst1;
    ArrayD2 cnc_dst2;
    ArrayD2 cnc_dst3;
    ArrayD2 cnc_dst4;
    AerosolConcentrations(int ncols);
  };

} // namespace ELM::aero_data

#include "aerosol_data_impl.hh"
