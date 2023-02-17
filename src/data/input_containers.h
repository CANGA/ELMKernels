
#pragma once

namespace ELM::atm_data {

  // convenience struct
  // holds atmospheric forcing file data interpolated in (time, lat, lon)
  template<typename ArrayD1>
  struct AtmosphereFileInput {
    ArrayD1 atm_tbot;
    ArrayD1 atm_pbot;
    ArrayD1 atm_qbot;
    ArrayD1 atm_flds;
    ArrayD1 atm_fsds;
    ArrayD1 atm_prec;
    ArrayD1 atm_wind;
    ArrayD1 atm_zbot;

    AtmosphereFileInput(int ncols);
    AtmosphereFileInput(const ArrayD1& atm_tbot_in, const ArrayD1& atm_pbot_in,
                        const ArrayD1& atm_qbot_in, const ArrayD1& atm_flds_in,
                        const ArrayD1& atm_fsds_in, const ArrayD1& atm_prec_in,
                        const ArrayD1& atm_wind_in, const ArrayD1& atm_zbot_in);
  };

} // namespace ELM::atm_data

namespace ELM::phen_data {

  // convenience struct
  // holds phenology forcing file data interpolated in (time, lat, lon)
  template<typename ArrayD1>
  struct PhenologyFileInput {
    ArrayD1 lai;
    ArrayD1 sai;
    ArrayD1 htop;
    ArrayD1 hbot;


    PhenologyFileInput(int ncols);
    PhenologyFileInput(const ArrayD1& lai_in, const ArrayD1& sai_in,
                        const ArrayD1& htop_in, const ArrayD1& hbot_in);
  };

  } // namespace ELM::phen_data

#include "input_containers_impl.hh"
