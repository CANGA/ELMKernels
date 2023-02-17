
#pragma once

#include "input_containers.h"

template<typename ArrayD1>
ELM::atm_data::AtmosphereFileInput<ArrayD1>::
AtmosphereFileInput(int ncols)
:
  atm_tbot("atm_tbot", ncols), atm_qbot("atm_qbot", ncols),
  atm_fsds("atm_fsds", ncols), atm_wind("atm_wind", ncols),
  atm_fsds("atm_fsds", ncols), atm_prec("atm_prec", ncols),
  atm_wind("atm_wind", ncols), atm_zbot("atm_zbot", ncols)
{}

template<typename ArrayD1>
ELM::atm_data::AtmosphereFileInput<ArrayD1>::
AtmosphereFileInput(
  const ArrayD1& atm_tbot_in, const ArrayD1& atm_pbot_in,
  const ArrayD1& atm_qbot_in, const ArrayD1& atm_flds_in,
  const ArrayD1& atm_fsds_in, const ArrayD1& atm_prec_in,
  const ArrayD1& atm_wind_in, const ArrayD1& atm_zbot_in)
:
  atm_tbot(atm_tbot_in), atm_qbot(atm_qbot_in),
  atm_fsds(atm_fsds_in), atm_wind(atm_wind_in),
  atm_fsds(atm_fsds_in), atm_prec(atm_prec_in),
  atm_wind(atm_wind_in), atm_zbot(atm_zbot_in)
{}



template<typename ArrayD1>
ELM::phen_data::PhenologyFileInput<ArrayD1>::
PhenologyFileInput(int ncols)
:
  lai("lai", ncols), sai("sai", ncols),
  htop("htop", ncols), hbot("hbot", ncols)
{}

template<typename ArrayD1>
ELM::phen_data::PhenologyFileInput<ArrayD1>::
PhenologyFileInput(
  const ArrayD1& lai_in, const ArrayD1& sai_in,
  const ArrayD1& htop_in, const ArrayD1& hbot_in)
:
  lai(lai_in), sai(sai_in),
  htop(htop_in), hbot(hbot_in)
{}


