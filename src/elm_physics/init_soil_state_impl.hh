
#pragma once

namespace ELM {

// from ColumnDataType.F90
//-----------------------------------------------------------------------
// set cold-start initial values for select members of col_es
//-----------------------------------------------------------------------
template <typename ArrayD1>
ACCELERATE
void init_soil_temp(const LandType& Land, const int& snl, ArrayD1 t_soisno, double& t_grnd)
{
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;
  using ELMdims::nlevurb;

  // Snow level temperatures - all land points
  if (snl > 0) {
    for (int i = nlevsno - snl; i < nlevsno; ++i) {
      t_soisno(i) = 250.0;
    }
  }

  // Below snow temperatures - nonlake points (lake points are set below)
  if (!Land.lakpoi) {
    if (Land.ltype == LND::istice || Land.ltype == LND::istice_mec) {
      for (int i = nlevsno; i < nlevgrnd + nlevsno; ++i) {
        t_soisno(i) = 250.0;
      }
    } else if (Land.ltype == LND::istwet) {
      for (int i = nlevsno; i < nlevgrnd + nlevsno; ++i) {
        t_soisno(i) = 277.0;
      }
    } else if (Land.urbpoi) {
      if (Land.ctype == LND::icol_road_perv || Land.ctype == LND::icol_road_imperv) {
        for (int i = nlevsno; i < nlevgrnd + nlevsno; ++i) {
          t_soisno(i) = 274.0;
        }
      } else if (Land.ctype == LND::icol_sunwall || Land.ctype == LND::icol_shadewall || Land.ctype == LND::icol_roof) {
        // Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
        // shock from large heating/air conditioning flux
        for (int i = nlevsno; i < nlevurb + nlevsno; ++i) {
          t_soisno(i) = 292.0;
        }
      }
    } else {
      for (int i = nlevsno; i < nlevgrnd + nlevsno; ++i) {
        t_soisno(i) = 274.0;
      }
    }
    t_grnd = t_soisno(nlevsno - snl);
  }
}

// from ColumnDataType.F90 and WaterStateType.F90
//--------------------------------------------
// Set soil water
//--------------------------------------------
// volumetric water is set first and liquid content and ice lens are obtained
// NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
// and urban pervious road (other urban columns have zero soil water)
template <typename ArrayD1>
ACCELERATE
void init_soilh2o_state(const LandType& Land, const int& snl, const ArrayD1 watsat, const ArrayD1 t_soisno,
                        const ArrayD1 dz, ArrayD1 h2osoi_vol, ArrayD1 h2osoi_liq, ArrayD1 h2osoi_ice)
{
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;
  using ELMdims::nlevsoi;
  using ELMdims::nlevbed;
  using ELMdims::nlevurb;
  using ELMconst::SPVAL;

  for (int i = 0; i < nlevgrnd; ++i) {
    h2osoi_vol(i) = SPVAL;
  }
  for (int i = 0; i < nlevgrnd + nlevsno; ++i) {
    h2osoi_liq(i) = SPVAL;
  }
  for (int i = 0; i < nlevgrnd + nlevsno; ++i) {
    h2osoi_ice(i) = SPVAL;
  }

  int nlevs = nlevgrnd;
  if (!Land.lakpoi) {

    if (Land.ltype == LND::istsoil || Land.ltype == LND::istcrop) {
      for (int i = 0; i < nlevgrnd; ++i) {
        if (i >= nlevbed) {
          h2osoi_vol(i) = 0.0;
        } else {
          // if (use_fates_planthydro .or. use_hydrstress) {
          //   h2osoi_vol(i) = 0.70_r8*watsat(i) !0.15_r8 to avoid very dry conditions that cause errors in FATES HYDRO
          // } else {
          h2osoi_vol(i) = 0.15;
        }
      }
    } else if (Land.urbpoi) {
      if (Land.ctype == LND::icol_road_perv) {
        for (int i = 0; i < nlevgrnd; ++i) {
          if (i < nlevbed) {
            h2osoi_vol(i) = 0.3;
          } else {
            h2osoi_vol(i) = 0.0;
          }
        }
      } else if (Land.ctype == LND::icol_road_imperv) {
        for (int i = 0; i < nlevgrnd; ++i) {
          h2osoi_vol(i) = 0.0;
        }
      } else {
       nlevs = nlevurb;
        for (int i = 0; i < nlevurb; ++i) {
          h2osoi_vol(i) = 0.0;
        }
      }
    } else if (Land.ltype == LND::istwet) {
      for (int i = 0; i < nlevgrnd; ++i) {
        if (i >= nlevbed) {
          h2osoi_vol(i) = 0.0;
        } else {
          h2osoi_vol(i) = 1.0;
        }
      }
    } else if (Land.ltype == LND::istice || Land.ltype == LND::istice_mec) {
      for (int i = 0; i < nlevgrnd; ++i) {
        h2osoi_vol(i) = 1.0;
      }
    }
    for (int i = 0; i < nlevs; ++i) {
      const int snw_offset{i + nlevsno};
      h2osoi_vol(i) = std::min(h2osoi_vol(i), watsat(i));
      if (t_soisno(snw_offset) <= ELMconst::TFRZ) {
        h2osoi_ice(snw_offset) = dz(snw_offset) * ELMconst::DENICE * h2osoi_vol(i);
        h2osoi_liq(snw_offset) = 0.0;
      } else {
        h2osoi_ice(snw_offset) = 0.0;
        h2osoi_liq(snw_offset) = dz(snw_offset) * ELMconst::DENH2O * h2osoi_vol(i);
      }
    }

    for (int i = 0; i < nlevsno; ++i) {
      if (i >= nlevsno - snl) {
        h2osoi_ice(i) = dz(i) * 250.0;
        h2osoi_liq(i) = 0.0;
      }
    }
  } else {
    //--------------------------------------------
    // Set Lake water
    //--------------------------------------------
    for (int i = 0; i < nlevsno; ++i) {
      if (i >= nlevsno - snl) {
        h2osoi_ice(i) = dz(i) * ELMconst::BDSNO;
        h2osoi_liq(i) = 0.0;
      }
    }
    for (int i = 0; i < nlevgrnd; ++i) {
      const int snw_offset{i + nlevsno};
      if (i < nlevsoi) { // soil
        h2osoi_vol(i) = watsat(i);
        h2osoi_liq(snw_offset) = SPVAL;
        h2osoi_ice(snw_offset) = SPVAL;
      } else { // bedrock
        h2osoi_vol(i) = 0.0;
      }
    }
  }

  //--------------------------------------------
  // For frozen layers !TODO - does the following make sense ???? it seems to overwrite everything
  //--------------------------------------------
  for (int i = 0; i < nlevgrnd; ++i) {
    const int snw_offset{i + nlevsno};
    if (t_soisno(snw_offset) <= ELMconst::TFRZ) {
      h2osoi_ice(snw_offset) = dz(snw_offset) * ELMconst::DENICE * h2osoi_vol(i);
      h2osoi_liq(snw_offset) = 0.0;
    } else {
      h2osoi_ice(snw_offset) = 0.0;
      h2osoi_liq(snw_offset) = dz(snw_offset) * ELMconst::DENH2O * h2osoi_vol(i);
    }
  }
  // h2osoi_liq_old(c,:) = h2osoi_liq(c,:)
  // h2osoi_ice_old(c,:) = h2osoi_ice(c,:)
}

template <typename ArrayD1>
ACCELERATE
void init_vegrootfr(const int& vtype, const double& roota_par, const double& rootb_par,
                    const ArrayD1 zi,ArrayD1 rootfr)
{
  using ELMdims::nlevsoi;
  using ELMdims::nlevgrnd;
  using ELMdims::nlevsno;
  using ELMdims::nlevbed;
  // (computing from surface, d is depth in meters: Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
  // Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs given in Zeng et al. (1998).

  for (int i = nlevsoi; i < nlevgrnd; ++i) {
    rootfr(i) = 0.0;
  }

  if (vtype != PFT::noveg) {
    double totrootfr = 0.0;
    for (int i = 0; i <nlevsoi - 1; i++) {
      rootfr(i) = 0.5 * (exp(-roota_par * zi(i +nlevsno)) + exp(-rootb_par * zi(i +nlevsno)) -
                         exp(-roota_par * zi(i + 1 +nlevsno)) - exp(-rootb_par * zi(i + 1 +nlevsno)));
      if (i <nlevbed) {
        totrootfr += rootfr(i);
      }
    }
    rootfr(nlevsoi - 1) =
        0.5 * (exp(-roota_par * zi(nlevsoi - 1 +nlevsno)) + exp(-rootb_par * zi(nlevsoi - 1 +nlevsno)));
  } else {
    for (int i = 0; i <nlevsoi; i++) {
      rootfr(i) = 0.0;
    }
  }
  for (int i = nlevsoi; i <nlevgrnd; i++) {
    rootfr(i) = 0.0;
  }
}

} // namespace ELM
