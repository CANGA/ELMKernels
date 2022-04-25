// functions derived from SurfaceRadiationMod.F90

#pragma once

namespace ELM::surface_radiation {

template <class ArrayD1>
ACCELERATED
void initialize_flux(const LandType& Land, double& sabg_soil, double& sabg_snow, double& sabg, double& sabv,
                     double& fsa, ArrayD1 sabg_lyr) {

  // Initialize fluxes
  if (!Land.urbpoi) {
    sabg_soil = 0.0;
    sabg_snow = 0.0;
    sabg = 0.0;
    sabv = 0.0;
    fsa = 0.0;

    for (int j = 0; j < nlevsno + 1; j++) {
      sabg_lyr[j] = 0.0;
    }
  }
}

template <class ArrayD1>
ACCELERATED
void total_absorbed_radiation(const LandType& Land, const int& snl, const ArrayD1 ftdd, const ArrayD1 ftid,
                              const ArrayD1 ftii, const ArrayD1 forc_solad, const ArrayD1 forc_solai,
                              const ArrayD1 fabd, const ArrayD1 fabi, const ArrayD1 albsod, const ArrayD1 albsoi,
                              const ArrayD1 albsnd_hst, const ArrayD1 albsni_hst, const ArrayD1 albgrd,
                              const ArrayD1 albgri, double& sabv, double& fsa, double& sabg, double& sabg_soil,
                              double& sabg_snow, double trd[numrad], double tri[numrad]) {

  double absrad, cad[numrad], cai[numrad];
  if (!Land.urbpoi) {
    for (int ib = 0; ib < numrad; ib++) {

      cad[ib] = forc_solad[ib] * fabd[ib];
      cai[ib] = forc_solai[ib] * fabi[ib];
      sabv += cad[ib] + cai[ib];
      fsa += cad[ib] + cai[ib];

      // Transmitted = solar fluxes incident on ground
      trd[ib] = forc_solad[ib] * ftdd[ib];
      tri[ib] = forc_solad[ib] * ftid[ib] + forc_solai[ib] * ftii[ib];
      // Solar radiation absorbed by ground surface
      // calculate absorbed solar by soil/snow separately
      absrad = trd[ib] * (1.0 - albsod[ib]) + tri[ib] * (1.0 - albsoi[ib]);
      sabg_soil += absrad;
      absrad = trd[ib] * (1.0 - albsnd_hst[ib]) + tri[ib] * (1.0 - albsni_hst[ib]);
      sabg_snow += absrad;
      absrad = trd[ib] * (1.0 - albgrd[ib]) + tri[ib] * (1.0 - albgri[ib]);
      sabg += absrad;
      fsa += absrad;

      if (snl == 0) {
        sabg_snow = sabg;
        sabg_soil = sabg;
      }

      if (subgridflag == 0) {
        sabg_snow = sabg;
        sabg_soil = sabg;
      }
    } // end of numrad
  }   // end if not urban
}

template <class ArrayD1>
ACCELERATED
void layer_absorbed_radiation(const LandType& Land, const int& snl, const double& sabg, const double& sabg_snow,
                              const double& snow_depth, const ArrayD1 flx_absdv, const ArrayD1 flx_absdn,
                              const ArrayD1 flx_absiv, const ArrayD1 flx_absin, const double trd[numrad],
                              const double tri[numrad], ArrayD1 sabg_lyr) {

  double err_sum = 0.0;
  double sabg_snl_sum;
  // compute absorbed flux in each snow layer and top soil layer,
  // based on flux factors computed in the radiative transfer portion of SNICAR.
  if (!Land.urbpoi) {

    sabg_snl_sum = 0.0;

    // CASE1: No snow layers: all energy is absorbed in top soil layer
    if (snl == 0) {
      for (int i = 0; i <= nlevsno; i++) {
        sabg_lyr[i] = 0.0;
      }
      sabg_lyr[nlevsno] = sabg;
      sabg_snl_sum = sabg_lyr[nlevsno];
    } else { // CASE 2: Snow layers present: absorbed radiation is scaled according to flux factors computed by SNICAR

      for (int i = 0; i < nlevsno + 1; i++) {
        sabg_lyr[i] = flx_absdv[i] * trd[0] + flx_absdn[i] * trd[1] + flx_absiv[i] * tri[0] + flx_absin[i] * tri[1];
        // summed radiation in active snow layers:
        // if snow layer is at or below snow surface
        if (i >= nlevsno - snl) {
          sabg_snl_sum += sabg_lyr[i];
        }
      }

      // Error handling: The situation below can occur when solar radiation is
      // NOT computed every timestep.
      // When the number of snow layers has changed in between computations of the
      // absorbed solar energy in each layer, we must redistribute the absorbed energy
      // to avoid physically unrealistic conditions. The assumptions made below are
      // somewhat arbitrary, but this situation does not arise very frequently.
      // This error handling is implemented to accomodate any value of the
      // radiation frequency.
      // change condition to match sabg_snow isntead of sabg

      if (std::abs(sabg_snl_sum - sabg_snow) > 0.00001) {
        if (snl == 0) {
          for (int j = 0; j < nlevsno; j++) {
            sabg_lyr[j] = 0.0;
          }
          sabg_lyr[nlevsno] = sabg;
        } else if (snl == 1) {
          for (int j = 0; j < nlevsno - 1; j++) {
            sabg_lyr[j] = 0.0;
          }
          sabg_lyr[nlevsno - 1] = sabg_snow * 0.6;
          sabg_lyr[nlevsno] = sabg_snow * 0.4;
        } else {
          for (int j = 0; j <= nlevsno; j++) {
            sabg_lyr[j] = 0.0;
          }
          sabg_lyr[nlevsno - snl] = sabg_snow * 0.75;
          sabg_lyr[nlevsno - snl + 1] = sabg_snow * 0.25;
        }
      }

      // If shallow snow depth, all solar radiation absorbed in top or top two snow layers
      // to prevent unrealistic timestep soil warming

      if (subgridflag == 0) {
        if (snow_depth < 0.1) {
          if (snl == 0) {
            for (int j = 0; j < nlevsno; j++) {
              sabg_lyr[j] = 0.0;
            }
            sabg_lyr[nlevsno] = sabg;
          } else if (snl == 1) {
            for (int j = 0; j < nlevsno - 1; j++) {
              sabg_lyr[j] = 0.0;
            }
            sabg_lyr[nlevsno - 1] = sabg;
            sabg_lyr[nlevsno] = 0.0;
          } else {
            for (int j = 0; j <= nlevsno; j++) {
              sabg_lyr[j] = 0.0;
            }
            sabg_lyr[nlevsno - snl] = sabg_snow * 0.75;
            sabg_lyr[nlevsno - snl + 1] = sabg_snow * 0.25;
          }
        }
      }
    }

    // Error check - This situation should not happen:
    for (int j = 0; j <= nlevsno; j++) {
      err_sum += sabg_lyr[j];
    }
    assert(!(std::abs(err_sum - sabg_snow) > 0.00001));
  }
}

template <class ArrayD1>
ACCELERATED
void reflected_radiation(const LandType& Land, const ArrayD1 albd, const ArrayD1 albi, const ArrayD1 forc_solad,
                         const ArrayD1 forc_solai, double& fsr) {

  double fsr_vis_d, fsr_nir_d, fsr_vis_i, fsr_nir_i, rvis, rnir;
  // Radiation diagnostics
  if (!Land.urbpoi) {
    // NDVI and reflected solar radiation
    rvis = albd[0] * forc_solad[0] + albi[0] * forc_solai[0];
    rnir = albd[1] * forc_solad[1] + albi[1] * forc_solai[1];
    fsr = rvis + rnir;
  } else {
    // Solar reflected per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
    fsr_vis_d = albd[0] * forc_solad[0];
    fsr_nir_d = albd[1] * forc_solad[1];
    fsr_vis_i = albi[0] * forc_solai[0];
    fsr_nir_i = albi[1] * forc_solai[1];

    fsr = fsr_vis_d + fsr_nir_d + fsr_vis_i + fsr_nir_i;
  }
}

template <class ArrayD1>
ACCELERATED
void canopy_sunshade_fractions(const LandType& Land, const int& nrad, const double& elai, const ArrayD1 tlai_z,
                               const ArrayD1 fsun_z, const ArrayD1 forc_solad, const ArrayD1 forc_solai,
                               const ArrayD1 fabd_sun_z, const ArrayD1 fabd_sha_z, const ArrayD1 fabi_sun_z,
                               const ArrayD1 fabi_sha_z, ArrayD1 parsun_z, ArrayD1 parsha_z, ArrayD1 laisun_z,
                               ArrayD1 laisha_z, double& laisun, double& laisha) {

  if (!Land.urbpoi) {
    int ipar = 0; // The band index for PAR
    for (int iv = 0; iv < nrad; iv++) {
      parsun_z[iv] = 0.0;
      parsha_z[iv] = 0.0;
      laisun_z[iv] = 0.0;
      laisha_z[iv] = 0.0;
    }
    // Loop over patches to calculate laisun_z and laisha_z for each layer.
    // Derive canopy laisun, laisha, from layer sums.
    // If sun/shade big leaf code, nrad=1 and fsun_z[0] and tlai_z[0] from
    // SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
    laisun = 0.0;
    laisha = 0.0;

    for (int iv = 0; iv < nrad; iv++) {
      laisun_z[iv] = tlai_z[iv] * fsun_z[iv];
      laisha_z[iv] = tlai_z[iv] * (1.0 - fsun_z[iv]);
      laisun += laisun_z[iv];
      laisha += laisha_z[iv];
    }

    // Absorbed PAR profile through canopy
    // If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
    // are canopy integrated so that layer values equal big leaf values.
    for (int iv = 0; iv < nrad; iv++) {
      parsun_z[iv] = forc_solad[ipar] * fabd_sun_z[iv] + forc_solai[ipar] * fabi_sun_z[iv];
      parsha_z[iv] = forc_solad[ipar] * fabd_sha_z[iv] + forc_solai[ipar] * fabi_sha_z[iv];
    }
  }
}

} // namespace ELM::surface_radiation
