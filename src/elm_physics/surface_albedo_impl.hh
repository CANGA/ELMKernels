// Functions derived from SurfaceAlbedoMod.F90

// call sequence -> SurfAlbInitTimestep() -> soil_albedo() -> SNICAR_AD_RT() (once for both wavebands) ->
// ground_albedo() -> flux_absorption_factor() CanopyLayers() -> two_stream_solver()

// will need to finish zenith angle ATS code, Aerosol functions

// need to initialize h2osoi_vol[nlevgrnd]
// need to read in soil color NetCDF, initialize isoicol, albsat, albdry
// need to figure out doalb - during nstep == 0 SurfaceAlbedo doesn't get called and values for the outputs are provided
// by SurfaceAlbedoType::InitCold
//  why does SurfAlb calc at nstep+1? there must be a reason
//  from CLM 4.5 tech note:
/*

The land model calculations are implemented in two steps. The land
model proceeds with the calculation of surface energy, constituent, momentum, and
radiative fluxes using the snow and soil hydrologic states from the previous time step.
The land model then updates the soil and snow hydrology calculations based on these
fluxes. These fields are passed to the atmosphere (Table 2.4). The albedos sent to the
atmosphere are for the solar zenith angle at the next time step but with surface conditions
from the current time step.

So it's for the ATM coupling
so we can probably calc at beginning of nstep with current solar zenith angle


these don't need to be persistent at the driver level, but will need to be passed from SNICAR_AD_RT() to
flux_absorption_factor() flx_absd_snw flx_absi_snw mss_cnc_aer_in_fdb - from SurfAlbInitTimestep() to SNICAR_AD_RT()

*/

#pragma once

namespace ELM::surface_albedo {

// not sure if necessary
// maybe for nstep == 0 this should be run, don't call surfalb, but call evrything else?
//  will need to investigate initialization of coupled system - ELM calls dummy nstep == 0 for a reason
// void SurfaceAlbedo_InitCold() {
// albgrd[numrad] = {0.2};
// albgri[numrad] = {0.2};
// albsod[numrad] = {0.2};
// albsoi[numrad] = {0.2};
// albd[numrad] = {0.2};
// albi[numrad] = {0.2};
// fabi[numrad] = {0.0};
// fabd[numrad] = {0.0};
// fabi_sun[numrad] = {0.0};
// fabd_sun[numrad] = {0.0};
// fabd_sha[numrad] = {0.0};
// fabi_sha[numrad] = {0.0};
// ftdd[numrad] = {1.0};
// ftid[numrad] = {0.0};
// ftii[numrad] = {1.0};
//}


// returns the cosine of the solar zenith angle.
// not currently used, i don't think
ACCELERATE
double calc_cosz(const double& jday, const double& lat, const double& lon, const double& declin)
{
  return sin(lat) * sin(declin) - cos(lat) * cos(declin) * cos((jday - floor(jday)) * 2.0 * ELM::constants::ELM_PI + lon);
}

ACCELERATE
bool vegsol(const LandType& Land, const double& coszen, const double& elai, const double& esai)
{
  if (!Land.urbpoi && coszen > 0.0 && (Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0) {
    return true;
  } else {
    return false;
  }
}

ACCELERATE
bool novegsol(const LandType& Land, const double& coszen, const double& elai, const double& esai)
{
  if (!Land.urbpoi && coszen > 0.0) {
    if (!((Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0)) {
      return true;
    }
  }
  return false;
}

template <class ArrayD1, class ArrayD2>
ACCELERATE
void init_timestep(const bool& urbpoi, const double& elai, const ArrayD1 mss_cnc_bcphi, const ArrayD1 mss_cnc_bcpho,
                   const ArrayD1 mss_cnc_dst1, const ArrayD1 mss_cnc_dst2, const ArrayD1 mss_cnc_dst3,
                   const ArrayD1 mss_cnc_dst4, double& vcmaxcintsun, double& vcmaxcintsha, ArrayD1 albsod,
                   ArrayD1 albsoi, ArrayD1 albgrd, ArrayD1 albgri, ArrayD1 albd, ArrayD1 albi, ArrayD1 fabd,
                   ArrayD1 fabd_sun, ArrayD1 fabd_sha, ArrayD1 fabi, ArrayD1 fabi_sun, ArrayD1 fabi_sha, ArrayD1 ftdd,
                   ArrayD1 ftid, ArrayD1 ftii, ArrayD1 flx_absdv, ArrayD1 flx_absdn, ArrayD1 flx_absiv,
                   ArrayD1 flx_absin, ArrayD2 mss_cnc_aer_in_fdb)
{
  // Initialize output because solar radiation only done if coszen > 0
  if (!urbpoi) {
    for (int ib = 0; ib < numrad; ++ib) {
      albsod(ib) = 0.0;
      albsoi(ib) = 0.0;
      albgrd(ib) = 0.0;
      albgri(ib) = 0.0;
      albd(ib) = 1.0;
      albi(ib) = 1.0;
      fabd(ib) = 0.0;
      fabd_sun(ib) = 0.0;
      fabd_sha(ib) = 0.0;
      fabi(ib) = 0.0;
      fabi_sun(ib) = 0.0;
      fabi_sha(ib) = 0.0;
      ftdd(ib) = 0.0;
      ftid(ib) = 0.0;
      ftii(ib) = 0.0;
    }
    for (int i = 0; i <= nlevsno; ++i) {
      flx_absdv(i) = 0.0;
      flx_absdn(i) = 0.0;
      flx_absiv(i) = 0.0;
      flx_absin(i) = 0.0;
    }
    if (nlevcan == 1) {
      vcmaxcintsun = 0.0;
      vcmaxcintsha = (1.0 - exp(-extkn * elai)) / extkn;
      if (elai > 0.0) {
        vcmaxcintsha /= elai;
      } else {
        vcmaxcintsha = 0.0;
      }
    } else if (nlevcan > 1) {
      vcmaxcintsun = 0.0;
      vcmaxcintsha = 0.0;
    }
  }
  // set soot and dust aerosol concentrations:
  for (int i = 0; i < nlevsno; ++i) {
    mss_cnc_aer_in_fdb(i, 0) = mss_cnc_bcphi(i);
    mss_cnc_aer_in_fdb(i, 1) = mss_cnc_bcpho(i);
    mss_cnc_aer_in_fdb(i, 2) = 0.0; // ignore OC concentrations due to poor constraint and negligible effect on snow
    mss_cnc_aer_in_fdb(i, 3) = 0.0; // ignore OC concentrations due to poor constraint and negligible effect on snow
    mss_cnc_aer_in_fdb(i, 4) = mss_cnc_dst1(i);
    mss_cnc_aer_in_fdb(i, 5) = mss_cnc_dst2(i);
    mss_cnc_aer_in_fdb(i, 6) = mss_cnc_dst3(i);
    mss_cnc_aer_in_fdb(i, 7) = mss_cnc_dst4(i);
  }
} // init_timestep

template <class ArrayD1>
ACCELERATE
void ground_albedo(const bool& urbpoi, const double& coszen, const double& frac_sno, const ArrayD1 albsod,
                   const ArrayD1 albsoi, const ArrayD1 albsnd, const ArrayD1 albsni, ArrayD1 albgrd, ArrayD1 albgri)
{
  if (!urbpoi && coszen > 0.0) {
    for (int ib = 0; ib < numrad; ++ib) {
      albgrd(ib) = albsod(ib) * (1.0 - frac_sno) + albsnd(ib) * frac_sno;
      albgri(ib) = albsoi(ib) * (1.0 - frac_sno) + albsni(ib) * frac_sno;
    }
  }
} // ground_albedo

template <class ArrayD1, class ArrayD2>
ACCELERATE
void flux_absorption_factor(const LandType& Land, const double& coszen, const double& frac_sno, const ArrayD1 albsod,
                            const ArrayD1 albsoi, const ArrayD1 albsnd, const ArrayD1 albsni,
                            const ArrayD2 flx_absd_snw, const ArrayD2 flx_absi_snw, ArrayD1 flx_absdv,
                            ArrayD1 flx_absdn, ArrayD1 flx_absiv, ArrayD1 flx_absin)
{
  // ground albedos and snow-fraction weighting of snow absorption factors
  if (!Land.urbpoi && coszen > 0.0) {
    for (int i = 0; i <= nlevsno; ++i) {
      for (int ib = 0; ib < numrad; ++ib) {
        // also in this loop (but optionally in a different loop for vectorized code)
        // weight snow layer radiative absorption factors based on snow fraction and soil albedo
        // (NEEDED FOR ENERGY CONSERVATION)
        if (subgridflag == 0 || Land.ltype == istdlak) {
          if (ib == 0) {
            flx_absdv(i) = flx_absd_snw(i, ib) * frac_sno +
                           ((1.0 - frac_sno) * (1.0 - albsod(ib)) * (flx_absd_snw(i, ib) / (1.0 - albsnd(ib))));
            flx_absiv(i) = flx_absi_snw(i, ib) * frac_sno +
                           ((1.0 - frac_sno) * (1.0 - albsoi(ib)) * (flx_absi_snw(i, ib) / (1.0 - albsni(ib))));
          } else if (ib == 1) {
            flx_absdn(i) = flx_absd_snw(i, ib) * frac_sno +
                           ((1.0 - frac_sno) * (1.0 - albsod(ib)) * (flx_absd_snw(i, ib) / (1.0 - albsnd(ib))));
            flx_absin(i) = flx_absi_snw(i, ib) * frac_sno +
                           ((1.0 - frac_sno) * (1.0 - albsoi(ib)) * (flx_absi_snw(i, ib) / (1.0 - albsni(ib))));
          }
        } else {
          if (ib == 0) {
            flx_absdv(i) = flx_absd_snw(i, ib) * (1.0 - albsnd(ib));
            flx_absiv(i) = flx_absi_snw(i, ib) * (1.0 - albsni(ib));
          } else if (ib == 1) {
            flx_absdn(i) = flx_absd_snw(i, ib) * (1.0 - albsnd(ib));
            flx_absin(i) = flx_absi_snw(i, ib) * (1.0 - albsni(ib));
          }
        } // if subgridflag
      }   // for numrad
    }     // for nlevsno +1 layers
  }       // if !Land.urbpoi && coszen > 0.0
} // flux_absorption_factor

template <class ArrayD1>
ACCELERATE
void canopy_layer_lai(const int& urbpoi, const double& elai, const double& esai, const double& tlai, const double& tsai,
                      int& nrad, int& ncan, ArrayD1 tlai_z, ArrayD1 tsai_z, ArrayD1 fsun_z, ArrayD1 fabd_sun_z,
                      ArrayD1 fabd_sha_z, ArrayD1 fabi_sun_z, ArrayD1 fabi_sha_z)
{
  if (!urbpoi) {
    if (nlevcan == 1) {
      nrad = 1;
      ncan = 1;
      tlai_z(0) = elai;
      tsai_z(0) = esai;
    } else if (nlevcan > 1) {
      if (elai + esai == 0.0) {
        nrad = 0;
      } else {
        double dincmax_sum = 0.0;
        for (int iv = 0; iv < nlevcan; ++iv) {
          dincmax_sum += dincmax;
          if (((elai + esai) - dincmax_sum) > mpe) {
            nrad = iv + 1;
            double dinc = dincmax;
            tlai_z(iv) = dinc * elai / std::max(elai + esai, mpe);
            tsai_z(iv) = dinc * esai / std::max(elai + esai, mpe);
          } else {
            nrad = iv + 1;
            double dinc = dincmax - (dincmax_sum - (elai + esai));
            tlai_z(iv) = dinc * elai / std::max(elai + esai, mpe);
            tsai_z(iv) = dinc * esai / std::max(elai + esai, mpe);
            break;
          }
        }
        // Mimumum of 4 canopy layers
        if (nrad < 4) {
          nrad = 4;
          for (int iv = 0; iv < nrad; ++iv) {
            tlai_z(iv) = elai / nrad;
            tsai_z(iv) = esai / nrad;
          }
        }
      }
    }

    // Error check: make sure cumulative of increments does not exceed total
    double laisum = 0.0;
    double saisum = 0.0;
    for (int iv = 0; iv < nrad; ++iv) {
      laisum += tlai_z(iv);
      saisum += tsai_z(iv);
    }

    if (std::abs(laisum - elai) > mpe || std::abs(saisum - esai) > mpe) {
      throw std::runtime_error("ELM ERROR: multi-layer canopy error 1 in SurfaceAlbedo");
    }

    // Repeat to find canopy layers buried by snow
    if (nlevcan > 1) {
      double blai = tlai - elai;
      double bsai = tsai - esai;
      if (blai + bsai == 0.0) {
        ncan = nrad;
      } else {
        double dincmax_sum = 0.0;
        for (int iv = nrad; iv < nlevcan; ++iv) {
          dincmax_sum += dincmax;
          if (((blai + bsai) - dincmax_sum) > mpe) {
            ncan = iv + 1;
            double dinc = dincmax;
            tlai_z(iv) = dinc * blai / std::max(blai + bsai, mpe);
            tsai_z(iv) = dinc * bsai / std::max(blai + bsai, mpe);
          } else {
            ncan = iv + 1;
            double dinc = dincmax - (dincmax_sum - (blai + bsai));
            tlai_z(iv) = dinc * blai / std::max(blai + bsai, mpe);
            tsai_z(iv) = dinc * bsai / std::max(blai + bsai, mpe);
            break;
          }
        }
      }

      // Error check: make sure cumulative of increments does not exceed total
      laisum = 0.0;
      saisum = 0.0;
      for (int iv = 0; iv < ncan; ++iv) {
        laisum += tlai_z(iv);
        saisum += tsai_z(iv);
      }
      if (std::abs(laisum - tlai) > mpe || std::abs(saisum - tsai) > mpe) {
        throw std::runtime_error("ELM ERROR: multi-layer canopy error 1 in SurfaceAlbedo");
      }
    }

    // Zero fluxes for active canopy layers
    for (int iv = 0; iv < nrad; ++iv) {
      fabd_sun_z(iv) = 0.0;
      fabd_sha_z(iv) = 0.0;
      fabi_sun_z(iv) = 0.0;
      fabi_sha_z(iv) = 0.0;
      fsun_z(iv) = 0.0;
    }
  } // if !urbpoi
} // canopy_layer_lai

template <class ArrayD1>
ACCELERATE
void two_stream_solver(const LandType& Land, const int& nrad, const double& coszen, const double& t_veg,
                       const double& fwet, const double& elai, const double& esai, const ArrayD1 tlai_z,
                       const ArrayD1 tsai_z, const ArrayD1 albgrd, const ArrayD1 albgri, const PFTDataAlb& alb_pft,
                       double& vcmaxcintsun, double& vcmaxcintsha, ArrayD1 albd, ArrayD1 ftid, ArrayD1 ftdd,
                       ArrayD1 fabd, ArrayD1 fabd_sun, ArrayD1 fabd_sha, ArrayD1 albi, ArrayD1 ftii, ArrayD1 fabi,
                       ArrayD1 fabi_sun, ArrayD1 fabi_sha, ArrayD1 fsun_z, ArrayD1 fabd_sun_z, ArrayD1 fabd_sha_z,
                       ArrayD1 fabi_sun_z, ArrayD1 fabi_sha_z)
{
  if (vegsol(Land, coszen, elai, esai)) {

    double omega[numrad]; // fraction of intercepted radiation that is scattered (0 to 1)
    double rho[numrad], tau[numrad];

    // Weight reflectance/transmittance by lai and sai
    const double wl = elai / std::max(elai + esai, mpe);
    const double ws = esai / std::max(elai + esai, mpe);

    // Calculate two-stream parameters that are independent of waveband:
    // chil, gdir, twostext, avmu, and temp0 and temp2 (used for asu)
    const double cosz = std::max(0.001, coszen);
    double chil = std::min(std::max(alb_pft.xl, -0.4), 0.6);
    if (std::abs(chil) <= 0.01) {
      chil = 0.01;
    }

    const double phi1 = 0.5 - 0.633 * chil - 0.330 * chil * chil;
    const double phi2 = 0.877 * (1.0 - 2.0 * phi1);
    const double gdir = phi1 + phi2 * cosz;
    const double twostext = gdir / cosz;
    const double avmu = (1.0 - phi1 / phi2 * std::log((phi1 + phi2) / phi1)) / phi2;
    const double temp0 = gdir + phi2 * cosz;
    const double temp1 = phi1 * cosz;
    const double temp2 = (1.0 - temp1 / temp0 * std::log((temp1 + temp0) / temp1));

    // Loop over all wavebands to calculate for the full canopy the scattered fluxes
    // reflected upward and transmitted downward by the canopy and the flux absorbed by the
    // canopy for a unit incoming direct beam and diffuse flux at the top of the canopy given
    // an underlying surface of known albedo.
    // Output:
    // ------------------
    // Direct beam fluxes
    // ------------------
    // albd       - Upward scattered flux above canopy (per unit direct beam flux)
    // ftid       - Downward scattered flux below canopy (per unit direct beam flux)
    // ftdd       - Transmitted direct beam flux below canopy (per unit direct beam flux)
    // fabd       - Flux absorbed by canopy (per unit direct beam flux)
    // fabd_sun   - Sunlit portion of fabd
    // fabd_sha   - Shaded portion of fabd
    // fabd_sun_z - absorbed sunlit leaf direct PAR (per unit sunlit lai+sai) for each canopy layer
    // fabd_sha_z - absorbed shaded leaf direct PAR (per unit shaded lai+sai) for each canopy layer
    // ------------------
    // Diffuse fluxes
    // ------------------
    // albi       - Upward scattered flux above canopy (per unit diffuse flux)
    // ftii       - Downward scattered flux below canopy (per unit diffuse flux)
    // fabi       - Flux absorbed by canopy (per unit diffuse flux)
    // fabi_sun   - Sunlit portion of fabi
    // fabi_sha   - Shaded portion of fabi
    // fabi_sun_z - absorbed sunlit leaf diffuse PAR (per unit sunlit lai+sai) for each canopy layer
    // fabi_sha_z - absorbed shaded leaf diffuse PAR (per unit shaded lai+sai) for each canopy layer

    for (int ib = 0; ib < numrad; ib++) {

      rho[ib] = std::max(alb_pft.rhol[ib] * wl + alb_pft.rhos[ib] * ws, mpe);
      tau[ib] = std::max(alb_pft.taul[ib] * wl + alb_pft.taus[ib] * ws, mpe);

      // Calculate two-stream parameters omega, betad, and betai.
      // Omega, betad, betai are adjusted for snow. Values for omega*betad
      // and omega*betai are calculated and then divided by the new omega
      // because the product omega*betai, omega*betad is used in solution.
      // Also, the transmittances and reflectances (tau, rho) are linear
      // weights of leaf and stem values.
      const double omegal = rho[ib] + tau[ib];
      const double asu = 0.5 * omegal * gdir / temp0 * temp2;
      const double betadl = (1.0 + avmu * twostext) / (omegal * avmu * twostext) * asu;
      const double betail = 0.5 * ((rho[ib] + tau[ib]) + (rho[ib] - tau[ib]) * pow(((1.0 + chil) / 2.0), 2.0)) / omegal;

      // Adjust omega, betad, and betai for intercepted snow
      double tmp0, tmp1, tmp2;
      if (t_veg > ELM::constants::TFRZ) { // no snow
        tmp0 = omegal;
        tmp1 = betadl;
        tmp2 = betail;
      } else {
        tmp0 = (1.0 - fwet) * omegal + fwet * omegas[ib];
        tmp1 = ((1.0 - fwet) * omegal * betadl + fwet * omegas[ib] * betads) / tmp0;
        tmp2 = ((1.0 - fwet) * omegal * betail + fwet * omegas[ib] * betais) / tmp0;
      }
      omega[ib] = tmp0;
      const double betad = tmp1;
      const double betai = tmp2;

      // Common terms
      const double b = 1.0 - omega[ib] + omega[ib] * betai;
      const double c1 = omega[ib] * betai;
      tmp0 = avmu * twostext;
      const double d = tmp0 * omega[ib] * betad;
      const double f = tmp0 * omega[ib] * (1.0 - betad);
      tmp1 = b * b - c1 * c1;
      const double h = sqrt(tmp1) / avmu;
      const double sigma = tmp0 * tmp0 - tmp1;
      const double p1 = b + avmu * h;
      const double p2 = b - avmu * h;
      const double p3 = b + tmp0;
      const double p4 = b - tmp0;

      // Absorbed, reflected, transmitted fluxes per unit incoming radiation for full canopy
      double t1 = std::min(h * (elai + esai), 40.0);
      double s1 = exp(-t1);
      t1 = std::min(twostext * (elai + esai), 40.0);
      double s2 = exp(-t1);

      // Direct beam
      double u1 = b - c1 / albgrd(ib);
      double u2 = b - c1 * albgrd(ib);
      double u3 = f + c1 * albgrd(ib);
      tmp2 = u1 - avmu * h;
      double tmp3 = u1 + avmu * h;
      double d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1;
      double tmp4 = u2 + avmu * h;
      double tmp5 = u2 - avmu * h;
      double d2 = tmp4 / s1 - tmp5 * s1;
      double h1 = -d * p4 - c1 * f;
      double tmp6 = d - h1 * p3 / sigma;
      double tmp7 = (d - c1 - h1 / sigma * (u1 + tmp0)) * s2;
      double h2 = (tmp6 * tmp2 / s1 - p2 * tmp7) / d1;
      double h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1;
      double h4 = -f * p3 - c1 * d;
      double tmp8 = h4 / sigma;
      double tmp9 = (u3 - tmp8 * (u2 - tmp0)) * s2;
      double h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2;
      double h6 = (tmp8 * tmp5 * s1 + tmp9) / d2;

      albd(ib) = h1 / sigma + h2 + h3;
      ftid(ib) = h4 * s2 / sigma + h5 * s1 + h6 / s1;
      ftdd(ib) = s2;
      fabd(ib) = 1.0 - albd(ib) - (1.0 - albgrd(ib)) * ftdd(ib) - (1.0 - albgri(ib)) * ftid(ib);

      double a1 = h1 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) + h2 * (1.0 - s2 * s1) / (twostext + h) +
                  h3 * (1.0 - s2 / s1) / (twostext - h);
      double a2 = h4 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) + h5 * (1.0 - s2 * s1) / (twostext + h) +
                  h6 * (1.0 - s2 / s1) / (twostext - h);

      fabd_sun(ib) = (1.0 - omega[ib]) * (1.0 - s2 + 1.0 / avmu * (a1 + a2));
      fabd_sha(ib) = fabd(ib) - fabd_sun(ib);

      // Diffuse
      u1 = b - c1 / albgri(ib);
      u2 = b - c1 * albgri(ib);
      tmp2 = u1 - avmu * h;
      tmp3 = u1 + avmu * h;
      d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1;
      tmp4 = u2 + avmu * h;
      tmp5 = u2 - avmu * h;
      d2 = tmp4 / s1 - tmp5 * s1;
      double h7 = (c1 * tmp2) / (d1 * s1);
      double h8 = (-c1 * tmp3 * s1) / d1;
      double h9 = tmp4 / (d2 * s1);
      double h10 = (-tmp5 * s1) / d2;

      albi(ib) = h7 + h8;
      ftii(ib) = h9 * s1 + h10 / s1;
      fabi(ib) = 1.0 - albi(ib) - (1.0 - albgri(ib)) * ftii(ib);

      a1 = h7 * (1.0 - s2 * s1) / (twostext + h) + h8 * (1.0 - s2 / s1) / (twostext - h);
      a2 = h9 * (1.0 - s2 * s1) / (twostext + h) + h10 * (1.0 - s2 / s1) / (twostext - h);

      fabi_sun(ib) = (1.0 - omega[ib]) / avmu * (a1 + a2);
      fabi_sha(ib) = fabi(ib) - fabi_sun(ib);

      // Repeat two-stream calculations for each canopy layer to calculate derivatives.
      // tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives are
      // calculated at the center of the layer. Derivatives are needed only for the
      // visible waveband to calculate absorbed PAR (per unit lai+sai) for each canopy layer.
      // Derivatives are calculated first per unit lai+sai and then normalized for sunlit
      // or shaded fraction of canopy layer.
      // Sun/shade big leaf code uses only one layer, with canopy integrated values from above
      // and also canopy-integrated scaling coefficients

      if (ib == 0) {
        double laisum;
        if (nlevcan == 1) {
          // sunlit fraction of canopy
          fsun_z(0) = (1.0 - s2) / t1;

          // absorbed PAR (per unit sun/shade lai+sai)
          laisum = elai + esai;
          fabd_sun_z(0) = fabd_sun(ib) / (fsun_z(0) * laisum);
          fabi_sun_z(0) = fabi_sun(ib) / (fsun_z(0) * laisum);
          fabd_sha_z(0) = fabd_sha(ib) / ((1.0 - fsun_z(0)) * laisum);
          fabi_sha_z(0) = fabi_sha(ib) / ((1.0 - fsun_z(0)) * laisum);

          // leaf to canopy scaling coefficients
          double extkb = twostext;
          vcmaxcintsun = (1.0 - exp(-(extkn + extkb) * elai)) / (extkn + extkb);
          vcmaxcintsha = (1.0 - exp(-extkn * elai)) / extkn - vcmaxcintsun;
          if (elai > 0.0) {
            vcmaxcintsun = vcmaxcintsun / (fsun_z(0) * elai);
            vcmaxcintsha = vcmaxcintsha / ((1.0 - fsun_z(0)) * elai);
          } else {
            vcmaxcintsun = 0.0;
            vcmaxcintsha = 0.0;
          }

        } else if (nlevcan > 1) {
          for (int iv = 0; iv < nrad; iv++) {
            // Cumulative lai+sai at center of layer
            if (iv == 0) {
              laisum = 0.5 * (tlai_z(iv) + tsai_z(iv));
            } else {
              laisum += 0.5 * ((tlai_z(iv - 1) + tsai_z(iv - 1)) + (tlai_z(iv) + tsai_z(iv)));
            }

            // Coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction
            t1 = std::min(h * laisum, 40.0);
            s1 = exp(-t1);
            t1 = std::min(twostext * laisum, 40.0);
            s2 = exp(-t1);
            fsun_z(iv) = s2;

            // Direct beam
            // Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai
            u1 = b - c1 / albgrd(ib);
            u2 = b - c1 * albgrd(ib);
            u3 = f + c1 * albgrd(ib);
            tmp2 = u1 - avmu * h;
            tmp3 = u1 + avmu * h;
            d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1;
            tmp4 = u2 + avmu * h;
            tmp5 = u2 - avmu * h;
            d2 = tmp4 / s1 - tmp5 * s1;
            h1 = -d * p4 - c1 * f;
            tmp6 = d - h1 * p3 / sigma;
            tmp7 = (d - c1 - h1 / sigma * (u1 + tmp0)) * s2;
            h2 = (tmp6 * tmp2 / s1 - p2 * tmp7) / d1;
            h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1;
            h4 = -f * p3 - c1 * d;
            tmp8 = h4 / sigma;
            tmp9 = (u3 - tmp8 * (u2 - tmp0)) * s2;
            h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2;
            h6 = (tmp8 * tmp5 * s1 + tmp9) / d2;

            a1 = h1 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) + h2 * (1.0 - s2 * s1) / (twostext + h) +
                 h3 * (1.0 - s2 / s1) / (twostext - h);
            a2 = h4 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) + h5 * (1.0 - s2 * s1) / (twostext + h) +
                 h6 * (1.0 - s2 / s1) / (twostext - h);

            // Derivatives for h2, h3, h5, h6 and a1, a2
            double v = d1;
            double dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1;
            double u = tmp6 * tmp2 / s1 - p2 * tmp7;
            double du = h * tmp6 * tmp2 / s1 + twostext * p2 * tmp7;
            double dh2 = (v * du - u * dv) / (v * v);
            u = -tmp6 * tmp3 * s1 + p1 * tmp7;
            du = h * tmp6 * tmp3 * s1 - twostext * p1 * tmp7;
            double dh3 = (v * du - u * dv) / (v * v);
            v = d2;
            dv = h * tmp4 / s1 + h * tmp5 * s1;
            u = -h4 / sigma * tmp4 / s1 - tmp9;
            du = -h * h4 / sigma * tmp4 / s1 + twostext * tmp9;
            double dh5 = (v * du - u * dv) / (v * v);
            u = h4 / sigma * tmp5 * s1 + tmp9;
            du = -h * h4 / sigma * tmp5 * s1 - twostext * tmp9;
            double dh6 = (v * du - u * dv) / (v * v);

            double da1 = h1 / sigma * s2 * s2 + h2 * s2 * s1 + h3 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh2 +
                         (1.0 - s2 / s1) / (twostext - h) * dh3;
            double da2 = h4 / sigma * s2 * s2 + h5 * s2 * s1 + h6 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh5 +
                         (1.0 - s2 / s1) / (twostext - h) * dh6;

            // Flux derivatives
            double d_ftid = -twostext * h4 / sigma * s2 - h * h5 * s1 + h * h6 / s1 + dh5 * s1 + dh6 / s1;
            double d_fabd = -(dh2 + dh3) + (1.0 - albgrd(ib)) * twostext * s2 - (1.0 - albgri(ib)) * d_ftid;
            double d_fabd_sun = (1.0 - omega[ib]) * (twostext * s2 + 1.0 / avmu * (da1 + da2));
            double d_fabd_sha = d_fabd - d_fabd_sun;
            fabd_sun_z(iv) = std::max(d_fabd_sun, 0.0);
            fabd_sha_z(iv) = std::max(d_fabd_sha, 0.0);

            // Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
            // to normalize derivatives by sunlit or shaded fraction to get
            // APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
            fabd_sun_z(iv) = fabd_sun_z(iv) / fsun_z(iv);
            fabd_sha_z(iv) = fabd_sha_z(iv) / (1.0 - fsun_z(iv));

            // Diffuse
            // Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai
            u1 = b - c1 / albgri(ib);
            u2 = b - c1 * albgri(ib);
            tmp2 = u1 - avmu * h;
            tmp3 = u1 + avmu * h;
            d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1;
            tmp4 = u2 + avmu * h;
            tmp5 = u2 - avmu * h;
            d2 = tmp4 / s1 - tmp5 * s1;
            h7 = (c1 * tmp2) / (d1 * s1);
            h8 = (-c1 * tmp3 * s1) / d1;
            h9 = tmp4 / (d2 * s1);
            h10 = (-tmp5 * s1) / d2;

            a1 = h7 * (1.0 - s2 * s1) / (twostext + h) + h8 * (1.0 - s2 / s1) / (twostext - h);
            a2 = h9 * (1.0 - s2 * s1) / (twostext + h) + h10 * (1.0 - s2 / s1) / (twostext - h);

            // Derivatives for h7, h8, h9, h10 and a1, a2
            v = d1;
            dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1;
            u = c1 * tmp2 / s1;
            du = h * c1 * tmp2 / s1;
            double dh7 = (v * du - u * dv) / (v * v);
            u = -c1 * tmp3 * s1;
            du = h * c1 * tmp3 * s1;
            double dh8 = (v * du - u * dv) / (v * v);
            v = d2;
            dv = h * tmp4 / s1 + h * tmp5 * s1;
            u = tmp4 / s1;
            du = h * tmp4 / s1;
            double dh9 = (v * du - u * dv) / (v * v);
            u = -tmp5 * s1;
            du = h * tmp5 * s1;
            double dh10 = (v * du - u * dv) / (v * v);

            da1 = h7 * s2 * s1 + h8 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh7 +
                  (1.0 - s2 / s1) / (twostext - h) * dh8;
            da2 = h9 * s2 * s1 + h10 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh9 +
                  (1.0 - s2 / s1) / (twostext - h) * dh10;

            // Flux derivatives
            double d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1;
            double d_fabi = -(dh7 + dh8) - (1.0 - albgri(ib)) * d_ftii;
            double d_fabi_sun = (1.0 - omega[ib]) / avmu * (da1 + da2);
            double d_fabi_sha = d_fabi - d_fabi_sun;
            fabi_sun_z(iv) = std::max(d_fabi_sun, 0.0);
            fabi_sha_z(iv) = std::max(d_fabi_sha, 0.0);

            // Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
            // to normalize derivatives by sunlit or shaded fraction to get
            // APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
            fabi_sun_z(iv) = fabi_sun_z(iv) / fsun_z(iv);
            fabi_sha_z(iv) = fabi_sha_z(iv) / (1.0 - fsun_z(iv));
          } // for nrad
        }   // if nlevcan > 1
      }     // if ib == 0
    }       // for numrad
  } else if (novegsol(Land, coszen, elai, esai)) {
    for (int ib = 0; ib < numrad; ++ib) {
      fabd(ib) = 0.0;
      fabd_sun(ib) = 0.0;
      fabd_sha(ib) = 0.0;
      fabi(ib) = 0.0;
      fabi_sun(ib) = 0.0;
      fabi_sha(ib) = 0.0;
      ftdd(ib) = 1.0;
      ftid(ib) = 0.0;
      ftii(ib) = 1.0;
      albd(ib) = albgrd(ib);
      albi(ib) = albgri(ib);
    }
  }
} // two_stream_solver

template <class ArrayD1>
void soil_albedo(const LandType& Land, const int& snl, const double& t_grnd, const double& coszen,
                 const ArrayD1 h2osoi_vol, const ArrayD1 albsat, const ArrayD1 albdry, ArrayD1 albsod, ArrayD1 albsoi)
{

  if (!Land.urbpoi) {
    if (coszen > 0.0) {
      for (int ib = 0; ib < numrad; ib++) {
        if (Land.ltype == istsoil || Land.ltype == istcrop) {
          double inc = std::max(0.11 - 0.40 * h2osoi_vol(0), 0.0); // soil water correction factor for soil albedo
          // double soilcol = isoicol;
          albsod(ib) = std::min(albsat(ib) + inc, albdry(ib));
          albsoi(ib) = albsod(ib);
        } else if (Land.ltype == istice || Land.ltype == istice_mec) {
          albsod(ib) = albice[ib];
          albsoi(ib) = albsod[ib];
          // comment out lake logic for now
          // } else if (t_grnd > ELM::constants::TFRZ || (lakepuddling && Land.ltype == istdlak && t_grnd == ELM::constants::TFRZ && lake_icefrac(0)
          // < 1.0 &&
          //                              lake_icefrac(1) > 0.0)) { // maybe get rid of lake logic?
          //   albsod(ib) = 0.05 / (std::max(0.001, coszen) + 0.15);
          //   // This expression is apparently from BATS according to Yongjiu Dai.
          //   // The diffuse albedo should be an average over the whole sky of an angular-dependent direct expression.
          //   // The expression above may have been derived to encompass both (e.g. Henderson-Sellers 1986),
          //   // but I'll assume it applies more appropriately to the direct form for now.
          //   // ZMS: Attn EK, currently restoring this for wetlands even though it is wrong in order to try to get
          //   // bfb baseline comparison when no lakes are present. I'm assuming wetlands will be phased out anyway.
          //   if (Land.ltype == istdlak) {
          //     albsoi(ib) = 0.10;
          //   } else {
          //     albsoi(ib) = albsod(ib);
          //   }

        } else {
          // frozen lake, wetland
          // Introduce crude surface frozen fraction according to D. Mironov (2010)
          // Attn EK: This formulation is probably just as good for "wetlands" if they are not phased out.
          // Tenatively I'm restricting this to lakes because I haven't tested it for wetlands. But if anything
          // the albedo should be lower when melting over frozen ground than a solid frozen lake.
          if (Land.ltype == istdlak && !lakepuddling && snl == 0) {
            // Need to reference snow layers here because t_grnd could be over snow or ice
            // but we really want the ice surface temperature with no snow
            double sicefr = 1.0 - exp(-calb * (ELM::constants::TFRZ - t_grnd) / ELM::constants::TFRZ);
            albsod(ib) =
                sicefr * alblak[ib] + (1.0 - sicefr) * std::max(alblakwi[ib], 0.05 / (std::max(0.001, coszen) + 0.15));
            albsoi(ib) = sicefr * alblak[ib] + (1.0 - sicefr) * std::max(alblakwi[ib], 0.10);
            // Make sure this is no less than the open water albedo above.
            // Setting alblakwi(:) = alblak(:) reverts the melting albedo to the cold
            // snow-free value.
          } else {
            albsod(ib) = alblak[ib];
            albsoi(ib) = albsod(ib);
          }
        }
      }
    }
  }
} // soil_albedo

} // namespace ELM::surface_albedo
