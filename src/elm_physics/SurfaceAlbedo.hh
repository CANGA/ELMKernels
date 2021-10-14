#include <cmath>
#include <stdexcept>

namespace ELM {

// need to make metric to mimic num_vegsol num_novegsol
// need to write coszen function
// need to initialize h2osoi_vol[nlevgrnd]
// need to read in soil color NetCDF, initialize isoicol


// initial values from SurfaceAlbedoType::InitCold
this%albgrd_col     (begc:endc, :) = 0.2_r8
this%albgri_col     (begc:endc, :) = 0.2_r8
this%albsod_col     (begc:endc, :) = 0.2_r8
this%albsoi_col     (begc:endc, :) = 0.2_r8
this%albsnd_hst_col (begc:endc, :) = 0.6_r8
this%albsni_hst_col (begc:endc, :) = 0.6_r8
this%albd_patch     (begp:endp, :) = 0.2_r8
this%albi_patch     (begp:endp, :) = 0.2_r8
this%albgrd_pur_col (begc:endc, :) = 0.2_r8
this%albgri_pur_col (begc:endc, :) = 0.2_r8
this%albgrd_bc_col  (begc:endc, :) = 0.2_r8
this%albgri_bc_col  (begc:endc, :) = 0.2_r8
this%albgrd_oc_col  (begc:endc, :) = 0.2_r8
this%albgri_oc_col  (begc:endc, :) = 0.2_r8
this%albgrd_dst_col (begc:endc, :) = 0.2_r8
this%albgri_dst_col (begc:endc, :) = 0.2_r8
this%fabi_patch     (begp:endp, :) = 0.0_r8
this%fabd_patch     (begp:endp, :) = 0.0_r8
this%fabi_sun_patch (begp:endp, :) = 0.0_r8
this%fabd_sun_patch (begp:endp, :) = 0.0_r8
this%fabd_sha_patch (begp:endp, :) = 0.0_r8
this%fabi_sha_patch (begp:endp, :) = 0.0_r8
this%ftdd_patch     (begp:endp, :) = 1.0_r8
this%ftid_patch     (begp:endp, :) = 0.0_r8
this%ftii_patch     (begp:endp, :) = 1.0_r8


// FUNCTION to return the cosine of the solar zenith angle.
// Assumes 365.0 days/year
//jday    Julian cal day (1.xx to 365.xx)
//lat     Centered latitude (radians)
//lon     Centered longitude (radians)
//declin  Solar declination (radians)
double cosz(const double jday, const double lat, const double lon, const double declin) {

  return orb_cosz = sin(lat) * sin(declin) - 
  cos(lat) * cos(declin) * cos((jday - floor(jday)) * 2.0 * pi + lon);
}

void SurfaceAlbedo() {

  //CosZen() - need to write - get coszen - elm gets for next time step 
  // figure out if there's a reason to do this in kernels - calling SurfaceAlbedo() at end of step and calculating for t+1
  // or if we can call SurfaceAlbedo() first and calc for t

  const double mpe = 1.e-06; // prevents overflow for division by zero

  // Cosine solar zenith angle for next time step
  double coszen = cosz(nextsw_cday, lat, lon, declinp1);

  // Initialize output because solar radiation only done if coszen > 0

  if (!Land.urbpoi) {
    for (int ib=0; ib < numrad; ++ib) {
      // column vars
      albsod[ib]     = 0.0;
      albsoi[ib]     = 0.0;
      albgrd[ib]     = 0.0;
      albgri[ib]     = 0.0;
      albgrd_pur[ib] = 0.0;
      albgri_pur[ib] = 0.0;
      albgrd_bc[ib]  = 0.0;
      albgri_bc[ib]  = 0.0;
      albgrd_oc[ib]  = 0.0;
      albgri_oc[ib]  = 0.0;
      albgrd_dst[ib] = 0.0;
      albgri_dst[ib] = 0.0;
      for (int i = 0; i <= nlevsno; ++i) {
        flx_absdv[i] = 0.0;
        flx_absdn[i] = 0.0;
        flx_absiv[i] = 0.0;
        flx_absin[i] = 0.0;
      }

      // patch vars
      albd[ib] = 1.0;
      albi[ib] = 1.0;
      fabd[ib] = 0.0;
      fabd_sun[ib] = 0.0;
      fabd_sha[ib] = 0.0;
      fabi[ib] = 0.0;
      fabi_sun[ib] = 0.0;
      fabi_sha[ib] = 0.0;
      ftdd[ib] = 0.0;
      ftid[ib] = 0.0;
      ftii[ib] = 0.0;
    }
  }


  // SoilAlbedo called before SNICAR_RT/SNICAR_AD_RT
  // so that reflectance of soil beneath snow column is known
  // ahead of time for snow RT calculation.

  // Snow albedos
  // Note that snow albedo routine will only compute nonzero snow albedos
  // where h2osno> 0 and coszen > 0

  // Ground surface albedos
  // Note that ground albedo routine will only compute nonzero snow albedos
  // where coszen > 0

  // Place SoilAlbedo call here. - maybe bring out to driver level so that the programming model can schedule it

  //call SoilAlbedo(bounds, &
  //     num_nourbanc, filter_nourbanc, &
  //  coszen_col(bounds%begc:bounds%endc), &
  //  albsnd(bounds%begc:bounds%endc, :), &
  //     albsni(bounds%begc:bounds%endc, :), &  
  //     lakestate_vars, temperature_vars, waterstate_vars, surfalb_vars)

  // set variables to pass to SNICAR.
  for (int i = 0; i < numrad; ++i) { albsfc[i]= albsoi[i]; }
  for (int i = 0; i < nlevsno; ++i) {
    h2osno_liq[i] = h2osoi_liq[i];
    h2osno_ice[i] = h2osoi_ice[i];
    snw_rds_in[i] = round(snw_rds[i]);
  }


  //// zero aerosol input arrays
  //do aer = 1, sno_nbr_aer
  //   do i = -nlevsno+1, 0
  //      do c = bounds%begc, bounds%endc
  //         mss_cnc_aer_in_frc_pur(c,i,aer) = 0._r8
  //         mss_cnc_aer_in_frc_bc(c,i,aer)  = 0._r8
  //         mss_cnc_aer_in_frc_oc(c,i,aer)  = 0._r8
  //         mss_cnc_aer_in_frc_dst(c,i,aer) = 0._r8
  //         mss_cnc_aer_in_fdb(c,i,aer)     = 0._r8
  //      end do
  //   end do
  //end do

  //! Set aerosol input arrays
  //! feedback input arrays have been zeroed
  //! set soot and dust aerosol concentrations:
  //if (DO_SNO_AER) then
  //   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)

  //   ! DO_SNO_OC is set in SNICAR_varpar. Default case is to ignore OC concentrations because:
  //   !  1) Knowledge of their optical properties is primitive
  //   !  2) When 'water-soluble' OPAC optical properties are applied to OC in snow,
  //   !     it has a negligible darkening effect.
  //   if (DO_SNO_OC) then
  //      mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
  //   endif

  //   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)
  //endif

  //! If radiative forcing is being calculated, first estimate clean-snow albedo

  //if (use_snicar_frc) then
  //   ! 1. BC input array:
  //   !  set dust and (optionally) OC concentrations, so BC_FRC=[(BC+OC+dust)-(OC+dust)]
  //   mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)
  //   if (DO_SNO_OC) then
  //      mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
  //   endif

  //   ! BC FORCING CALCULATIONS
  //      flg_slr = 1; ! direct-beam
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsnd_bc(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsnd_bc(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad

  //      flg_slr = 2; ! diffuse
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsni_bc(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsni_bc(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad

  //   ! 2. OC input array:
  //   !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
  //   if (DO_SNO_OC) then
  //      mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)

  //   ! OC FORCING CALCULATIONS
  //      flg_slr = 1; ! direct-beam
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsnd_oc(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsnd_oc(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad

  //      flg_slr = 2; ! diffuse
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsni_oc(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsni_oc(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad
  //   endif

  //   ! 3. DUST input array:
  //   ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
  //   mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
  //   mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)
  //   if (DO_SNO_OC) then
  //      mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
  //      mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
  //   endif

  //   ! DUST FORCING CALCULATIONS
  //      flg_slr = 1; ! direct-beam
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsnd_dst(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsnd_dst(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad

  //      flg_slr = 2; ! diffuse
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsni_dst(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsni_dst(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad


  //   ! 4. ALL AEROSOL FORCING CALCULATION
  //   ! (pure snow albedo)
  //      flg_slr = 1; ! direct-beam
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsnd_pur(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsnd_pur(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad

  //      flg_slr = 2; ! diffuse
  //      if (use_snicar_ad) then
  //          call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                            coszen_col(bounds%begc:bounds%endc), &
  //                            flg_slr, &
  //                            h2osno_liq(bounds%begc:bounds%endc, :), &
  //                            h2osno_ice(bounds%begc:bounds%endc, :), &
  //                            snw_rds_in(bounds%begc:bounds%endc, :), &
  //                            mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
  //                            albsfc(bounds%begc:bounds%endc, :), &
  //                            albsni_pur(bounds%begc:bounds%endc, :), &
  //                            foo_snw(bounds%begc:bounds%endc, :, :), &
  //                            waterstate_vars)
  //      else
  //          call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                         coszen_col(bounds%begc:bounds%endc), &
  //                         flg_slr, &
  //                         h2osno_liq(bounds%begc:bounds%endc, :), &
  //                         h2osno_ice(bounds%begc:bounds%endc, :), &
  //                         snw_rds_in(bounds%begc:bounds%endc, :), &
  //                         mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
  //                         albsfc(bounds%begc:bounds%endc, :), &
  //                         albsni_pur(bounds%begc:bounds%endc, :), &
  //                         foo_snw(bounds%begc:bounds%endc, :, :), &
  //                         waterstate_vars)
  //      endif ! end if use_snicar_ad
  //end if !end if use_snicar_frc

  // CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
  flg_slr = 0; // direct-beam
  //if (use_snicar_ad) then
  //    call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                      coszen_col(bounds%begc:bounds%endc), &
  //                      flg_slr, &
  //                      h2osno_liq(bounds%begc:bounds%endc, :), &
  //                      h2osno_ice(bounds%begc:bounds%endc, :), &
  //                      snw_rds_in(bounds%begc:bounds%endc, :), &
  //                      mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
  //                      albsfc(bounds%begc:bounds%endc, :), &
  //                      albsnd(bounds%begc:bounds%endc, :), &
  //                      flx_absd_snw(bounds%begc:bounds%endc, :, :), &
  //                      waterstate_vars)
  //else
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsnd(bounds%begc:bounds%endc, :), &
                      flx_absd_snw(bounds%begc:bounds%endc, :, :), &
                      waterstate_vars)
  //endif ! end if use_snicar_ad

  flg_slr = 1; // diffuse
  //if (use_snicar_ad) then
  //    call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
  //                      coszen_col(bounds%begc:bounds%endc), &
  //                      flg_slr, &
  //                      h2osno_liq(bounds%begc:bounds%endc, :), &
  //                      h2osno_ice(bounds%begc:bounds%endc, :), &
  //                      snw_rds_in(bounds%begc:bounds%endc, :), &
  //                      mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
  //                      albsfc(bounds%begc:bounds%endc, :), &
  //                      albsni(bounds%begc:bounds%endc, :), &
  //                      flx_absi_snw(bounds%begc:bounds%endc, :, :), &
  //                      waterstate_vars)
  //else
       call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                      coszen_col(bounds%begc:bounds%endc), &
                      flg_slr, &
                      h2osno_liq(bounds%begc:bounds%endc, :), &
                      h2osno_ice(bounds%begc:bounds%endc, :), &
                      snw_rds_in(bounds%begc:bounds%endc, :), &
                      mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                      albsfc(bounds%begc:bounds%endc, :), &
                      albsni(bounds%begc:bounds%endc, :), &
                      flx_absi_snw(bounds%begc:bounds%endc, :, :), &
                      waterstate_vars)
  //endif ! end if use_snicar_ad

  // ground albedos and snow-fraction weighting of snow absorption factors
  if (!Land.urbpoi) {
    for (int ib = 0; ib < nband; ++ib) {
      if (coszen > 0.0) {
        // ground albedo was originally computed in SoilAlbedo, but is now computed here
        // because the order of SoilAlbedo and SNICAR_RT/SNICAR_AD_RT was switched for SNICAR/SNICAR_AD_RT.
        albgrd[ib] = albsod[ib] * (1.0 - frac_sno) + albsnd[ib] * frac_sno;
        albgri[ib] = albsoi[ib] * (1.0 - frac_sno) + albsni[ib] * frac_sno;

        // also in this loop (but optionally in a different loop for vectorized code)
        // weight snow layer radiative absorption factors based on snow fraction and soil albedo
        // (NEEDED FOR ENERGY CONSERVATION)
        for (int i = 0; i <= nlevsno; ++i) {
          if (subgridflag == 0 || Land.ltype == istdlak) {
            if (ib == 0) {
              flx_absdv[i] = flx_absd_snw[i][ib] * frac_sno +
                    ((1.0 - frac_sno) * (1.0 - albsod[ib]) * (flx_absd_snw[i][ib] / (1.0 - albsnd[ib])));
              flx_absiv[i] = flx_absi_snw[i][ib] * frac_sno +
                    ((1.0 - frac_sno) * (1.0 - albsoi[ib]) * (flx_absi_snw[i][ib] / (1.0 - albsni[ib])));
            } else if (ib == 1) {
              flx_absdn[i] = flx_absd_snw[i][ib] * frac_sno +
                    ((1.0 - frac_sno) * (1.0 - albsod[ib]) * (flx_absd_snw[i][ib] / (1.0 - albsnd[ib])));
              flx_absin[i] = flx_absi_snw[i][ib] * frac_sno +
                    ((1.0 - frac_sno) * (1.0 - albsoi[ib]) * (flx_absi_snw[i][ib] / (1.0 - albsni[ib])));
                }
          } else {
            if (ib == 0) {
              flx_absdv[i] = flx_absd_snw[i][ib] * (1.0 - albsnd[ib]);
              flx_absiv[i] = flx_absi_snw[i][ib] * (1.0 - albsni[ib]);
            } else if (ib == 1) {
              flx_absdn[i] = flx_absd_snw[i][ib] * (1.0 - albsnd[ib]);
              flx_absin[i] = flx_absi_snw[i][ib] * (1.0 - albsni[ib]);
            }
          } // if subgridflag
        } // for nlensno +1 layers
      } // if coszen > 0
    } // for nband
  } // if !Land.urbpoi

  // Create solar-vegetated filter for the following calculations
  num_vegsol = 0;
  num_novegsol = 0;
  if (!Land.urbpoi) {
    if (coszen > 0.0) {
      if ((Land.ltype == istsoil || Land.ltype == istcrop) && (elai + esai) > 0.0) {
        num_vegsol = 1;
      } else {
        num_novegsol = 1;
      }
    }
  }


  // Weight reflectance/transmittance by lai and sai
  // Only perform on vegetated patches where coszen > 0
  if (num_vegsol) {
    wl = elai / max(elai + esai, mpe);
    ws = esai / max(elai + esai, mpe);

    for (int ib = 0; ib < numrad; ++ib) {
      ho[ib] = max(rhol[Land.vtype][ib] * wl + rhos[Land.vtype][ib] * ws, mpe);
      tau[ib] = max(taul[Land.vtype][ib] * wl + taus[Land.vtype][ib] * ws, mpe);
    }
  }


  // Diagnose number of canopy layers for radiative transfer, in increments of dincmax.
  // Add to number of layers so long as cumulative leaf+stem area does not exceed total
  // leaf+stem area. Then add any remaining leaf+stem area to next layer and exit the loop.
  // Do this first for elai and esai (not buried by snow) and then for the part of the
  // canopy that is buried by snow.
  // ------------------
  // tlai_z = leaf area increment for a layer
  // tsai_z = stem area increment for a layer
  // nrad   = number of canopy layers above snow
  // ncan   = total number of canopy layers
  //
  // tlai_z summed from 1 to nrad = elai
  // tlai_z summed from 1 to ncan = tlai

  // tsai_z summed from 1 to nrad = esai
  // tsai_z summed from 1 to ncan = tsai
  // ------------------
  //
  // Canopy layering needs to be done for all "num_nourbanp" not "num_vegsol"
  // because layering is needed for all time steps regardless of radiation
  //
  // Sun/shade big leaf code uses only one layer (nrad = ncan = 1), triggered by
  // nlevcan = 1

  const double dincmax = 0.25;
  if (!Land.urbpoi) {

    if (nlevcan == 1) {
      nrad = 1;
      ncan = 1;
      tlai_z[0] = elai;
      tsai_z[0] = esai;
    } else if (nlevcan > 1) {
      if (elai + esai == 0.0) {
        nrad = 0;
      } else {
        dincmax_sum = 0.0;
        for (int iv = 0; iv < nlevcan; ++iv) {
          dincmax_sum += dincmax;
          if (((elai + esai) - dincmax_sum) > 1.e-06) {
            nrad = iv + 1;
            dinc = dincmax;
            tlai_z[iv] = dinc * elai / max(elai + esai, mpe);
            tsai_z[iv] = dinc * esai / max(elai + esai, mpe);
          } else {
            nrad = iv + 1;
            dinc = dincmax - (dincmax_sum - (elai + esai));
            tlai_z[iv] = dinc * elai / max(elai + esai, mpe);
            tsai_z[iv] = dinc * esai / max(elai + esai, mpe);
            break;
          }
        }
        // Mimumum of 4 canopy layers
        if (nrad < 4) {
          nrad = 4;
          for (int iv = 0; iv < nrad; ++iv) {
            tlai_z[iv] = elai / nrad;
            tsai_z[iv] = esai / nrad;
          }
        }
      }
    }

    // Error check: make sure cumulative of increments does not exceed total
    laisum = 0.0;
    saisum = 0.0;
    for (int iv = 0; iv < nrad; ++iv) {
      laisum += tlai_z[iv];
      saisum += tsai_z[iv];
    }
    
    if (abs(laisum - elai) > 1.e-06 || abs(saisum - esai) > 1.e-06) {
      throw std::runtime_error("ELM ERROR: multi-layer canopy error 1 in SurfaceAlbedo");
    }

    // Repeat to find canopy layers buried by snow
    if (nlevcan > 1) {
      blai = tlai - elai;
      bsai = tsai - esai;
      if (blai + bsai == 0.0) {
         ncan = nrad;
      } else {
        dincmax_sum = 0.0;
        for (int iv = nrad; iv < nlevcan; ++iv) {
          dincmax_sum += dincmax;
          if (((blai + bsai) - dincmax_sum) > 1.e-06) {
            ncan = iv + 1;
            dinc = dincmax;
            tlai_z[iv] = dinc * blai / max(blai + bsai, mpe);
            tsai_z[iv] = dinc * bsai / max(blai + bsai, mpe);
          } else {
            ncan = iv + 1;
            dinc = dincmax - (dincmax_sum - (blai + bsai));
            tlai_z[iv] = dinc * blai / max(blai + bsai, mpe);
            tsai_z[iv] = dinc * bsai / max(blai + bsai, mpe);
            break;
          }
        }
      }

      // Error check: make sure cumulative of increments does not exceed total
      laisum = 0.0;
      saisum = 0.0;
      for (int iv = 0; iv < ncan; ++iv) {
        laisum += tlai_z[iv];
        saisum += tsai_z[iv];
      }
      do iv = 1, ncan(p)
         laisum = laisum + tlai_z(p,iv)
         saisum = saisum + tsai_z(p,iv)
      end do
      if (abs(laisum - tlai) > 1.e-06 || abs(saisum - tsai) > 1.e-06) {
        throw std::runtime_error("ELM ERROR: multi-layer canopy error 1 in SurfaceAlbedo");
      }
    }


    for (int iv = 0; iv < nrad; ++iv) {
      fabd_sun_z[iv] = 0.0;
      fabd_sha_z[iv] = 0.0;
      fabi_sun_z[iv] = 0.0;
      fabi_sha_z[iv] = 0.0;
      fsun_z[iv] = 0.0;
    }

    // Default leaf to canopy scaling coefficients, used when coszen <= 0.
    // This is the leaf nitrogen profile integrated over the full canopy.
    // Integrate exp(-kn*x) over x=0 to x=elai and assign to shaded canopy,
    // because sunlit fraction is 0. Canopy scaling coefficients are set in
    // TwoStream for coszen > 0. So kn must be set here and in TwoStream.

    const double extkn = 0.30;
    if (nlevcan == 1) {
      vcmaxcintsun = 0.0;
      vcmaxcintsha = (1.0 - exp(-extkn * elai)) / extkn;
      if (elai > 0.0) { 
        vcmaxcintsha = vcmaxcintsha / elai;
      } else {
        vcmaxcintsha = 0.0;
      }
    } else if (nlevcan > 1) {
       vcmaxcintsun = 0.0;
       vcmaxcintsha = 0.0;
    }


    // Calculate surface albedos and fluxes
    // Only perform on vegetated patches where coszen > 0
    //call TwoStream (bounds, filter_vegsol, num_vegsol, &
    //     coszen_patch(bounds%begp:bounds%endp), &
    //     rho(bounds%begp:bounds%endp, :), &
    //     tau(bounds%begp:bounds%endp, :), &
    //     canopystate_vars, temperature_vars, waterstate_vars, surfalb_vars)


    // Determine values for non-vegetated patches where coszen > 0

    do ib = 1,numrad
    if (num_novegsol) {
      for (int ib = 0; ib < numrad; ++ib) {
        fabd[ib] = 0.0;
        fabd_sun[ib] = 0.0;
        fabd_sha[ib] = 0.0;
        fabi[ib] = 0.0;
        fabi_sun[ib] = 0.0;
        fabi_sha[ib] = 0.0;
        ftdd[ib] = 1.0;
        ftid[ib] = 0.0;
        ftii[ib] = 1.0;
        albd[ib] = albgrd[ib];
        albi[ib] = albgri[ib];
      }
    }
  }
}























/*
DESCRIPTION:
Two-stream fluxes for canopy radiative transfer
Use two-stream approximation of Dickinson (1983) Adv Geophysics
25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
to calculate fluxes absorbed by vegetation, reflected by vegetation,
and transmitted through vegetation for unit incoming direct or diffuse
flux given an underlying surface with known albedo.
Calculate sunlit and shaded fluxes as described by
Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
a multi-layer canopy to calculate APAR profile

*/
void TwoStream() {
  const int omegas[numrad] = {0.8, 0.4};

  if (num_vegsol) {
    // note that the following limit only acts on cosz values > 0 and less than 0.001, not on values cosz = 0, since
    // these zero have already been filtered out in filter_vegsol
    cosz = std::max(0.001, coszen);
    chil = std::min(std::max(xl[vtype], -0.4), 0.6);
    if (std::abs(chil) <= 0.01) {
      chil = 0.01;
    }
    phi1 = 0.5 - 0.633 * chil - 0.330 * chil * chil;
    phi2 = 0.877 * (1.0 - 2.0 * phi1);
    gdir = phi1 + phi2 * cosz;
    twostext = gdir / cosz;
    avmu = (1.0 - phi1 / phi2 * std::log((phi1 + phi2) / phi1)) / phi2;
    temp0 = gdir + phi2 * cosz;
    temp1 = phi1 * cosz;
    temp2 = (1.0 - temp1 / temp0 * std::log((temp1 + temp0) / temp1));
  }

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

  if (num_vegsol) {
    for (int ib = 0; ib < numrad; ib++) {
      // Calculate two-stream parameters omega, betad, and betai.
      // Omega, betad, betai are adjusted for snow. Values for omega*betad
      // and omega*betai are calculated and then divided by the new omega
      // because the product omega*betai, omega*betad is used in solution.
      // Also, the transmittances and reflectances (tau, rho) are linear
      // weights of leaf and stem values.
      omegal = rho[ib] + tau[ib];
      asu = 0.5 * omegal * gdir / temp0 * temp2;
      betadl = (1.0 + avmu * twostext) / (omegal * avmu * twostext) * asu;
      betail = 0.5 * ((rho[ib] + tau[ib]) + (rho[ib] - tau[ib]) * pow(((1.0 + chil) / 2.0), 2.0)) / omegal;

      // Adjust omega, betad, and betai for intercepted snow
      if (t_veg > tfrz) { // no snow
        tmp0 = omegal;
        tmp1 = betadl;
        tmp2 = betail;
      } else {
        tmp0 = (1.0 - fwet) * omegal + fwet * omegas[ib];
        tmp1 = ((1.0 - fwet) * omegal * betadl + fwet * omegas[ib] * betads) / tmp0;
        tmp2 = ((1.0 - fwet) * omegal * betail + fwet * omegas[ib] * betais) / tmp0;
      }
      omega[ib] = tmp0;
      betad = tmp1;
      betai = tmp2;

      // Common terms
      b = 1.0 - omega[ib] + omega[ib] * betai;
      c1 = omega[ib] * betai;
      tmp0 = avmu * twostext;
      d = tmp0 * omega[ib] * betad;
      f = tmp0 * omega[ib] * (1.0 - betad);
      tmp1 = b * b - c1 * c1;
      h = sqrt(tmp1) / avmu;
      sigma = tmp0 * tmp0 - tmp1;
      p1 = b + avmu * h;
      p2 = b - avmu * h;
      p3 = b + tmp0;
      p4 = b - tmp0;

      // Absorbed, reflected, transmitted fluxes per unit incoming radiation for full canopy
      t1 = std::min(h * (elai + esai), 40.0);
      s1 = exp(-t1);
      t1 = std::min(twostext * (elai + esai), 40.0);
      s2 = exp(-t1);

      // Direct beam
      u1 = b - c1 / albgrd[ib] u2 = b - c1 * albgrd[ib] u3 = f + c1 * albgrd[ib] tmp2 = u1 - avmu * h;
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

      albd[ib] = h1 / sigma + h2 + h3;
      ftid[ib] = h4 * s2 / sigma + h5 * s1 + h6 / s1;
      ftdd[ib] = s2;
      fabd[ib] = 1.0 - albd[ib] - (1.0 - albgrd[ib]) * ftdd[ib] - (1.0 - albgri[ib]) * ftid[ib];

      a1 = h1 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) + h2 * (1.0 - s2 * s1) / (twostext + h) +
           h3 * (1.0 - s2 / s1) / (twostext - h);
      a2 = h4 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) + h5 * (1.0 - s2 * s1) / (twostext + h) +
           h6 * (1.0 - s2 / s1) / (twostext - h);

      fabd_sun[ib] = (1.0 - omega[ib]) * (1.0 - s2 + 1.0 / avmu * (a1 + a2));
      fabd_sha[ib] = fabd[ib] - fabd_sun[ib];

      // Diffuse
      u1 = b - c1 / albgri[ib];
      u2 = b - c1 * albgri[ib];
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

      albi[ib] = h7 + h8;
      ftii[ib] = h9 * s1 + h10 / s1;
      fabi[ib] = 1.0 - albi[ib] - (1.0 - albgri[ib]) * ftii[ib];

      a1 = h7 * (1.0 - s2 * s1) / (twostext + h) + h8 * (1.0 - s2 / s1) / (twostext - h);
      a2 = h9 * (1.0 - s2 * s1) / (twostext + h) + h10 * (1.0 - s2 / s1) / (twostext - h);

      fabi_sun[ib] = (1.0 - omega[ib]) / avmu * (a1 + a2);
      fabi_sha[ib] = fabi[ib] - fabi_sun[ib];

      // Repeat two-stream calculations for each canopy layer to calculate derivatives.
      // tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives are
      // calculated at the center of the layer. Derivatives are needed only for the
      // visible waveband to calculate absorbed PAR (per unit lai+sai) for each canopy layer.
      // Derivatives are calculated first per unit lai+sai and then normalized for sunlit
      // or shaded fraction of canopy layer.
      // Sun/shade big leaf code uses only one layer, with canopy integrated values from above
      // and also canopy-integrated scaling coefficients

      if (ib == 0) {
        if (nlevcan == 1) {
          // sunlit fraction of canopy
          fsun_z[0] = (1.0 - s2) / t1;

          // absorbed PAR (per unit sun/shade lai+sai)
          laisum = elai + esai;
          fabd_sun_z[0] = fabd_sun[ib] / (fsun_z[0] * laisum);
          fabi_sun_z[0] = fabi_sun[ib] / (fsun_z[0] * laisum);
          fabd_sha_z[0] = fabd_sha[ib] / ((1.0 - fsun_z[0]) * laisum);
          fabi_sha_z[0] = fabi_sha[ib] / ((1.0 - fsun_z[0]) * laisum);

          // leaf to canopy scaling coefficients
          extkn = 0.30;
          extkb = twostext;
          vcmaxcintsun = (1.0 - exp(-(extkn + extkb) * elai)) / (extkn + extkb);
          vcmaxcintsha = (1.0 - exp(-extkn * elai)) / extkn - vcmaxcintsun;
          if (elai > 0.0) {
            vcmaxcintsun = vcmaxcintsun / (fsun_z[0] * elai);
            vcmaxcintsha = vcmaxcintsha / ((1.0 - fsun_z[0]) * elai);
          } else {
            vcmaxcintsun = 0.0;
            vcmaxcintsha = 0.0;
          }

        } else if (nlevcan > 1) {
          for (int iv = 0; iv < nrad; iv++) {
            // Cumulative lai+sai at center of layer
            if (iv == 0) {
              laisum = 0.5 * (tlai_z[iv] + tsai_z[iv]);
            } else {
              laisum = laisum + 0.5 * ((tlai_z[iv - 1] + tsai_z[iv - 1]) + (tlai_z[iv] + tsai_z[iv]));
            }

            // Coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction
            t1 = std::min(h * laisum, 40.0);
            s1 = exp(-t1);
            t1 = std::min(twostext * laisum, 40.0);
            s2 = exp(-t1);
            fsun_z[iv] = s2;

            // Direct beam
            // Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai
            u1 = b - c1 / albgrd[ib];
            u2 = b - c1 * albgrd[ib];
            u3 = f + c1 * albgrd[ib];
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
            v = d1;
            dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1;
            u = tmp6 * tmp2 / s1 - p2 * tmp7;
            du = h * tmp6 * tmp2 / s1 + twostext * p2 * tmp7;
            dh2 = (v * du - u * dv) / (v * v);
            u = -tmp6 * tmp3 * s1 + p1 * tmp7;
            du = h * tmp6 * tmp3 * s1 - twostext * p1 * tmp7;
            dh3 = (v * du - u * dv) / (v * v);
            v = d2;
            dv = h * tmp4 / s1 + h * tmp5 * s1;
            u = -h4 / sigma * tmp4 / s1 - tmp9;
            du = -h * h4 / sigma * tmp4 / s1 + twostext * tmp9;
            dh5 = (v * du - u * dv) / (v * v);
            u = h4 / sigma * tmp5 * s1 + tmp9;
            du = -h * h4 / sigma * tmp5 * s1 - twostext * tmp9;
            dh6 = (v * du - u * dv) / (v * v);

            da1 = h1 / sigma * s2 * s2 + h2 * s2 * s1 + h3 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh2 +
                  (1.0 - s2 / s1) / (twostext - h) * dh3;
            da2 = h4 / sigma * s2 * s2 + h5 * s2 * s1 + h6 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh5 +
                  (1.0 - s2 / s1) / (twostext - h) * dh6;

            // Flux derivatives
            d_ftid = -twostext * h4 / sigma * s2 - h * h5 * s1 + h * h6 / s1 + dh5 * s1 + dh6 / s1;
            d_fabd = -(dh2 + dh3) + (1.0 - albgrd[ib]) * twostext * s2 - (1.0 - albgri[ib]) * d_ftid;
            d_fabd_sun = (1.0 - omega[ib]) * (twostext * s2 + 1.0 / avmu * (da1 + da2));
            d_fabd_sha = d_fabd - d_fabd_sun;
            fabd_sun_z[iv] = std::max(d_fabd_sun, 0.0);
            fabd_sha_z[iv] = std::max(d_fabd_sha, 0.0);

            // Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
            // to normalize derivatives by sunlit or shaded fraction to get
            // APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
            fabd_sun_z[iv] = fabd_sun_z[iv] / fsun_z[iv];
            fabd_sha_z[iv] = fabd_sha_z[iv] / (1.0 - fsun_z[iv]);

            // Diffuse
            // Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai
            u1 = b - c1 / albgri[ib];
            u2 = b - c1 * albgri[ib];
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
            dh7 = (v * du - u * dv) / (v * v);
            u = -c1 * tmp3 * s1;
            du = h * c1 * tmp3 * s1;
            dh8 = (v * du - u * dv) / (v * v);
            v = d2;
            dv = h * tmp4 / s1 + h * tmp5 * s1;
            u = tmp4 / s1;
            du = h * tmp4 / s1;
            dh9 = (v * du - u * dv) / (v * v);
            u = -tmp5 * s1;
            du = h * tmp5 * s1;
            dh10 = (v * du - u * dv) / (v * v);

            da1 = h7 * s2 * s1 + h8 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh7 +
                  (1.0 - s2 / s1) / (twostext - h) * dh8;
            da2 = h9 * s2 * s1 + h10 * s2 / s1 + (1.0 - s2 * s1) / (twostext + h) * dh9 +
                  (1.0 - s2 / s1) / (twostext - h) * dh10;

            // Flux derivatives
            d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1;
            d_fabi = -(dh7 + dh8) - (1.0 - albgri[ib]) * d_ftii;
            d_fabi_sun = (1.0 - omega[ib]) / avmu * (da1 + da2);
            d_fabi_sha = d_fabi - d_fabi_sun;
            fabi_sun_z[iv] = std::max(d_fabi_sun, 0.0);
            fabi_sha_z[iv] = std::max(d_fabi_sha, 0.0);

            // Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
            // to normalize derivatives by sunlit or shaded fraction to get
            // APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
            fabi_sun_z[iv] = fabi_sun_z[iv] / fsun_z[iv];
            fabi_sha_z[iv] = fabi_sha_z[iv] / (1.0 - fsun_z[iv]);
          }
        }
      }
    }
  }
}

void SoilAlbedo() {
  const double albice[numrad] = {0.8, 0.55};

  if (!urbpoi) {
    for (int ib = 0; ib < nband; ib++) {
      if (coszen > 0.0) {
        if (ltype == istsoil || ltype == istcrop) {
          inc = std::max(0.11 - 0.40 * h2osoi_vol[0], 0.0);
          soilcol = isoicol;
          albsod[ib] = std::min(albsat[soilcol][ib] + inc, albdry[soilcol][ib]);
          albsoi[ib] = albsod[ib];
        } else if (ltype == istice || ltype == istice_mec) {
          albsod[ib] = albice[ib] albsoi[ib] = albsod[ib]
        } else if (t_grnd > tfrz || (lakepuddling && ltype == istdlak && t_grnd == tfrz && lake_icefrac[0] < 1.0 &&
                                     lake_icefrac[1] > 0.0)) { // maybe get rid of lake logic?
          albsod[ib] = 0.05 / (std::max(0.001, coszen) + 0.15);
          // This expression is apparently from BATS according to Yongjiu Dai.
          // The diffuse albedo should be an average over the whole sky of an angular-dependent direct expression.
          // The expression above may have been derived to encompass both (e.g. Henderson-Sellers 1986),
          // but I'll assume it applies more appropriately to the direct form for now.
          // ZMS: Attn EK, currently restoring this for wetlands even though it is wrong in order to try to get
          // bfb baseline comparison when no lakes are present. I'm assuming wetlands will be phased out anyway.
          if (litype == istdlak) {
            albsoi[ib] = 0.10;
          } else {
            albsoi[ib] = albsod[ib];
          }

        } else {
          // frozen lake, wetland
          // Introduce crude surface frozen fraction according to D. Mironov (2010)
          // Attn EK: This formulation is probably just as good for "wetlands" if they are not phased out.
          // Tenatively I'm restricting this to lakes because I haven't tested it for wetlands. But if anything
          // the albedo should be lower when melting over frozen ground than a solid frozen lake.
          if (ltype == istdlak && !lakepuddling && snl == 0) {
            // Need to reference snow layers here because t_grnd could be over snow or ice
            // but we really want the ice surface temperature with no snow
            sicefr = 1.0 - exp(-calb * (tfrz - t_grnd) / tfrz);
            albsod[ib] =
                sicefr * alblak[ib] + (1.0 - sicefr) * std::max(alblakwi[ib], 0.05 / (std::max(0.001, coszen) + 0.15));
            albsoi[ib] = sicefr * alblak[ib] + (1.0 - sicefr) * std::max(alblakwi[ib], 0.10);
            // Make sure this is no less than the open water albedo above.
            // Setting lake_melt_icealb(:) = alblak(:) in namelist reverts the melting albedo to the cold
            // snow-free value.
          } else {
            albsod[ib] = alblak[ib];
            albsoi[ib] = albsod[ib];
          }
        }
      }
    }
  }
}

} // namespace ELM
