subroutine CanopyHydrologyKern1( dtime, &  
     forc_rain, forc_snow, irrig_rate, &
     ltype, ctype, urbpoi, do_capsnow, &
     elai, esai, dewmx, frac_veg_nosno, &
     h2ocan, n_irrig_steps_left, &
     qflx_prec_intr, qflx_irrig, qflx_prec_grnd, &
     qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd_patch, qflx_rain_grnd)

  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4, &
       bool => shr_kind_bool

  use column_varcon      , only : icol_sunwall, icol_shadewall
  use landunit_varcon    , only : istcrop, istice, istwet, istsoil, istice_mec

  implicit none 

  real(r8), intent(in)  :: dtime 
  integer(i4), intent(in)   :: ltype , ctype
  integer(i4), intent(inout) ::  n_irrig_steps_left
  logical(bool), intent(in)   :: urbpoi, do_capsnow 
  real(r8), intent(out) :: qflx_prec_intr 
  integer(i4), intent(in)   :: frac_veg_nosno 
  real(r8), intent(in)  :: forc_rain 
  real(r8), intent(in)  :: forc_snow 
  real(r8), intent(in)  :: dewmx 
  real(r8), intent(in)  :: elai
  real(r8), intent(in)  :: esai 
  real(r8), intent(inout) :: h2ocan
  real(r8), intent(out) :: qflx_irrig 
  real(r8), intent(out)    :: qflx_prec_grnd
  real(r8), intent(out) :: qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd_patch, qflx_rain_grnd
  real(r8), intent(in)  ::   irrig_rate

  !local variables  
  real(r8) :: fpi, xrun, h2ocanmx   
  real(r8) :: qflx_candrip, qflx_through_snow, qflx_through_rain 
  real(r8) :: qflx_prec_grnd_snow
  real(r8) :: qflx_prec_grnd_rain 
  real(r8) :: fracsnow 
  real(r8) :: fracrain 


  ! Canopy interception and precipitation onto ground surface
  ! Add precipitation to leaf water

  if (ltype==istsoil .or. ltype==istwet .or. urbpoi .or. &
       ltype==istcrop) then

     qflx_candrip = 0._r8      ! rate of canopy runoff
     qflx_through_snow = 0._r8 ! rain precipitation direct through canopy
     qflx_through_rain = 0._r8 ! snow precipitation direct through canopy
     qflx_prec_intr = 0._r8    ! total intercepted precipitation
     fracsnow = 0._r8          ! fraction of input precip that is snow
     fracrain = 0._r8          ! fraction of input precip that is rain


     if (ctype /= icol_sunwall .and. ctype /= icol_shadewall) then
        if (frac_veg_nosno == 1 .and. (forc_rain + forc_snow) > 0._r8) then
           ! determine fraction of input precipitation that is snow and rain
           fracsnow = forc_snow/(forc_snow + forc_rain)
           fracrain = forc_rain/(forc_snow + forc_rain)

           ! The leaf water capacities for solid and liquid are different,
           ! generally double for snow, but these are of somewhat less
           ! significance for the water budget because of lower evap. rate at
           ! lower temperature.  Hence, it is reasonable to assume that
           ! vegetation storage of solid water is the same as liquid water.
           h2ocanmx = dewmx * (elai + esai)

           ! Coefficient of interception
           ! set fraction of potential interception to max 0.25
           fpi = 0.25_r8*(1._r8 - exp(-0.5_r8*(elai + esai)))

           ! Direct throughfall
           qflx_through_snow = forc_snow * (1._r8-fpi)
           qflx_through_rain = forc_rain * (1._r8-fpi)

           ! Intercepted precipitation [mm/s]
           qflx_prec_intr = (forc_snow + forc_rain) * fpi

           ! Water storage of intercepted precipitation and dew
           h2ocan = max(0._r8, h2ocan + dtime*qflx_prec_intr)

           ! Initialize rate of canopy runoff and snow falling off canopy
           qflx_candrip = 0._r8

           ! Excess water that exceeds the leaf capacity
           xrun = (h2ocan - h2ocanmx)/dtime

           ! Test on maximum dew on leaf
           ! Note if xrun > 0 then h2ocan must be at least h2ocanmx
           if (xrun > 0._r8) then
              qflx_candrip = xrun
              h2ocan = h2ocanmx
           end if

        end if
     end if

  else if (ltype==istice .or. ltype==istice_mec) then
     
     h2ocan            = 0._r8
     qflx_candrip      = 0._r8
     qflx_through_snow = 0._r8
     qflx_through_rain = 0._r8
     qflx_prec_intr    = 0._r8
     fracsnow          = 0._r8
     fracrain          = 0._r8

  end if

  ! Precipitation onto ground (kg/(m2 s))

  if (ctype /= icol_sunwall .and. ctype /= icol_shadewall) then
     if (frac_veg_nosno == 0) then
        qflx_prec_grnd_snow = forc_snow
        qflx_prec_grnd_rain = forc_rain
     else
        qflx_prec_grnd_snow = qflx_through_snow + (qflx_candrip * fracsnow)
        qflx_prec_grnd_rain = qflx_through_rain + (qflx_candrip * fracrain)
     end if
     ! Urban sunwall and shadewall have no intercepted precipitation
  else
     qflx_prec_grnd_snow = 0.
     qflx_prec_grnd_rain = 0.
  end if

  ! Determine whether we're irrigating here; set qflx_irrig appropriately
  if (n_irrig_steps_left > 0) then
     qflx_irrig         = irrig_rate
     n_irrig_steps_left = n_irrig_steps_left - 1
  else
     qflx_irrig = 0._r8
  end if

  ! Add irrigation water directly onto ground (bypassing canopy interception)
  ! Note that it's still possible that (some of) this irrigation water will runoff (as runoff is computed later)
  qflx_prec_grnd_rain = qflx_prec_grnd_rain + qflx_irrig

  ! Done irrigation

  qflx_prec_grnd = qflx_prec_grnd_snow + qflx_prec_grnd_rain

  if (do_capsnow) then
     qflx_snwcp_liq = qflx_prec_grnd_rain
     qflx_snwcp_ice = qflx_prec_grnd_snow

     qflx_snow_grnd_patch = 0._r8
     qflx_rain_grnd = 0._r8
  else

     qflx_snwcp_liq = 0._r8
     qflx_snwcp_ice = 0._r8
     qflx_snow_grnd_patch = qflx_prec_grnd_snow           ! ice onto ground (mm/s)
     qflx_rain_grnd     = qflx_prec_grnd_rain           ! liquid water onto ground (mm/s)
  end if

  return 

end subroutine CanopyHydrologyKern1

