subroutine CanopyHydrologyKern2( dtime, oldfflag, newnode, &  
     ctype, qflx_floodc, qflx_snow_h2osfc, snow_depth, snl, & 
     swe_old, h2osoi_liq, h2osoi_ice, dz, z, zi, & 
     t_soisno, frac_iceold, &  
     do_capsnow, frac_h2osfc, qflx_snow_grnd_col, frac_sno, int_snow, forc_t, & 
     h2osno,qflx_snow_melt, n_melt, frac_sno_eff, t_grnd, & 
     qflx_floodg,  &  
     ltype, urbpoi) bind(C)

  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4, &
       bool => shr_kind_bool
  use column_varcon      , only : icol_sunwall, icol_shadewall
  use landunit_varcon    , only : istcrop, istice, istwet, istsoil, istice_mec
  use clm_varpar         , only : nlevsno
  !use column_varcon      , only : icol_sunwall, icol_shadewall
  !use landunit_varcon    , only : istcrop, istice, istwet, istsoil, istice_mec 
  !use clm_varcon         , only : zlnd, rpi,  tfrz
  use clm_varctl         , only : subgridflag

  
  implicit none 

  real(r8), parameter :: rpi=4.0d0*atan(1.0d0)  
  real(r8), parameter :: tfrz=273.15
  real(r8) :: zlnd = 0.01_r8


  real(r8), intent(in)    :: dtime 
  integer, intent(in)     :: oldfflag, ctype , ltype
  integer, intent(inout)  :: snl 
  integer, intent(out)  :: newnode 
  real(r8), intent(out) :: qflx_floodc 
  real(r8), intent(in)  :: qflx_floodg 
  real(r8), intent(out)  :: qflx_snow_h2osfc 
  real(r8), intent(out), dimension(-nlevsno+1:0)  :: swe_old 
  real(r8), intent(inout) :: snow_depth , h2osno 
  real(r8), intent(inout), dimension(-nlevsno+1:0) :: h2osoi_liq, h2osoi_ice
  real(r8), intent(inout), dimension(-nlevsno+1:0)  :: dz
  real(r8), intent(out), dimension(-nlevsno+1:0)  :: z, zi
  real(r8), intent(out), dimension(-nlevsno+1:0)  :: t_soisno, frac_iceold
  logical, intent(in) :: do_capsnow, urbpoi 
  real(r8), intent(in) :: frac_h2osfc, qflx_snow_grnd_col, forc_t , qflx_snow_melt, n_melt
  real(r8), intent(out) :: int_snow, frac_sno_eff, frac_sno
  real(r8), intent(in) :: t_grnd

  ! local variables 
  real(r8) :: temp_intsnow, temp_snow_depth, z_avg, fmelt, dz_snowf, snowmelt
  real(r8) :: newsnow, bifall, fsnow_new, accum_factor, fsno_new, smr 
  integer :: j 

  ! apply gridcell flood water flux to non-lake columns
  if (ctype /= icol_sunwall .and. ctype /= icol_shadewall) then      
     qflx_floodc = qflx_floodg
  else
     qflx_floodc = 0._r8
  endif

  ! Determine snow height and snow water

  ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
  ! U.S.Department of Agriculture Forest Service, Project F,
  ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

  qflx_snow_h2osfc = 0._r8
  ! set temporary variables prior to updating
  temp_snow_depth=snow_depth
  ! save initial snow content
  do j= -nlevsno+1,snl
     swe_old(j) = 0.0_r8
  end do
  do j= snl+1,0
     swe_old(j)=h2osoi_liq(j)+h2osoi_ice(j)
  enddo

  if (do_capsnow) then
     dz_snowf = 0._r8
     newsnow = (1._r8 - frac_h2osfc) * qflx_snow_grnd_col * dtime
     frac_sno=1._r8
     int_snow = 5.e2_r8
  else
     if (forc_t > tfrz + 2._r8) then
        bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
     else if (forc_t > tfrz - 15._r8) then
        bifall=50._r8 + 1.7_r8*(forc_t - tfrz + 15._r8)**1.5_r8
     else
        bifall=50._r8
     end if

     ! newsnow is all snow that doesn't fall on h2osfc
     newsnow = (1._r8 - frac_h2osfc) * qflx_snow_grnd_col * dtime

     ! update int_snow
     int_snow = max(int_snow,h2osno) !h2osno could be larger due to frost

     ! snowmelt from previous time step * dtime
     snowmelt = qflx_snow_melt * dtime

     ! set shape factor for accumulation of snow
     accum_factor=0.1

     if (h2osno > 0.0) then

        !======================  FSCA PARAMETERIZATIONS  ======================
        ! fsca parameterization based on *changes* in swe
        ! first compute change from melt during previous time step
        if(snowmelt > 0._r8) then

           smr=min(1._r8,(h2osno)/(int_snow))

           frac_sno = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(n_melt)

        endif

        ! update fsca by new snow event, add to previous fsca
        if (newsnow > 0._r8) then
           fsno_new = 1._r8 - (1._r8 - tanh(accum_factor*newsnow))*(1._r8 - frac_sno)
           frac_sno = fsno_new

           ! reset int_snow after accumulation events
           temp_intsnow= (h2osno + newsnow) &
                / (0.5*(cos(rpi*(1._r8-max(frac_sno,1e-6_r8))**(1./n_melt))+1._r8))
           int_snow = min(1.e8_r8,temp_intsnow)
        endif

        !====================================================================

        ! for subgrid fluxes
        if (subgridflag ==1 .and. .not. urbpoi) then
           if (frac_sno > 0._r8)then
              snow_depth=snow_depth + newsnow/(bifall * frac_sno)
           else
              snow_depth=0._r8
           end if
        else
           ! for uniform snow cover
           snow_depth=snow_depth+newsnow/bifall
        endif

        ! use original fsca formulation (n&y 07)
        if (oldfflag == 1) then 
           ! snow cover fraction in Niu et al. 2007
           if(snow_depth > 0.0_r8)  then
              frac_sno = tanh(snow_depth/(2.5_r8*zlnd* &
                   (min(800._r8,(h2osno+ newsnow)/snow_depth)/100._r8)**1._r8) )
           endif
           if(h2osno < 1.0_r8)  then
              frac_sno=min(frac_sno,h2osno)
           endif
        endif

     else !h2osno == 0
        ! initialize frac_sno and snow_depth when no snow present initially
        if (newsnow > 0._r8) then 
           z_avg = newsnow/bifall
           fmelt=newsnow
           frac_sno = tanh(accum_factor*newsnow)

           ! make int_snow consistent w/ new fsno, h2osno
           int_snow = 0. !reset prior to adding newsnow below
           temp_intsnow= (h2osno + newsnow) &
                / (0.5*(cos(rpi*(1._r8-max(frac_sno,1e-6_r8))**(1./n_melt))+1._r8))
           int_snow = min(1.e8_r8,temp_intsnow)

           ! update snow_depth and h2osno to be consistent with frac_sno, z_avg
           if (subgridflag ==1 .and. .not. urbpoi) then
              snow_depth=z_avg/frac_sno
           else
              snow_depth=newsnow/bifall
           endif
           ! use n&y07 formulation
           if (oldfflag == 1) then 
              ! snow cover fraction in Niu et al. 2007
              if(snow_depth > 0.0_r8)  then
                 frac_sno = tanh(snow_depth/(2.5_r8*zlnd* &
                      (min(800._r8,newsnow/snow_depth)/100._r8)**1._r8) )
              endif
           endif
        else
           z_avg = 0._r8
           snow_depth = 0._r8
           frac_sno = 0._r8
        endif
     endif ! end of h2osno > 0

     ! snow directly falling on surface water melts, increases h2osfc
     qflx_snow_h2osfc = frac_h2osfc*qflx_snow_grnd_col

     ! update h2osno for new snow
     h2osno = h2osno + newsnow 
     int_snow = int_snow + newsnow

     ! update change in snow depth
     dz_snowf = (snow_depth - temp_snow_depth) / dtime

  end if !end of do_capsnow construct

  ! set frac_sno_eff variable
  if (ltype == istsoil .or. ltype == istcrop) then
     if (subgridflag ==1) then 
        frac_sno_eff = frac_sno
     else
        frac_sno_eff = 1._r8
     endif
  else
     frac_sno_eff = 1._r8
  endif

  if (ltype==istwet .and. t_grnd>tfrz) then
     h2osno=0._r8
     snow_depth=0._r8
  end if

  ! When the snow accumulation exceeds 10 mm, initialize snow layer
  ! Currently, the water temperature for the precipitation is simply set
  ! as the surface air temperature
  newnode = 0    ! flag for when snow node will be initialized
  if (snl == 0 .and. qflx_snow_grnd_col > 0.0_r8 .and. frac_sno*snow_depth >= 0.01_r8) then
     newnode = 1
     snl = -1
     dz(0) = snow_depth                       ! meter
     z(0) = -0.5_r8*dz(0)
     zi(-1) = -dz(0)
     t_soisno(0) = min(tfrz, forc_t)      ! K
     h2osoi_ice(0) = h2osno               ! kg/m2
     h2osoi_liq(0) = 0._r8                   ! kg/m2
     frac_iceold(0) = 1._r8
  end if

  ! The change of ice partial density of surface node due to precipitation.
  ! Only ice part of snowfall is added here, the liquid part will be added
  ! later.
  if (snl < 0 .and. newnode == 0) then
     h2osoi_ice(snl+1) = h2osoi_ice(snl+1)+newsnow
     dz(snl+1) = dz(snl+1)+dz_snowf*dtime
  end if

end subroutine CanopyHydrologyKern2

