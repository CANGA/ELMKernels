program test_CanopyHydrology_module
  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4, &
       bool => shr_kind_bool
  implicit none

  character (len = *), dimension(2),parameter :: SURFDAT_FILE_NAME = (/"../links/surfacedataWBW.nc", "../links/surfacedataBRW.nc"/) 
  character (len = 100)  :: FORCDAT_FILE_NAME 
  character (len = *), parameter :: FORCDAT_BASE = "../links/forcing"

  integer(i4) :: frac_veg_nosno = 1 
  integer(i4), parameter :: npfts=17, nmonths=12, ngrcs=24  
  integer(i4), parameter :: ltype=1, ctype=1 
  logical(bool) :: urbpoi = .false.  
  logical(bool) :: do_capsnow = .false. 
  integer(i4) :: n_irrig_steps_left = 0 
  real(r8) :: dtime = 1800. ! time step in seconds 
  real(r8) :: dewmx = 0.1 
  real(r8) :: forc_irrig = 0.0d0 
  real(r8) :: qflx_floodc = 0.0d0 
  real(r8) :: qflx_floodg = 0.0d0 


  real(r8) :: elai, esai, qflx_prec_intr  
  real(r8) :: forc_snow, forc_rain, forc_t, t_grnd 
  real(r8) :: fwet, fdry 
  real(r8) :: qflx_irrig, qflx_prec_grnd, qflx_snwcp_liq, qflx_snwcp_ice
  real(r8) :: qflx_rain_grnd, h2ocan 
  real(r8) :: qflx_snow_grnd_patch, qflx_snow_grnd_col  
  integer(i4) :: itime, p, g

  real(r8), dimension(npfts,ngrcs) :: h2ocan_pft  

  real(r8), dimension(:,:,:), allocatable  :: total_precip
  real(r8), dimension(:,:,:), allocatable  :: tair  

  ! surface grid cell type - contains forcing and phenology variables   
  type surface_cell 
     real(r8), dimension(:), pointer  :: forc_rain  
     real(r8), dimension(:), pointer  :: forc_snow  
     real(r8), dimension(:), pointer  :: forc_t
     real(r8), dimension(:), pointer  :: lai  
     real(r8), dimension(:), pointer  :: sai 
  end type surface_cell

  type(surface_cell), dimension(ngrcs) :: surfdata

  ! snow related variables for kern2 
  ! assumes one column only 
  integer(i4),parameter  :: oldfflag=0 
  integer(i4) :: newnode 
  integer(i4), parameter :: nlevsno=5 
  integer(i4) :: snl 
  real(r8), dimension(-nlevsno+1:0) :: swe_old=0., h2osoi_liq=0., h2osoi_ice=0. 
  real(r8), dimension(-nlevsno+1:0) :: dz=0., zi=0., z=0. 
  real(r8), dimension(-nlevsno+1:0) :: t_soisno, frac_iceold 
  real(r8) :: frac_sno, frac_sno_eff, int_snow, h2osno 
  real(r8) :: qflx_snow_h2osfc , snow_depth
  real(r8) :: frac_h2osfc
  real(r8) :: qflx_snow_melt = 0.0d0 ! think about updating this 
  real(r8) :: n_melt = 0.7 

  integer(i4), dimension(ngrcs) :: snl_grc = 0
  real(r8), dimension(ngrcs)  :: snow_depth_grc=0.0d0 
  real(r8), dimension(ngrcs) :: h2osno_grc=0.0d0 
  real(r8), dimension(ngrcs,-nlevsno+1:0) :: h2osoi_liq_grc=0.0d0, h2osoi_ice_grc=0.0d0 
  real(r8), dimension(ngrcs,-nlevsno+1:0) :: dz_grc=0.0d0 

  ! variables for fracH2Osfc  
  real(r8) :: micro_sigma = 0.1
  real(r8) :: min_h2osfc = 1.0d-8 
  real(r8) :: h2osfc = 0.1 !arbitrary 
  real(r8) :: qflx_h2osfc2topsoi


  ! load lai and sai 
  call get_phenology(SURFDAT_FILE_NAME(1), surfdata(1:12), npfts, nmonths )    ! walker branch watershed 
  call get_phenology(SURFDAT_FILE_NAME(2), surfdata(12+1:24), npfts, nmonths ) ! barrow AK  

  ! now get forcing data 
  do g=1,ngrcs 

     write(FORCDAT_FILE_NAME,"(a,i0,a)") FORCDAT_BASE,g,".nc"
     FORCDAT_FILE_NAME = trim( adjustl(FORCDAT_FILE_NAME) ) 

     call get_forcing_data(FORCDAT_FILE_NAME, surfdata(g)) 

  end do !month loop 


  ! end setup 


  h2ocan_pft = 0.0d0 
  frac_h2osfc = 0.0d0
  print*, "Time", "Total Canopy Water", "Min Water", "Max Water", "Total Snow", &
  "Min Snow", "Max Snow", "Avg Frac Sfc", "Min Frac Sfc", "Max Frac Sfc"
  print*, 0, sum(h2ocan_pft), minval(h2ocan_pft), maxval(h2ocan_pft), sum(h2osno_grc), &
  minval(h2osno_grc), maxval(h2osno_grc), sum(h2osoi_liq), minval(h2osoi_liq), maxval(h2osoi_liq)

  do itime=1,28*48  ! February is shortest month 

     do g=1,ngrcs ! grid cell loop 
        forc_rain = surfdata(g)%forc_rain(itime) 
        forc_snow = surfdata(g)%forc_snow(itime) 
        forc_t = surfdata(g)%forc_t(itime) 
        t_grnd = forc_t 


        qflx_snow_grnd_col=0.0d0 
        do p = 1, npfts !pft loop  

           elai = surfdata(g)%lai(p)  
           esai = surfdata(g)%sai(p)  
           h2ocan = h2ocan_pft(p,g) 

           call CanopyHydrology_Interception( dtime, &
                forc_rain, forc_snow, forc_irrig, &
                ltype, ctype, urbpoi, do_capsnow, &
                elai, esai, dewmx, frac_veg_nosno, &
                h2ocan, n_irrig_steps_left, &
                qflx_prec_intr, qflx_irrig, qflx_prec_grnd, &
                qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd_patch, qflx_rain_grnd)

           call CanopyHydrology_FracWet(frac_veg_nosno, h2ocan, elai, esai, dewmx, fwet, fdry)

           h2ocan_pft(p,g) = h2ocan
           qflx_snow_grnd_col = qflx_snow_grnd_col + qflx_snow_grnd_patch 

        end do

        snl =  snl_grc(g) 
        snow_depth = snow_depth_grc(g) 
        h2osoi_liq(:) = h2osoi_liq_grc(g,:) 
        h2osoi_ice(:) = h2osoi_ice_grc(g,:) 
        dz(:) = dz_grc(g,:) 
        h2osno = h2osno_grc(g)
        int_snow = 0.
        
        call CanopyHydrology_SnowWater( dtime, &
             qflx_floodg, &
             ltype, ctype, urbpoi, do_capsnow, oldfflag, &
             forc_t, t_grnd, qflx_snow_grnd_col, qflx_snow_melt, n_melt, frac_h2osfc, & ! forcing 
             snow_depth, h2osno, int_snow, swe_old, &
             h2osoi_liq, h2osoi_ice, t_soisno, frac_iceold, & ! state
             snl, dz, z, zi, newnode, & ! snow mesh for initialization
             qflx_floodc, qflx_snow_h2osfc, &
             frac_sno, frac_sno_eff)

        ! FIXME: Fortran black magic... h2osoi_liq is a vector, but
        ! the interface specifies a single double.  --etc
        call CanopyHydrology_FracH2OSfc( dtime, min_h2osfc, &
             ltype, micro_sigma, h2osno, &
             h2osfc, h2osoi_liq, frac_sno, frac_sno_eff, &
             qflx_h2osfc2topsoi, frac_h2osfc)

        snow_depth_grc(g) = snow_depth  
        snl_grc(g) = snl 
        h2osno_grc(g) = h2osno 

        h2osoi_liq_grc(g,:) = h2osoi_liq 
        h2osoi_ice_grc(g,:) = h2osoi_ice 
        dz_grc(g,:) = dz 

     end do !mesh loop 
     print*, itime, sum(h2ocan_pft), minval(h2ocan_pft), maxval(h2ocan_pft), &
     sum(h2osno_grc), minval(h2osno_grc), maxval(h2osno_grc), &
     sum(h2osoi_liq), minval(h2osoi_liq), maxval(h2osoi_liq)

  end do ! time loop  


  stop 
contains 

  subroutine err_handle( status )
    use netcdf
    integer :: status 
    print *, 'error ', status, nf90_strerror(status)
    stop 
  end subroutine err_handle

  subroutine get_forcing_data( filename , asurfcell)
    use netcdf
    character(len=*), intent(in) :: filename 
    type(surface_cell), intent(out) :: asurfcell
    ! local variables 
    integer(i4), dimension(3) :: start3, count3 
    integer(i4) :: ncid, varid, status , dimid 
    integer(i4) :: ntimes 

    start3=(/1,1,1/) 


    status = nf90_open(filename, NF90_NOWRITE, ncid) 
    if( status .ne. 0) call err_handle( status ) 
    status = nf90_inq_dimid(ncid, "time", dimid) 
    if( status .ne. 0) call err_handle( status ) 
    status = nf90_inquire_dimension(ncid, dimid, len=ntimes) 
    if( status .ne. 0) call err_handle( status ) 
    count3=(/1,1,ntimes/) 

    allocate(asurfcell%forc_rain(ntimes) ) 
    allocate(asurfcell%forc_snow(ntimes) ) 
    allocate(asurfcell%forc_t(ntimes) ) 

    ! total precip 
    status = nf90_inq_varid(ncid, "PRECTmms", varid) 
    if( status .ne. 0) call err_handle( status ) 
    allocate( total_precip(1,1,ntimes) ) 
    status = nf90_get_var(ncid, varid, total_precip, start=start3, count=count3) 
    if( status .ne. 0) call err_handle( status ) 

    ! air temperature 
    status = nf90_inq_varid(ncid, "TBOT", varid) 
    if( status .ne. 0) call err_handle( status ) 
    allocate( tair(1,1,ntimes) ) 
    status = nf90_get_var(ncid, varid, tair, start=start3, count=count3) 
    if( status .ne. 0) call err_handle( status ) 

    do itime=1,ntimes 
       if( tair(1,1,itime) .ge. 273.15) then 
          asurfcell%forc_rain(itime) = total_precip(1,1,itime) 
          asurfcell%forc_snow(itime) = 0.0 
       else 
          asurfcell%forc_rain(itime) = 0.0 
          asurfcell%forc_snow(itime) = total_precip(1,1,itime) 
       end if
       asurfcell%forc_t(itime) = tair(1,1,itime) 
    end do
    deallocate(total_precip) 
    deallocate(tair) 

  end subroutine get_forcing_data

  subroutine get_phenology( filename, smesh, npfts,nmonths) 
    use netcdf

    character(len=*), intent(in) :: filename 
    type(surface_cell), dimension(nmonths) :: smesh 
    integer(i4), intent(in) :: npfts, nmonths 
    real(r8), dimension(npfts,nmonths) :: lai
    real(r8), dimension(npfts,nmonths) :: sai
    integer(i4), dimension(4) ::  count
    integer(i4), dimension(4) ::  start=(/1,1,1,1/)    
    integer(i4) :: ncid, varid, status

    count=(/1,1,npfts,nmonths/)  ! netcdf files are row major damn it 

    ! open surface data file 
    status = nf90_open(filename, NF90_NOWRITE, ncid) 
    if( status .ne. 0) call err_handle( status ) 

    ! get monthly leaf area index and stem area index 

    status = nf90_inq_varid(ncid, "MONTHLY_LAI", varid) 
    status = nf90_get_var(ncid, varid, lai,start=start,count=count) 
    if( status .ne. 0) call err_handle( status ) 

    status = nf90_inq_varid(ncid, "MONTHLY_SAI", varid) 
    status = nf90_get_var(ncid, varid, sai,start=start,count=count) 
    if( status .ne. 0) call err_handle( status ) 

    do g=1,nmonths 
       allocate(smesh(g)%lai(npfts) ) 
       allocate(smesh(g)%sai(npfts) ) 
       smesh(g)%lai(:) = lai(:,g) 
       smesh(g)%sai(:) = sai(:,g) 
    end do

    status = nf90_close(ncid) 
    if( status .ne. 0) call err_handle( status ) 

  end subroutine get_phenology



end program test_CanopyHydrology_module
