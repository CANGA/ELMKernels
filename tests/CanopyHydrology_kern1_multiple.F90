program CanopyHydrology_kern1_multiple
  use netcdf 
  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4, &
       bool => shr_kind_bool
  implicit none

  character (len = *), dimension(2),parameter :: SURFDAT_FILE_NAME = (/"links/surfacedataWBW.nc", "links/surfacedataBRW.nc"/) 
  character (len = 100)  :: FORCDAT_FILE_NAME 
  character (len = *), parameter :: FORCDAT_BASE = "links/forcing"

  integer(i4) :: frac_veg_nosno = 1 
  integer(i4), parameter :: npfts=17, nmonths=12, ngrcs=24  
  integer(i4), parameter :: ltype=1, ctype=1 
  logical(bool) :: urbpoi = .false.  
  logical(bool) :: do_capsnow = .false. 
  integer(i4) :: n_irrig_steps_left = 0 
  real(r8) :: dtime = 1800. ! time step in seconds 
  real(r8) :: dewmx = 0.1 
  real(r8) :: irrig_rate = 0.0d0 

  real(r8) :: elai, esai, qflx_prec_intr  
  real(r8) :: forc_snow, forc_rain
  real(r8) :: qflx_irrig, qflx_prec_grnd, qflx_snwcp_liq, qflx_snwcp_ice
  real(r8) :: qflx_rain_grnd, h2ocan 
  real(r8) :: qflx_snow_grnd_patch, qflx_snow_grnd_col  
  integer :: itime, p, g
  real(r8) :: forc_t

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

  ! load lai and sai 
  call get_phenology(SURFDAT_FILE_NAME(1), surfdata(1:12), npfts, nmonths )    ! walker branch watershed 
  call get_phenology(SURFDAT_FILE_NAME(2), surfdata(12+1:24), npfts, nmonths ) ! barrow AK  

  ! now get forcing data 
  do g=1,ngrcs 
     write(FORCDAT_FILE_NAME,"(a,i0,a)") FORCDAT_BASE,g,".nc"
     FORCDAT_FILE_NAME = trim( adjustl(FORCDAT_FILE_NAME) ) 

     call get_forcing_data(FORCDAT_FILE_NAME, surfdata(g)) 
  end do !month loop 

  ! ------------------------------
  ! end setup stage
  ! -----------------------------
  h2ocan_pft = 0.0d0 
  print*, "Time", "Total Canopy Water", "Min Water", "MaxWater"
  print*, 0, sum(h2ocan_pft), minval(h2ocan_pft), maxval(h2ocan_pft)
  
  do itime=1,28*48  ! February is shortest month 

     do g=1,ngrcs ! grid cell loop 
        forc_rain = surfdata(g)%forc_rain(itime) 
        forc_snow = surfdata(g)%forc_snow(itime) 
        forc_t = surfdata(g)%forc_t(itime) 

        qflx_snow_grnd_col=0.0d0 
        do p = 1, npfts !pft loop  

           elai = surfdata(g)%lai(p)  
           esai = surfdata(g)%sai(p)  
           h2ocan = h2ocan_pft(p,g) 

           call CanopyHydrologyKern1( dtime, &
                forc_rain, forc_snow, irrig_rate, &
                ltype, ctype, urbpoi, do_capsnow, &
                elai, esai, dewmx, frac_veg_nosno, &
                h2ocan, n_irrig_steps_left, &
                qflx_prec_intr, qflx_irrig, qflx_prec_grnd, &
                qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd_patch, qflx_rain_grnd)

           !print*, g, p, forc_rain, forc_snow, elai, esai, h2ocan, qflx_prec_intr
           
           h2ocan_pft(p,g) = h2ocan
           qflx_snow_grnd_col = qflx_snow_grnd_col + qflx_snow_grnd_patch 
        end do ! PFT loop

     end do ! grid cell loop

     print*, itime, sum(h2ocan_pft), minval(h2ocan_pft), maxval(h2ocan_pft)
     
  end do ! time loop  

  ! print*, ""
  ! print*, "Final canopy water:"
  ! print*, "Grid Cell", "PFT", "H2O"
  ! do g=1,ngrcs
  !    do p=1,npfts
  !       print*, g, p, h2ocan_pft(p,g)
  !    end do
  ! end do  
  stop 
contains 

  subroutine err_handle( status ) 
    integer :: status 
    print *, 'error ', status 
    stop 
  end subroutine err_handle

  subroutine get_forcing_data( filename , asurfcell)  
    character(len=*), intent(in) :: filename 
    type(surface_cell), intent(out) :: asurfcell
    ! local variables 
    integer, dimension(3) :: start3, count3 
    integer :: ncid, varid, status , dimid 
    integer :: ntimes 

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

    character(len=*), intent(in) :: filename 
    type(surface_cell), dimension(nmonths) :: smesh 
    integer, intent(in) :: npfts, nmonths 
    real(r8), dimension(npfts,nmonths) :: lai
    real(r8), dimension(npfts,nmonths) :: sai
    integer, dimension(4) ::  count
    integer, dimension(4) ::  start=(/1,1,1,1/)    
    integer :: ncid, varid, status

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



end program CanopyHydrology_kern1_multiple
