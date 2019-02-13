program CanopyHydrology_kern1_single
  use netcdf
  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4, &
       bool => shr_kind_bool

  implicit none

  character (len = *), parameter :: SURFDAT_FILE_NAME = "../links/surfacedataWBW.nc"
  character (len = *), parameter :: FORCDAT_FILE_NAME = "../links/forcing7.nc"

  integer(i4), parameter :: npfts=17, nmonths=12
  integer(i4) :: ltype=1, ctype=1 
  logical(bool) :: urbpoi, do_capsnow 
  integer(i4) :: n_irrig_steps_left = 0 
  integer(i4) :: frac_veg_nosno = 1 
  real(r8) :: elai, esai, dewmx, forc_snow, forc_rain, qflx_prec_intr  
  real(r8) :: qflx_irrig, qflx_prec_grnd, qflx_snwcp_liq, qflx_snwcp_ice
  real(r8) :: qflx_snow_grnd_patch, qflx_rain_grnd,irrig_rate, h2ocan 
  integer(i4), dimension(4) :: start, count 
  integer(i4), dimension(3) :: start3, count3 
  integer(i4) :: itime,ntimes 
  real(r8), dimension(npfts,nmonths) :: monthly_lai 
  real(r8), dimension(npfts,nmonths) :: monthly_sai 
  real(r8), dimension(:,:,:), allocatable  :: total_precip
  real(r8) :: dtime = 1800.

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid, dimid, status

  status = nf90_open(SURFDAT_FILE_NAME, NF90_NOWRITE, ncid) 
  if( status .ne. 0) call err_handle( status ) 
   

  status = nf90_inq_varid(ncid, "MONTHLY_LAI", varid) 
  count=(/1,1,npfts,nmonths/)  ! netcdf files are row major damn it 
  start=(/1,1,1,1/)    ! 
  status = nf90_get_var(ncid, varid, monthly_lai,start=start,count=count) 
  if( status .ne. 0) call err_handle( status ) 
  !do ipft=1,npfts 
  ! print *, monthly_lai(ipft,12), monthly_lai(ipft,3), monthly_lai(ipft,6) 
  !end do 

  status = nf90_inq_varid(ncid, "MONTHLY_SAI", varid) 
  status = nf90_get_var(ncid, varid, monthly_sai,start=start,count=count) 
  if( status .ne. 0) call err_handle( status ) 

  status = nf90_close(ncid) 

  status = nf90_open(FORCDAT_FILE_NAME, NF90_NOWRITE, ncid) 
  if( status .ne. 0) call err_handle( status ) 
  status = nf90_inq_dimid(ncid, "time", dimid) 
  if( status .ne. 0) call err_handle( status ) 
  status = nf90_inquire_dimension(ncid, dimid, len=ntimes) 
  if( status .ne. 0) call err_handle( status ) 
  status = nf90_inq_varid(ncid, "PRECTmms", varid) 
  if( status .ne. 0) call err_handle( status ) 
  start3=(/1,1,1/) 
  count3=(/1,1,ntimes/) 
  allocate( total_precip(1,1,ntimes) ) 
  status = nf90_get_var(ncid, varid, total_precip, start=start3, count=count3) 
  if( status .ne. 0) call err_handle( status ) 
  status = nf90_close(ncid) 

  h2ocan = 0.0d0 
  print *, "Timestep, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr"
  do itime=1,ntimes 


   forc_rain = total_precip(1,1,itime) 
   forc_snow = 0.d0 
   dewmx = 0.1 
   elai = monthly_lai(8,6)  
   esai = monthly_sai(8,6)  
   frac_veg_nosno = 1 
   irrig_rate = 0.0d0 
   n_irrig_steps_left = 0
   urbpoi=.false. 
   do_capsnow = .false. 
   call CanopyHydrology_Interception( dtime, &
     forc_rain, forc_snow, irrig_rate, &
     ltype, ctype, urbpoi, do_capsnow, &
     elai, esai, dewmx, frac_veg_nosno, &
     h2ocan, n_irrig_steps_left, &
     qflx_prec_intr, qflx_irrig, qflx_prec_grnd, &
     qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd_patch, qflx_rain_grnd)

   print *, itime, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr 
  end do 


  stop 
  contains 

    subroutine err_handle( status ) 
    integer :: status 
   
    print *, 'error ', status 

    stop 

    end subroutine err_handle 

end program 
