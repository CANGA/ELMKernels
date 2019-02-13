subroutine CanopyHydrology_FracH2OSfc( dtime, min_h2osfc, &
     itype, micro_sigma, h2osno, &
     h2osfc, h2osoi_liq, frac_sno, frac_sno_eff, &
     qflx_h2osfc2topsoi, frac_h2osfc, &
     no_update)

  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4
  use landunit_varcon    , only : istcrop, istice, istwet, istsoil, istice_mec
  implicit none 

  real(r8), intent(in) :: dtime, min_h2osfc, micro_sigma, h2osno
  integer(i4), intent(in) :: itype
  integer(i4), intent(in), optional :: no_update

  real(r8), intent(inout) :: h2osfc, h2osoi_liq, frac_sno, frac_sno_eff
  real(r8), intent(out) :: qflx_h2osfc2topsoi, frac_h2osfc

  !local
  real(r8), parameter :: shr_const_pi=4.0d0*atan(1.0d0)

  integer :: l
  real(r8):: d,fd,dfdd      ! temporary variable for frac_h2oscs iteration
  real(r8):: sigma          ! microtopography pdf sigma in mm

  qflx_h2osfc2topsoi = 0._r8
  ! h2osfc only calculated for soil vegetated land units
  if ( itype  == istsoil .or. itype == istcrop) then

     !  Use newton-raphson method to iteratively determine frac_h20sfc
     !  based on amount of surface water storage (h2osfc) and
     !  microtopography variability (micro_sigma)

     if (h2osfc > min_h2osfc) then
        ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
        d=0.0

        sigma=1.0e3 * micro_sigma ! convert to mm
        do l=1,10
           fd = 0.5*d*(1.0_r8+erf(d/(sigma*sqrt(2.0)))) &
                +sigma/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*sigma**2)) &
                -h2osfc
           dfdd = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

           d = d - fd/dfdd
        enddo
        !--  update the submerged areal fraction using the new d value
        frac_h2osfc = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

     else
        frac_h2osfc = 0._r8
        h2osoi_liq = h2osoi_liq + h2osfc
        qflx_h2osfc2topsoi = h2osfc/dtime
        h2osfc=0._r8
     endif

     if (.not. present(no_update)) then

        ! adjust fh2o, fsno when sum is greater than zero
        if (frac_sno > (1._r8 - frac_h2osfc) .and. h2osno > 0) then

           if (frac_h2osfc > 0.01_r8) then
              frac_h2osfc = max(1.0_r8 - frac_sno,0.01_r8)
              frac_sno = 1.0_r8 - frac_h2osfc
           else
              frac_sno = 1.0_r8 - frac_h2osfc
           endif
           frac_sno_eff=frac_sno

        endif

     endif ! end of no_update construct

  else !if landunit not istsoil/istcrop, set frac_h2osfc to zero

     frac_h2osfc = 0._r8

  endif

  return
end subroutine CanopyHydrology_FracH2OSfc
