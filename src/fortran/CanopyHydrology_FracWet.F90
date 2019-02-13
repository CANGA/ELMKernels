subroutine CanopyHydrology_FracWet( frac_veg_nosno, h2ocan, elai, esai, dewmx, fwet, fdry )
  use shr_kind_mod, only : &
       r8 => shr_kind_r8, &
       i4 => shr_kind_i4
  implicit none 


  integer(i4), intent(in)  :: frac_veg_nosno
  real(r8), intent(in) :: h2ocan, elai, esai, dewmx
  real(r8), intent(out) :: fwet, fdry

  !LOCAL VARIABLES:
  real(r8) :: vegt             ! lsai
  real(r8) :: dewmxi           ! inverse of maximum allowed dew [1/mm]

  if (frac_veg_nosno == 1) then
     if (h2ocan > 0._r8) then
        vegt    = frac_veg_nosno*(elai + esai)
        dewmxi  = 1.0_r8/dewmx
        fwet = ((dewmxi/vegt)*h2ocan)**0.666666666666_r8
        fwet = min (fwet,1.0_r8)   ! Check for maximum limit of fwet
        fdry = (1._r8-fwet)*elai/(elai+esai) ! moved here from below 
     else
        fwet = 0._r8
        fdry = 0._r8 
     end if
  else
     fwet = 0._r8
     fdry = 0._r8
  end if

  return

end subroutine CanopyHydrology_FracWet
