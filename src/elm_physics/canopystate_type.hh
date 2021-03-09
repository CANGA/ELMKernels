/*
DESCRIPTION:
Store and initialize canopy state variables

*/

struct CanopyState_type {

  frac_veg_nosno_patch = 0._r8 tlai_patch = 0._r8 tsai_patch = 0._r8 elai_patch = 0._r8 esai_patch = 0._r8 htop_patch =
      0._r8 hbot_patch = 0._r8 dewmx_patch = 0.1_r8 vegwp_patch[4] = {-2.5e4, -2.5e4, -2.5e4, -2.5e4};
};

subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: p,l,c,g 
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       l = veg_pp%landunit(p)

       this%frac_veg_nosno_patch(p) = 0._r8
       this%tlai_patch(p)       = 0._r8
       this%tsai_patch(p)       = 0._r8
       this%elai_patch(p)       = 0._r8
       this%esai_patch(p)       = 0._r8
       this%htop_patch(p)       = 0._r8
       this%hbot_patch(p)       = 0._r8
       this%dewmx_patch(p)      = 0.1_r8
       this%vegwp_patch(p,:)    = -2.5e4_r8

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%laisun_patch(p) = 0._r8
          this%laisha_patch(p) = 0._r8
       end if

       ! needs to be initialized to spval to avoid problems when averaging for the accum
       ! field
       this%fsun_patch(p) = spval
    end do

    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%alt_col(c)               = 0._r8 !iniitialized to spval for all columns
          this%altmax_col(c)            = 0._r8 !iniitialized to spval for all columns
          this%altmax_lastyear_col(c)   = 0._r8 !iniitialized to spval for all columns
          this%alt_indx_col(c)          = 0     !initiialized to huge  for all columns
          this%altmax_indx_col(c)       = 0     !initiialized to huge  for all columns
          this%altmax_lastyear_indx_col = 0     !initiialized to huge  for all columns
       end if
    end do

  end subroutine InitCold