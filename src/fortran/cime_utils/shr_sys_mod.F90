!
! Mocked cime shr_sys_mod
!
! NOTE this abort is simpler than CIME's real abort -- just die
! already and don't require too many other libraries.
!----------------------------------------------------------------

module shr_sys_mod

  use shr_kind_mod, only : shr_kind_in, shr_kind_cx
  implicit none

! PUBLIC: Public interfaces

   private

   ! Imported from shr_abort_mod and republished with renames. Other code that wishes to
   ! use these routines should use these shr_sys names rather than directly using the
   ! routines from shr_abort_abort. (This is for consistency with older code, from when
   ! these routines were defined in shr_sys_mod.)
   public :: shr_sys_abort     ! abort a program

!===============================================================================
contains
!===============================================================================

  !===============================================================================
  subroutine shr_sys_abort(string,rc)
    ! Consistent stopping mechanism

    !----- arguments -----
    character(len=*)    , intent(in), optional :: string  ! error message string
    integer(shr_kind_in), intent(in), optional :: rc      ! error code

    ! Local version of the string.
    ! (Gets a default value if string is not present.)
    character(len=shr_kind_cx) :: local_string
    !-------------------------------------------------------------------------------

    if (present(string)) then
       local_string = trim(string)
    else
       local_string = "Unknown error submitted to shr_abort_abort."
    end if

    ! A compiler's abort method may print a backtrace or do other nice
    ! things, but in fact we can rarely leverage this, because MPI_Abort
    ! usually sends SIGTERM to the process, and we don't catch that signal.
    call abort()

  end subroutine shr_sys_abort

end module shr_sys_mod
